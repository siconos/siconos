/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <math.h>

#include "GenericMechanical_Solvers.h"
#include "GenericMechanical_cst.h"

#include "NonSmoothDrivers.h"

#include "GMPReduced.h"
#include "SiconosBlas.h"
#include "numerics_verbose.h"

#include "LCP_Solvers.h"

#include "FrictionContactProblem.h"
#include "Friction_cst.h"
#include "fc3d_compute_error.h"
#include "fc3d_unitary_enumerative.h"
#include "fc3d_onecontact_nonsmooth_Newton_solvers.h"

#include "Relay_Solvers.h"
#include "RelayProblem.h"
#include "relay_cst.h"

#include "NumericsMatrix.h"

#include "numerics_verbose.h"
/* #define GENERICMECHANICAL_DEBUG  */
/* #define GENERICMECHANICAL_DEBUG2  */
/* #define GENERICMECHANICAL_DEBUG_CMP */
/* #define GENERICMECHANICAL_DEBUG_COMPUTE_ERROR */
/* #define GENERICMECHANICAL_FC3D */
/* #define GMP_WRITE_FAILED_PRB */
/* #define GMP_WRITE_PRB */


const char* const   SICONOS_GENERIC_MECHANICAL_NSGS_STR = "GMP_NSGS";


//#define DEBUG_MESSAGES
//#define DEBUG_STDOUT
#include "debug.h"
#include "LinearComplementarityProblem.h"
#include "lcp_cst.h"

int gmp_compute_error(GenericMechanicalProblem* pGMP, double *reaction , double *velocity, double tol, SolverOptions* options, double * err)
{
  listNumericsProblem * curProblem = pGMP->firstListElem;
  int storageType = pGMP->M->storageType;
  NumericsMatrix* numMat = pGMP->M;
  int currentRowNumber = 0;
  int ii;
  int  posInX = 0;
  int curSize = 0;
  *err = 0.0;
  double localError = 0;
  double * bufForLocalProblemDense = (storageType == 0) ? (double*) malloc(pGMP->maxLocalSize * pGMP->maxLocalSize * sizeof(double)) : 0;

#ifdef GENERICMECHANICAL_DEBUG_COMPUTE_ERROR
  printf("GenericMechanical compute_error BEGIN:\n");
#endif
  /*update localProblem->q and compute V = M*R+Q of the GMP */
  posInX = 0;
  while (curProblem)
  {
    curSize = curProblem->size;

    /*localproblem->q <-- GMP->q */
    memcpy(curProblem->q, &(pGMP->q[posInX]), curSize * sizeof(double));
#ifdef GENERICMECHANICAL_DEBUG_COMPUTE_ERROR
    printDenseMatrice("q", 0, curProblem->q, curSize, 1);
    printDenseMatrice("reaction", 0, reaction, pGMP->size, 1);
#endif
    /*computation of the localproblem->q*/
    NM_row_prod_no_diag(pGMP->size, curSize, currentRowNumber, posInX, numMat, reaction, curProblem->q, NULL, 0);
#ifdef GENERICMECHANICAL_DEBUG_COMPUTE_ERROR
    printDenseMatrice("qnodiag", 0, curProblem->q, curSize, 1);
#endif
    /*computation of the velocity of the GMP*/
    memcpy(velocity + posInX, curProblem->q, curSize * sizeof(double));
    /*add the missing product to the velocity: the diagonal one*/

    double * diagBlock = 0;
    if (storageType == 0) /*dense*/
    {
      NM_extract_diag_block(numMat, currentRowNumber, posInX, curSize, &bufForLocalProblemDense);
      diagBlock = bufForLocalProblemDense;
    }
    else
    {
      NM_extract_diag_block(numMat, currentRowNumber, posInX, curSize, &diagBlock);
    }
#ifdef GENERICMECHANICAL_DEBUG_COMPUTE_ERROR
    printDenseMatrice("diagBlock", 0, diagBlock, curSize, curSize);
    printDenseMatrice("Rlocal", 0, reaction + posInX, curSize, 1);
#endif
    cblas_dgemv(CblasColMajor,CblasNoTrans, curSize, curSize, 1.0, diagBlock, curSize, reaction + posInX, 1, 1.0, velocity + posInX, 1);
#ifdef GENERICMECHANICAL_DEBUG_COMPUTE_ERROR
    printDenseMatrice("velocity", 0, velocity + posInX, curSize, 1);
#endif
    /*next*/
    posInX += curProblem->size;
    curProblem = curProblem->nextProblem;
    currentRowNumber++;
  }


  /*For each sub-problem, call the corresponding function computing the error.*/
  posInX = 0;
  currentRowNumber = 0;
  curProblem = pGMP->firstListElem;
  while (curProblem)
  {
    curSize = curProblem->size;
    double * Vl = velocity + posInX;
    double * Rl = reaction + posInX;
    for (ii = 0; ii < curSize; ii++)
      if (isnan(Vl[ii]) || isnan(Rl[ii]))
      {
        *err = 10;
        return 1;
      }
    switch (curProblem->type)
    {
    case SICONOS_NUMERICS_PROBLEM_EQUALITY:
    {
      double * w = velocity + posInX;
      localError = 0.;
      for (ii = 0; ii < curSize; ii++)
      {
        if (fabs(w[ii]) > *err)
          *err = fabs(w[ii]);
        if (fabs(w[ii]) > localError)
          localError = *err;

      }
#ifdef GENERICMECHANICAL_DEBUG_COMPUTE_ERROR
      printf("GenericMechanical_driver, localerror of linearSystem: %e\n", localError);
#endif
      break;
    }
    case SICONOS_NUMERICS_PROBLEM_LCP:
    {
      //LinearComplementarityProblem* lcpproblem = (LinearComplementarityProblem*) curProblem->problem;

      lcp_compute_error_only(curSize, reaction + posInX, velocity + posInX, &localError);
      localError = localError / (1 + cblas_dnrm2(curSize , curProblem->q , 1));
      if (localError > *err)
        *err = localError ;
#ifdef GENERICMECHANICAL_DEBUG_COMPUTE_ERROR
      printf("GenericMechanical_driver, localerror of lcp: %e\n", localError);
#endif
      break;
    }
    case SICONOS_NUMERICS_PROBLEM_RELAY:
    {
      relay_compute_error((RelayProblem*) curProblem->problem,
                          reaction + posInX, velocity + posInX,
                          options->dparam[0], &localError);

      localError = localError / (1 + cblas_dnrm2(curSize , curProblem->q , 1));
      if (localError > *err)
        *err = localError ;
#ifdef GENERICMECHANICAL_DEBUG_COMPUTE_ERROR
      printf("GenericMechanical_driver, localerror of lcp: %e\n", localError);
#endif
      break;
    }
    case SICONOS_NUMERICS_PROBLEM_FC3D:
    {
      FrictionContactProblem * fcProblem = (FrictionContactProblem *)curProblem->problem;
      localError = 0.;
      double worktmp[3];
      fc3d_unitary_compute_and_add_error(reaction + posInX, velocity + posInX, fcProblem->mu[0], &localError, worktmp);
      localError = sqrt(localError) / (1 + cblas_dnrm2(curSize , curProblem->q , 1));
      if (localError > *err)
        *err = localError ;
#ifdef GENERICMECHANICAL_DEBUG_COMPUTE_ERROR
      printf("GenericMechanical_driver FC3D, Local Velocity v_n=%e v_t1=%e v_t2=%e localerror=%e\n", *(velocity + posInX), *(velocity + posInX + 1), *(velocity + posInX + 2), localError);
#endif
      break;
    }
    default:
      printf("Numerics : gmp_gauss_seidel unknown problem type %d.\n", curProblem->type);
    }
    /*next*/
    posInX += curProblem->size;
    curProblem = curProblem->nextProblem;
    currentRowNumber++;
  }
#ifdef GENERICMECHANICAL_DEBUG_COMPUTE_ERROR
  if (*err > tol)
    printf("GenericMechanical_driver compute_error END:, err>tol: error : %e\n", *err);
  else
    printf("GenericMechanical_driver compute_error END:, err<tol: error : %e\n", *err);
#endif
  if (storageType == 0)
    free(bufForLocalProblemDense);

  if (*err > tol)
    return 1;
  else
    return 0;
}
#ifdef GENERICMECHANICAL_DEBUG_CMP
static int SScmp = 0;
static int SScmpTotal = 0;
#endif
//#define GMP_WRITE_PRB
//static double sCoefLS=1.0;
void gmp_gauss_seidel(GenericMechanicalProblem* pGMP, double * reaction, double * velocity, int * info,
				 SolverOptions* options)
{
#ifdef GMP_WRITE_PRB
  FILE * toto1  = fopen("GMP_CURRENT.txt", "w");
  genericMechanicalProblem_printInFile(pGMP, toto1);
  fclose(toto1);
#endif

#ifdef GENERICMECHANICAL_DEBUG_CMP
  SScmp++;
#endif
  listNumericsProblem * curProblem = 0;
  int storageType = pGMP->M->storageType;
  NumericsMatrix* numMat = pGMP->M;
  int iterMax = options->iparam[0];
  int it = 0;
  int currentRowNumber = 0;
  //  int diagBlockNumber =0;
  double tol = options->dparam[0];
  double * err = &(options->dparam[2]);
  double * errLS = &(options->dparam[3]);
  int tolViolate = 1;
  int tolViolateLS = 1;
  double * sol = 0;
  double * w = 0;
  int resLocalSolver = 0;
  int local_solver_error_occurred = 0;
  //printf("gmp_gauss_seidel \n");
  //genericMechanicalProblem_display(pGMP);
  double * pPrevReaction = NULL;
  double * pBuffVelocity = NULL;
  int withLS = options->iparam[1];
  double * pCoefLS = &(options->dparam[1]);
  double * bufForLocalProblemDense = (storageType == 0) ? (double*) malloc(pGMP->maxLocalSize * pGMP->maxLocalSize * sizeof(double)) : 0;

  if (options->dWork)
  {
    pPrevReaction = options->dWork;
  }
  else
  {
    pPrevReaction = (double *) malloc(gmp_get_nb_dwork(pGMP, options) * sizeof(double));
  }
  pBuffVelocity = pPrevReaction + pGMP->size;
  while (it < iterMax && tolViolate)
  {
#ifdef GENERICMECHANICAL_DEBUG_CMP
    SScmpTotal++;
#endif
    memcpy(pPrevReaction, reaction, pGMP->size * sizeof(double));
    currentRowNumber = 0;
    curProblem =  pGMP->firstListElem;
    int  posInX = 0;
    size_t curSize = 0;

    DEBUG_PRINTF("GS it %d, initial value:\n", it);
    DEBUG_EXPR(
      for (int ii = 0; ii < pGMP->size; ii++)
        printf("R[%i]=%e | V[%i]=%e \n", ii, reaction[ii], ii, velocity[ii]);
      );

    while (curProblem)
    {
      //if (currentRowNumber){
      //  posInX = m->blocksize0[currentRowNumber-1];
      //}
      //curSize=m->blocksize0[currentRowNumber] - posInX;
      curSize = curProblem->size;
      curProblem->error = 0;
      /*about the diagonal block:*/
      //diagBlockNumber = NM_extract_diag_blockPos(m,currentRowNumber);
      //diagBlockNumber = NM_extract_diag_blockPos(numMat,currentRowNumber,posInX,size);
      double * diagBlock = 0;
      if (storageType == 0) /*dense*/
      {
        NM_extract_diag_block(numMat, currentRowNumber, posInX, curSize, &bufForLocalProblemDense);
        diagBlock = bufForLocalProblemDense;
      }
      else
      {
        NM_extract_diag_block(numMat, currentRowNumber, posInX, curSize, &diagBlock);

      }

      sol = reaction + posInX;
      w = velocity + posInX;

      switch (curProblem->type)
      {
      case SICONOS_NUMERICS_PROBLEM_EQUALITY:
      {
        NumericsMatrix M;
        NM_fill(&M, NM_DENSE, curSize, curSize, diagBlock);

        memcpy(curProblem->q, &(pGMP->q[posInX]), curSize * sizeof(double));
        NM_row_prod_no_diag(pGMP->size, curSize, currentRowNumber, posInX, numMat, reaction, curProblem->q, NULL, 0);
        for (size_t i = 0; i < curSize; ++i) sol[i] = -curProblem->q[i];

        resLocalSolver = NM_gesv(&M, sol, true);

        M.matrix0 = NULL;
        NM_clear(&M);
        break;
      }
      case SICONOS_NUMERICS_PROBLEM_LCP:
      {
        /*Mz*/
        LinearComplementarityProblem* lcpProblem = (LinearComplementarityProblem*) curProblem->problem;
        lcpProblem->M->matrix0 = diagBlock;
        /*about q.*/
        memcpy(curProblem->q, &(pGMP->q[posInX]), curSize * sizeof(double));
        NM_row_prod_no_diag(pGMP->size, curSize, currentRowNumber, posInX, numMat, reaction, lcpProblem->q, NULL, 0);
        resLocalSolver = linearComplementarity_driver(lcpProblem, sol, w, options->internalSolvers);
        break;
      }
      case SICONOS_NUMERICS_PROBLEM_RELAY:
      {
        /*Mz*/
        RelayProblem* relayProblem = (RelayProblem*) curProblem->problem;
        relayProblem->M->matrix0 = diagBlock;
        /*about q.*/
        memcpy(curProblem->q, &(pGMP->q[posInX]), curSize * sizeof(double));
        NM_row_prod_no_diag(pGMP->size, curSize, currentRowNumber, posInX, numMat, reaction, relayProblem->q, NULL, 0);
        resLocalSolver = relay_driver(relayProblem, sol, w, &options->internalSolvers[2]);
        break;
      }
      case SICONOS_NUMERICS_PROBLEM_FC3D:
      {
        FrictionContactProblem * fcProblem = (FrictionContactProblem *)curProblem->problem;
        assert(fcProblem);
        assert(fcProblem->M);
        assert(fcProblem->q);
        fcProblem->M->matrix0 = diagBlock;
        memcpy(curProblem->q, &(pGMP->q[posInX]), curSize * sizeof(double));

        DEBUG_EXPR_WE(for (int i =0 ; i < 3; i++) printf("curProblem->q[%i]= %12.8e,\t fcProblem->q[%i]= %12.8e,\n",i,curProblem->q[i],i,fcProblem->q[i]););

        NM_row_prod_no_diag(pGMP->size, curSize, currentRowNumber, posInX, numMat, reaction, fcProblem->q, NULL, 0);

        DEBUG_EXPR_WE(for (int i =0 ; i < 3; i++)  printf("reaction[%i]= %12.8e,\t fcProblem->q[%i]= %12.8e,\n",i,reaction[i],i,fcProblem->q[i]););

        /* We call the generic driver (rather than the specific) since we may choose between various local solvers */
        resLocalSolver = fc3d_driver(fcProblem, sol, w, &options->internalSolvers[1]);
        //resLocalSolver=fc3d_unitary_enumerative_solve(fcProblem,sol,&options->internalSolvers[1]);
        break;
      }
      default:
        printf("genericMechanical_GS Numerics : gmp_gauss_seidel unknown problem type %d.\n", curProblem->type);
      }
      if (resLocalSolver) {
        curProblem->error = 1;
        local_solver_error_occurred = 1;
      }

      DEBUG_PRINTF("GS it %d, the line number is %d:\n", it, currentRowNumber);
      /* DEBUG_EXPR(for (int ii = 0; ii < pGMP->size; ii++) */
      /*              printf("R[%d]=%e | V[]=%e \n", ii, reaction[ii], velocity[ii]); */
      /*            if (resLocalSolver) */
      /*              printf("Numerics:GenericMechanical_drivers Local solver failed\n"); */
      /*   ); */
      posInX += curProblem->size;
      curProblem = curProblem->nextProblem;
      currentRowNumber++;

    }
    /*compute global error.*/

    if (withLS)
    {
      tolViolate = gmp_compute_error(pGMP, reaction, pBuffVelocity, tol, options, err);
      for (int i = 0; i < pGMP->size; i++)
        pPrevReaction[i] = reaction[i] + (*pCoefLS) * (reaction[i] - pPrevReaction[i]);
      tolViolateLS = gmp_compute_error(pGMP, pPrevReaction, velocity, tol, options, errLS);

      DEBUG_PRINTF("GMP :noscale error=%e error LS=%e\n", *err, *errLS);
      DEBUG_PRINTF("GMP :scale coeff=%e\n", *pCoefLS);

      if (*errLS < *err)
      {
        if ((*pCoefLS) < 10.0)
          (*pCoefLS) = 1.0 + (*pCoefLS);
        memcpy(reaction, pPrevReaction, pGMP->size * sizeof(double));
        tolViolate = tolViolateLS;
        *err = *errLS;
      }
      else
      {
        *pCoefLS = 1.0;
        memcpy(velocity, pBuffVelocity, pGMP->size * sizeof(double));
      }



    }
    else
    {
      tolViolate = gmp_compute_error(pGMP, reaction, velocity, tol, options, err);
    }
    if (verbose > 0)
      printf("--------------- GMP - GS - Iteration %i Residual = %14.7e <= %7.3e\n", it, *err, options->dparam[0]);

    //tolViolate=gmp_compute_error(pGMP,reaction,velocity,tol,options,&err);
    /*next GS it*/
    it++;
    //printf("---GenericalMechanical_drivers,  IT=%d, err=%e.\n",it,err);
  }
  /*
  FILE * titi  = fopen("GMP_FAILED_scilab.txt", "w");
  FILE * tata  = fopen("SBM.txt", "w");
  printf("GMP_drivers, print file SBM\n");
  SBM_write_in_file(pGMP->M->matrix1,tata);
  fclose(tata);
  SBM_write_in_fileForScilab(pGMP->M->matrix1,titi);
  fclose(titi);
  */
  options->iparam[3] = it;
#ifdef GMP_WRITE_FAILED_PRB
  FILE * toto  = fopen("GMP_NOT_FAILED.txt", "w");
  genericMechanicalProblem_printInFile(pGMP, toto);
  fclose(toto);
#endif
  if (tolViolate)
  {
    if (verbose > 0)
      printf("gmp_gauss_seidel failed with Iteration %i Residual = %14.7e <= %7.3e\n", it, *err, options->dparam[0]);

#ifdef GMP_WRITE_FAILED_PRB
    FILE * toto  = fopen("GMP_FAILED.txt", "w");
    genericMechanicalProblem_printInFile(pGMP, toto);
    fclose(toto);
#endif

    //    if (storageType==0)
    //      free(bufForLocalProblemDense);
    //    exit(0);
  }
  else
  {
    ;
#ifdef GENERICMECHANICAL_DEBUG_CMP
    printf("---GenericalMechanical_drivers, CV %d at it=%d, itTotal=%d.\n", SScmp, it, SScmpTotal);
#endif
  }

  if (local_solver_error_occurred) {
    currentRowNumber = 0;
    curProblem =  pGMP->firstListElem;
    while (curProblem) {
      if (curProblem->error && verbose)
        printf("genericMechanical_GS Numerics : Local solver FAILED row %d of type %s\n",
               currentRowNumber, ns_problem_id_to_name(curProblem->type));
      curProblem = curProblem->nextProblem;
      currentRowNumber++;
    }
  }

  //printf("---GenericalMechanical_drivers,  IT=%d, err=%e.\n",it,*err);
  if (! options->dWork)
    free(pPrevReaction);
  *info = tolViolate;
  if (storageType == 0)
    free(bufForLocalProblemDense);
}

/*
 * options->iparam[2] == 0 then GS on all blocks
 * options->iparam[2] == 1 The equalities are substituated
 * options->iparam[2] == 2 Equalities are assemblated in one block
 * options->iparam[2] == 3 Try to solve like a MLCP (==> No FC3d)
 */
int gmp_driver(GenericMechanicalProblem* problem, double *reaction , double *velocity,
                             SolverOptions* options)
{
  // if (options == NULL )
  //  numerics_error("fc3d_driver", "null input for solver options");

  /* If the options for solver have not been set, read default values in .opt file */

  int info = 0;
  DEBUG_EXPR(
    NM_display(problem->M);
    genericMechanicalProblem_display(problem);
    );
  if (!options->iparam[2])
  {
    gmp_gauss_seidel(problem, reaction, velocity, &info, options);
  }
  else if (options->iparam[2] == 1)
  {
    gmp_reduced_solve(problem, reaction, velocity, &info, options);
  }
  else if (options->iparam[2] == 2)
  {
    gmp_reduced_equality_solve(problem, reaction, velocity, &info, options);
  }
  else if (options->iparam[2] == 3)
  {
    gmp_as_mlcp(problem, reaction, velocity, &info, options);
  }
  else
  {
    printf("gmp_driver error, options->iparam[2] wrong value. (0, 1 or 2).\n");
  }

  return info;

}




void gmp_setDefaultSolverOptions(SolverOptions* options, int id)
{
  options->solverId = SICONOS_GENERIC_MECHANICAL_NSGS;
  options->iSize = 15;
  options->dSize = 15;
  options->numberOfInternalSolvers = 3;
  options->dWork = NULL;
  solver_options_nullify(options);
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  options->iparam[0] = 10000;
  /*with Line search 1 without 0.*/
  options->iparam[1] = 0;

  options->dparam[0] = 1e-4;
  /*Useful parameter for LS*/
  options->dparam[1] = 1.0;
  options->dparam[2] = 1e-7;
  options->dparam[3] = 1e-7;

  options->internalSolvers = (SolverOptions *)malloc(3 * sizeof(SolverOptions));;

  linearComplementarity_setDefaultSolverOptions(0, options->internalSolvers, SICONOS_LCP_LEMKE);
  relay_setDefaultSolverOptions(0, &options->internalSolvers[2], SICONOS_RELAY_LEMKE);

  switch (id)
  {
  case SICONOS_FRICTION_3D_ONECONTACT_QUARTIC:
  case SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU:
    fc3d_unitary_enumerative_setDefaultSolverOptions(&options->internalSolvers[1]);
    break;
  case SICONOS_FRICTION_3D_ONECONTACT_NSN:
    fc3d_onecontact_nonsmooth_Newton_setDefaultSolverOptions(&options->internalSolvers[1]);
    (&options->internalSolvers[1])->solverId=SICONOS_FRICTION_3D_ONECONTACT_NSN;
    //(&options->internalSolvers[1])->iparam[10]=1; /* VA 26/11/2015 For robustness reasons on mechanisms, we choose the JeanMoreau formulation */
    break;
  case SICONOS_FRICTION_3D_ONECONTACT_NSN_GP:
  case SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID:
    fc3d_onecontact_nonsmooth_Newton_setDefaultSolverOptions(&options->internalSolvers[1]);
    //(&options->internalSolvers[1])->iparam[10]=1; /* VA 26/11/2015 For robustness reasons on mechanisms, we choose the JeanMoreau formulation */
    break;
  case SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration:
    (&options->internalSolvers[1])->solverId=SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration;
    (&options->internalSolvers[1])->numberOfInternalSolvers = 0;
    (&options->internalSolvers[1])->isSet = 1;
    (&options->internalSolvers[1])->filterOn = 1;
    (&options->internalSolvers[1])->iSize = 5;
    (&options->internalSolvers[1])->dSize = 5;
    (&options->internalSolvers[1])->iparam = (int *)malloc(options->iSize * sizeof(int));
    (&options->internalSolvers[1])->dparam = (double *)malloc(options->dSize * sizeof(double));
    solver_options_nullify((&options->internalSolvers[1]));
    for (int i = 0; i < 5; i++)
    {
      (&options->internalSolvers[1])->iparam[i] = 0;
      (&options->internalSolvers[1])->dparam[i] = 0.0;
    }
    (&options->internalSolvers[1])->iparam[0] = 100;
    (&options->internalSolvers[1])->dparam[0] = 1e-12;
    break;
  default:
    printf("gmp_setDefaultSolverOptions : fc3d_solverId unknown :%d\n", id);
  }


  //fc3d_AlartCurnierNewton_setDefaultSolverOptions(&options->internalSolvers[1]);
}

/*Alloc memory iff options->iWork options->dWork and are  null.
 Return 0 if the memory is not allocated. else return 1.*/
int gmp_working_memory_alloc(GenericMechanicalProblem* pGMP, SolverOptions* options)
{
  if (options->dWork)
    return 0;

  options->dWork = (double *) malloc(gmp_get_nb_dwork(pGMP, options) * sizeof(double));
  return 1;

}
/*free the Work memory, and set pointer to zero.*/
void gmp_working_memory_free(GenericMechanicalProblem* pGMP, SolverOptions* options)
{
  if (options->dWork)
    free(options->dWork);
  options->dWork = NULL;
}
int gmp_get_nb_dwork(GenericMechanicalProblem* pGMP, SolverOptions* options)
{
  return 2 * pGMP->size;
}
