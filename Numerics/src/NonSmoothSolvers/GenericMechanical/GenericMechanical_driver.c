/* Siconos-Numerics, Copyright INRIA 2005-2010.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include "LA.h"
#include "NumericsOptions.h"
#include "GenericMechanical_Solvers.h"
#include "NonSmoothDrivers.h"
#include "FrictionContact3D_compute_error.h"
#include "FrictionContact3D_unitary_enumerative.h"
#include "GMPReduced.h"
/* #define GENERICMECHANICAL_DEBUG */
/* #define GENERICMECHANICAL_DEBUG2 */
/* #define GENERICMECHANICAL_DEBUG_CMP */
/* #define GENERICMECHANICAL_DEBUG_COMPUTE_ERROR */
/* #define GENERICMECHANICAL_FC3D */
/* #define GMP_WRITE_FAILED_PRB */
/* #define GMP_WRITE_PRB */
int GenericMechanical_compute_error(GenericMechanicalProblem* pGMP, double *reaction , double *velocity, double tol, SolverOptions* options, double * err)
{
  listNumericsProblem * curProblem = pGMP->firstListElem;
  SparseBlockStructuredMatrix* m = pGMP->M->matrix1;
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
    if (storageType == 0)
    {
      rowProdNoDiag(pGMP->size, curSize, posInX, numMat, reaction, curProblem->q, 0);
    }
    else
      rowProdNoDiagSBM(pGMP->size, curSize, currentRowNumber, m, reaction, curProblem->q, 0);
#ifdef GENERICMECHANICAL_DEBUG_COMPUTE_ERROR
    printDenseMatrice("qnodiag", 0, curProblem->q, curSize, 1);
#endif
    /*computation of the velocity of the GMP*/
    memcpy(velocity + posInX, curProblem->q, curSize * sizeof(double));
    /*add the missing product to the velocity: the diagonal one*/

    double * diagBlock = 0;
    if (storageType == 0) /*dense*/
    {
      getDiagonalBlock(numMat, currentRowNumber, posInX, curSize, &bufForLocalProblemDense);
      diagBlock = bufForLocalProblemDense;
    }
    else
    {
      getDiagonalBlock(numMat, currentRowNumber, posInX, curSize, &diagBlock);
    }
#ifdef GENERICMECHANICAL_DEBUG_COMPUTE_ERROR
    printDenseMatrice("diagBlock", 0, diagBlock, curSize, curSize);
    printDenseMatrice("Rlocal", 0, reaction + posInX, curSize, 1);
#endif
    DGEMV(LA_NOTRANS, curSize, curSize, 1.0, diagBlock, curSize, reaction + posInX, 1, 1.0, velocity + posInX, 1);
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
      //  LinearSystemProblem* linearProblem = (LinearSystemProblem*) curProblem->problem;
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
      localError = localError / (1 + DNRM2(curSize , curProblem->q , 1));
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
      FrictionContact3D_unitary_compute_and_add_error(reaction + posInX, velocity + posInX, fcProblem->mu[0], &localError);
      localError = sqrt(localError) / (1 + DNRM2(curSize , curProblem->q , 1));
      if (localError > *err)
        *err = localError ;
#ifdef GENERICMECHANICAL_DEBUG_COMPUTE_ERROR
      printf("GenericMechanical_driver FC3D, Local Velocity v_n=%e v_t1=%e v_t2=%e localerror=%e\n", *(velocity + posInX), *(velocity + posInX + 1), *(velocity + posInX + 2), localError);
#endif
      break;
    }
    default:
      printf("Numerics : genericMechanicalProblem_GS unknown problem type %d.\n", curProblem->type);
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
void genericMechanicalProblem_GS(GenericMechanicalProblem* pGMP, double * reaction, double * velocity, int * info, SolverOptions* options)
{
#ifdef GMP_WRITE_PRB
  FILE * toto1  = fopen("GMP_CURRENT.txt", "w");
  genericMechnical_printInFile(pGMP, toto1);
  fclose(toto1);
#endif

#ifdef GENERICMECHANICAL_DEBUG_CMP
  SScmp++;
#endif
  listNumericsProblem * curProblem = 0;
  int storageType = pGMP->M->storageType;
  SparseBlockStructuredMatrix* m = pGMP->M->matrix1;
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
  //printf("genericMechanicalProblem_GS \n");
  //displayGMP(pGMP);
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
    pPrevReaction = (double *) malloc(genericMechnical_getNbDWork(pGMP, options) * sizeof(double));
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
    int curSize = 0;
#ifdef GENERICMECHANICAL_DEBUG
    printf("GS it %d, initial value:\n", it);
    for (int ii = 0; ii < pGMP->size; ii++)
      printf("R[%d]=%e | V[]=%e \n", ii, reaction[ii], velocity[ii]);
#endif

    while (curProblem)
    {
      //if (currentRowNumber){
      //  posInX = m->blocksize0[currentRowNumber-1];
      //}
      //curSize=m->blocksize0[currentRowNumber] - posInX;
      curSize = curProblem->size;
      /*about the diagonal block:*/
      //diagBlockNumber = getDiagonalBlockPos(m,currentRowNumber);
      //diagBlockNumber = getDiagonalBlockPos(numMat,currentRowNumber,posInX,size);
      double * diagBlock = 0;
      if (storageType == 0) /*dense*/
      {
        getDiagonalBlock(numMat, currentRowNumber, posInX, curSize, &bufForLocalProblemDense);
        diagBlock = bufForLocalProblemDense;
      }
      else
      {
        getDiagonalBlock(numMat, currentRowNumber, posInX, curSize, &diagBlock);

      }

      sol = reaction + posInX;
      w = velocity + posInX;

      switch (curProblem->type)
      {
      case SICONOS_NUMERICS_PROBLEM_EQUALITY:
      {
        /*Mz*/
        LinearSystemProblem* linearProblem = (LinearSystemProblem*) curProblem->problem;
        linearProblem->M->matrix0 = diagBlock;
        /*about q.*/
        memcpy(linearProblem->q, &(pGMP->q[posInX]), curSize * sizeof(double));
        if (storageType == 0)
        {
          rowProdNoDiag(pGMP->size, curSize, posInX, numMat, reaction, linearProblem->q, 0);
        }
        else
          rowProdNoDiagSBM(pGMP->size, curSize, currentRowNumber, m, reaction, linearProblem->q, 0);

        resLocalSolver = LinearSystem_driver(linearProblem, sol, w, 0);

        break;
      }
      case SICONOS_NUMERICS_PROBLEM_LCP:
      {
        /*Mz*/
        LinearComplementarityProblem* lcpProblem = (LinearComplementarityProblem*) curProblem->problem;
        lcpProblem->M->matrix0 = diagBlock;
        /*about q.*/
        memcpy(curProblem->q, &(pGMP->q[posInX]), curSize * sizeof(double));
        if (storageType == 0)
        {
          rowProdNoDiag(pGMP->size, curSize, posInX, numMat, reaction, lcpProblem->q, 0);
        }
        else
          rowProdNoDiagSBM(pGMP->size, curSize, currentRowNumber, m, reaction, lcpProblem->q, 0);
        resLocalSolver = linearComplementarity_driver(lcpProblem, sol, w, options->internalSolvers, 0);

        break;
      }
      case SICONOS_NUMERICS_PROBLEM_FC3D:
      {
        FrictionContactProblem * fcProblem = (FrictionContactProblem *)curProblem->problem;
        fcProblem->M->matrix0 = diagBlock;
        memcpy(curProblem->q, &(pGMP->q[posInX]), curSize * sizeof(double));
        if (storageType == 0)
        {
          rowProdNoDiag(pGMP->size, curSize, posInX, numMat, reaction, fcProblem->q, 0);
        }
        else
          rowProdNoDiagSBM(pGMP->size, curSize, currentRowNumber, m, reaction, fcProblem->q, 0);
        resLocalSolver = frictionContact3D_driver(fcProblem, sol, w, &options->internalSolvers[1], NULL);
        //resLocalSolver=frictionContact3D_unitary_enumerative_solve(fcProblem,sol,&options->internalSolvers[1]);
        break;
      }
      default:
        printf("genericMechanical_GS Numerics : genericMechanicalProblem_GS unknown problem type %d.\n", curProblem->type);
      }
      if (resLocalSolver)
        printf("Local solver FAILED, GS continue\n");
#ifdef GENERICMECHANICAL_DEBUG
      printf("GS it %d, the line number is %d:\n", it, currentRowNumber);
      for (int ii = 0; ii < pGMP->size; ii++)
        printf("R[%d]=%e | V[]=%e \n", ii, reaction[ii], velocity[ii]);
      if (resLocalSolver)
        printf("Numerics:GenericMechanical_drivers Local solver failed\n");
#endif
      posInX += curProblem->size;
      curProblem = curProblem->nextProblem;
      currentRowNumber++;

    }
    /*compute global error.*/

    if (withLS)
    {
      tolViolate = GenericMechanical_compute_error(pGMP, reaction, pBuffVelocity, tol, options, err);
      for (int i = 0; i < pGMP->size; i++)
        pPrevReaction[i] = reaction[i] + (*pCoefLS) * (reaction[i] - pPrevReaction[i]);
      tolViolateLS = GenericMechanical_compute_error(pGMP, pPrevReaction, velocity, tol, options, errLS);
#ifdef GENERICMECHANICAL_DEBUG
      printf("GMD :noscale error=%e error LS=%e\n", *err, *errLS);
      printf("GMD :scale ceoef=%e\n", *pCoefLS);
#endif
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
      tolViolate = GenericMechanical_compute_error(pGMP, reaction, velocity, tol, options, err);
    }
    //tolViolate=GenericMechanical_compute_error(pGMP,reaction,velocity,tol,options,&err);
    /*next GS it*/
    it++;
    //printf("---GenericalMechanical_drivers,  IT=%d, err=%e.\n",it,err);
  }
  /*
  FILE * titi  = fopen("GMP_FAILED_scilab.txt", "w");
  FILE * tata  = fopen("SBM.txt", "w");
  printf("GMP_drivers, print file SBM\n");
  printInFileSBM(pGMP->M->matrix1,tata);
  fclose(tata);
  printInFileSBMForScilab(pGMP->M->matrix1,titi);
  fclose(titi);
  */
  options->iparam[3] = it;
#ifdef GMP_WRITE_FAILED_PRB
  FILE * toto  = fopen("GMP_NOT_FAILED.txt", "w");
  genericMechnical_printInFile(pGMP, toto);
  fclose(toto);
#endif
  if (tolViolate)
  {
    printf("---GenericalMechanical_drivers, FAILED***************************************\n");
#ifdef GMP_WRITE_FAILED_PRB
    FILE * toto  = fopen("GMP_FAILED.txt", "w");
    genericMechnical_printInFile(pGMP, toto);
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
int genericMechanical_driver(GenericMechanicalProblem* problem, double *reaction , double *velocity, SolverOptions* options)
{
  // if (options == NULL )
  //  numericsError("FrictionContact3D_driver", "null input for solver options");

  /* If the options for solver have not been set, read default values in .opt file */

  int info = 0;
#ifdef GENERICMECHANICAL_DEBUG
  display(problem->M);
  displayGMP(problem);
#endif
  if (!options->iparam[2])
  {
    genericMechanicalProblem_GS(problem, reaction, velocity, &info, options);
  }
  else if (options->iparam[2] == 1)
  {
    GMPReducedSolve(problem, reaction, velocity, &info, options);
  }
  else if (options->iparam[2] == 2)
  {
    GMPReducedEqualitySolve(problem, reaction, velocity, &info, options);
  }
  else if (options->iparam[2] == 3)
  {
    GMPasMLCP(problem, reaction, velocity, &info, options);
  }
  else
  {
    printf("genericMechanical_driver error, options->iparam[2] wrong value. (0, 1 or 2).\n");
  }

  return info;

}




void genericMechnicalProblem_setDefaultSolverOptions(SolverOptions* options, int id)
{
  options->iSize = 5;
  options->dSize = 5;
  options->numberOfInternalSolvers = 2;
  options->dWork = NULL;
  options->iWork = NULL;
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
  options->dparam[2] = 1e-7;
  options->dparam[3] = 1e-7;
  options->internalSolvers = (SolverOptions *)malloc(2 * sizeof(SolverOptions));;
  linearComplementarity_setDefaultSolverOptions(0, options->internalSolvers, SICONOS_LCP_LEMKE);
  switch (id)
  {
  case SICONOS_FRICTION_3D_QUARTIC:
  case SICONOS_FRICTION_3D_QUARTIC_NU:
    frictionContact3D_unitary_enumerative_setDefaultSolverOptions(&options->internalSolvers[1]);
    break;
  case SICONOS_FRICTION_3D_AlartCurnierNewton:
  case SICONOS_FRICTION_3D_DampedAlartCurnierNewton:
    frictionContact3D_AlartCurnierNewton_setDefaultSolverOptions(&options->internalSolvers[1]);
    break;
  default:
    printf("FC3D_solverId unknown :%d\n", id);
  }
  //frictionContact3D_AlartCurnierNewton_setDefaultSolverOptions(&options->internalSolvers[1]);
}

/*Alloc memory iff options->iWork options->dWork and are  null.
 Return 0 if the memory is not allocated. else return 1.*/
int genericMechnical_alloc_working_memory(GenericMechanicalProblem* pGMP, SolverOptions* options)
{
  if (options->dWork)
    return 0;

  options->dWork = (double *) malloc(genericMechnical_getNbDWork(pGMP, options) * sizeof(double));
  return 1;

}
/*free the Work memory, and set pointer to zero.*/
void genericMechnical_free_working_memory(GenericMechanicalProblem* pGMP, SolverOptions* options)
{
  if (options->dWork)
    free(options->dWork);
  options->dWork = NULL;
}
int genericMechnical_getNbDWork(GenericMechanicalProblem* pGMP, SolverOptions* options)
{
  return 2 * pGMP->size;
}
