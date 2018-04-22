/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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


#include "fc3d_Solvers.h"
#include "projectionOnCone.h"
#include "projectionOnCylinder.h"
#include "fc3d_compute_error.h"

#include "SiconosBlas.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "sanitizer.h"
#include "numerics_verbose.h"

/* #define DEBUG_NOCOLOR */
/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "debug.h"

/* Static variables */

/* The global problem of size n= 3*nc, nc being the number of contacts, is locally saved in MGlobal and qGlobal */
/* mu corresponds to the vector of friction coefficients */
/* note that either MGlobal or MBGlobal is used, depending on the chosen storage */
/* static int n=0; */
/* static const NumericsMatrix* MGlobal = NULL; */
/* static const double* qGlobal = NULL; */
/* static const double* mu = NULL; */

/* Local problem operators */
/* static const int nLocal = 3; */
/* static double* MLocal; */
/* static int isMAllocatedIn = 0; /\* True if a malloc is done for MLocal, else false *\/ */
/* static double qLocal[3]; */
/* static double mu_i = 0.0; */


void fc3d_projection_initialize(FrictionContactProblem * problem, FrictionContactProblem * localproblem)
{

}

void fc3d_projection_update(int contact, FrictionContactProblem* problem, FrictionContactProblem* localproblem, double* reaction, SolverOptions* options)
{
  /* Build a local problem for a specific contact
     reaction corresponds to the global vector (size n) of the global problem.
  */

  /* Call the update function which depends on the storage for MGlobal/MBGlobal */
  /* Build a local problem for a specific contact
   reaction corresponds to the global vector (size n) of the global problem.
  */

  /* The part of MGlobal which corresponds to the current block is copied into MLocal */
  fc3d_local_problem_fill_M(problem, localproblem, contact);

  /****  Computation of qLocal = qBlock + sum over a row of blocks in MGlobal of the products MLocal.reactionBlock,
     excluding the block corresponding to the current contact. ****/
  fc3d_local_problem_compute_q(problem, localproblem, reaction, contact);

  /* Friction coefficient for current block*/
  localproblem->mu[0] = problem->mu[contact];

}

void fc3d_projectionWithDiagonalization_update(int contact, FrictionContactProblem* problem, FrictionContactProblem* localproblem,  double* reaction, SolverOptions* options)
{
  /* Build a local problem for a specific contact
     reaction corresponds to the global vector (size n) of the global problem.
  */

  /* Call the update function which depends on the storage for MGlobal/MBGlobal */
  /* Build a local problem for a specific contact
   reaction corresponds to the global vector (size n) of the global problem.
  */

  /* The part of MGlobal which corresponds to the current block is copied into MLocal */
  fc3d_local_problem_fill_M(problem, localproblem, contact);

  /****  Computation of qLocal = qBlock + sum over a row of blocks in MGlobal of the products MLocal.reactionBlock,
     excluding the block corresponding to the current contact. ****/

  NumericsMatrix * MGlobal = problem->M;
  double * MLocal =  localproblem->M->matrix0;


  double *qLocal = localproblem->q;
  double * qGlobal = problem->q;
  int n = 3 * problem->numberOfContacts;


  int in = 3 * contact, it = in + 1, is = it + 1;
  /* reaction current block set to zero, to exclude current contact block */
  /*   double rin= reaction[in] ; double rit= reaction[it] ; double ris= reaction[is] ;  */
  /* qLocal computation*/
  qLocal[0] = qGlobal[in];
  qLocal[1] = qGlobal[it];
  qLocal[2] = qGlobal[is];

  if (MGlobal->storageType == NM_DENSE)
  {
    double * MM = MGlobal->matrix0;
    int incx = n, incy = 1;
    qLocal[0] += cblas_ddot(n , &MM[in] , incx , reaction , incy);
    qLocal[1] += cblas_ddot(n , &MM[it] , incx , reaction , incy);
    qLocal[2] += cblas_ddot(n , &MM[is] , incx , reaction , incy);
    // Substract diagonal term
    qLocal[0] -= MM[in + n * in] * reaction[in];
    qLocal[1] -= MM[it + n * it] * reaction[it];
    qLocal[2] -= MM[is + n * is] * reaction[is];
  }
  else if (MGlobal->storageType == NM_SPARSE_BLOCK)
  {
    /* qLocal += rowMB * reaction
       with rowMB the row of blocks of MGlobal which corresponds to the current contact
    */
    SBM_row_prod(n, 3, contact, MGlobal->matrix1, reaction, qLocal, 0);
    // Substract diagonal term
    qLocal[0] -= MLocal[0] * reaction[in];
    qLocal[1] -= MLocal[4] * reaction[it];
    qLocal[2] -= MLocal[8] * reaction[is];

  }
  else
  {
    fprintf(stderr, "fc3d_projectionWithDiagonalization_update :: Unsupported matrix storage)");
    exit(EXIT_FAILURE);
  }
  /*   reaction[in] = rin; reaction[it] = rit; reaction[is] = ris; */

  /* Friction coefficient for current block*/
  localproblem->mu[0] = problem->mu[contact];
}


void fc3d_projection_initialize_with_regularization(FrictionContactProblem * problem, FrictionContactProblem * localproblem)
{
  if (!localproblem->M->matrix0)
    localproblem->M->matrix0 = (double*)calloc(9, sizeof(double));
}

void fc3d_projection_update_with_regularization(int contact, FrictionContactProblem * problem, FrictionContactProblem * localproblem, double* reaction, SolverOptions* options)
{


  /* Build a local problem for a specific contact
     reaction corresponds to the global vector (size n) of the global problem.
  */

  /* Call the update function which depends on the storage for MGlobal/MBGlobal */
  /* Build a local problem for a specific contact
   reaction corresponds to the global vector (size n) of the global problem.
  */

  /* The part of MGlobal which corresponds to the current block is copied into MLocal */

  NM_copy_diag_block3(problem->M, contact, &localproblem->M->matrix0);

  /****  Computation of qLocal = qBlock + sum over a row of blocks in MGlobal of the products MLocal.reactionBlock,
     excluding the block corresponding to the current contact. ****/
  fc3d_local_problem_compute_q(problem, localproblem, reaction, contact);

  double rho = options->dparam[3];
  for (int i = 0 ; i < 3 ; i++) localproblem->M->matrix0[i + 3 * i] += rho ;

  double *qLocal = localproblem->q;
  int in = 3 * contact, it = in + 1, is = it + 1;

  /* qLocal computation*/
  qLocal[0] -= rho * reaction[in];
  qLocal[1] -= rho * reaction[it];
  qLocal[2] -= rho * reaction[is];

  /* Friction coefficient for current block*/
  localproblem->mu[0] = problem->mu[contact];


}

int fc3d_projectionWithDiagonalization_solve(FrictionContactProblem* localproblem, double* reaction, SolverOptions * options)
{



  /* Current block position */

  /* Builds local problem for the current contact */
  /*  fc3d_projection_update(contact, reaction); */
  /*  fc3d_projectionWithDiagonalization_update(contact, reaction);  */


  double * MLocal = localproblem->M->matrix0;
  double * qLocal = localproblem->q;
  double mu_i = localproblem->mu[0];
  int nLocal = 3;

  double mrn, num, mu2 = mu_i * mu_i;


  /* projection */
  if (qLocal[0] > 0.)
  {
    reaction[0] = 0.;
    reaction[1] = 0.;
    reaction[2] = 0.;
  }
  else
  {
    if (MLocal[0] < DBL_EPSILON || MLocal[nLocal + 1] < DBL_EPSILON || MLocal[2 * nLocal + 2] < DBL_EPSILON)
    {
      fprintf(stderr, "fc3d_projection error: null term on MLocal diagonal.\n");
      exit(EXIT_FAILURE);
    }

    reaction[0] = -qLocal[0] / MLocal[0];
    reaction[1] = -qLocal[1] / MLocal[nLocal + 1];
    reaction[2] = -qLocal[2] / MLocal[2 * nLocal + 2];

    mrn = reaction[1] * reaction[1] + reaction[2] * reaction[2];

    if (mrn > mu2 * reaction[0]*reaction[0])
    {
      num = mu_i * reaction[0] / sqrt(mrn);
      reaction[1] = reaction[1] * num;
      reaction[2] = reaction[2] * num;
    }
  }
  return 0;
}

void fc3d_projectionOnConeWithLocalIteration_initialize(FrictionContactProblem * problem, FrictionContactProblem * localproblem, SolverOptions* localsolver_options )
{
  int nc = problem->numberOfContacts;
  /* printf("fc3d_projectionOnConeWithLocalIteration_initialize. Allocation of dwork\n"); */
  if (!localsolver_options->dWork
      || localsolver_options->dWorkSize < nc)
  {
    localsolver_options->dWork = (double *)realloc(localsolver_options->dWork,
                                                   nc * sizeof(double));
    localsolver_options->dWorkSize = nc ;
  }
  for (int i = 0; i < nc; i++)
  {
    localsolver_options->dWork[i]=1.0;
  }
}

void fc3d_projectionOnConeWithLocalIteration_free(FrictionContactProblem * problem, FrictionContactProblem * localproblem, SolverOptions* localsolver_options )
{
  free(localsolver_options->dWork);
  localsolver_options->dWork=NULL;
}

int fc3d_projectionOnConeWithLocalIteration_solve(FrictionContactProblem* localproblem, double* reaction, SolverOptions* options)
{
  DEBUG_BEGIN("fc3d_projectionOnConeWithLocalIteration_solve(...)\n");
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;

  double * MLocal = localproblem->M->matrix0;
  double * qLocal = localproblem->q;
  double mu_i = localproblem->mu[0];
  /* int nLocal = 3; */


  /*   /\* Builds local problem for the current contact *\/ */
  /*   fc3d_projection_update(localproblem, reaction); */


  /*double an = 1./(MLocal[0]);*/
  /*   double alpha = MLocal[nLocal+1] + MLocal[2*nLocal+2]; */
  /*   double det = MLocal[1*nLocal+1]*MLocal[2*nLocal+2] - MLocal[2*nLocal+1] + MLocal[1*nLocal+2]; */
  /*   double beta = alpha*alpha - 4*det; */
  /*   double at = 2*(alpha - beta)/((alpha + beta)*(alpha + beta)); */

  /* double an = 1. / (MLocal[0]); */

  /* double at = 1.0 / (MLocal[4] + mu_i); */
  /* double as = 1.0 / (MLocal[8] + mu_i); */
  /* at = an; */
  /* as = an; */
  double rho=   options->dWork[options->iparam[SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_CONTACTNUMBER]] , rho_k;
  DEBUG_PRINTF (" Contact options->iparam[4] = %i\n",options->iparam[4] );
  DEBUG_PRINTF("saved rho = %14.7e\n",rho );
  assert(rho >0);



  /* int incx = 1, incy = 1; */
  int i ;


  double velocity[3],velocity_k[3],reaction_k[3],worktmp[3];
  double normUT;
  double localerror = 1.0;
  //printf ("localerror = %14.7e\n",localerror );
  int localiter = 0;
  double localtolerance = dparam[0];


  /* Variable for Line_search */
  double a1,a2;
  int success = 0;
  double localerror_k;
  int ls_iter = 0;
  int ls_itermax = 10;
  /* double tau=dparam[4], tauinv=dparam[5], L= dparam[6], Lmin = dparam[7]; */
  double tau=2.0/3.0, tauinv = 3.0/2.0,  L= 0.9, Lmin =0.3;




  /*     printf ("localtolerance = %14.7e\n",localtolerance ); */
  while ((localerror > localtolerance) && (localiter < iparam[0]))
  {
    localiter ++;

    /*    printf ("reaction[0] = %14.7e\n",reaction[0]); */
    /*    printf ("reaction[1] = %14.7e\n",reaction[1]); */
    /*    printf ("reaction[2] = %14.7e\n",reaction[2]); */

    /* Store the error */
    localerror_k = localerror;

    /* store the reaction at the beginning of the iteration */
    /* cblas_dcopy(nLocal , reaction , 1 , reaction_k, 1); */

    reaction_k[0]=reaction[0];
    reaction_k[1]=reaction[1];
    reaction_k[2]=reaction[2];
    
    /* /\* velocity_k <- q  *\/ */
    /* cblas_dcopy_msan(nLocal , qLocal , 1 , velocity_k, 1); */
    /* /\* velocity_k <- q + M * reaction  *\/ */
    /* cblas_dgemv(CblasColMajor,CblasNoTrans, nLocal, nLocal, 1.0, MLocal, 3, reaction, incx, 1.0, velocity_k, incy); */
    for (i = 0; i < 3; i++) velocity_k[i] = MLocal[i + 0 * 3] * reaction[0] + qLocal[i]
                              + MLocal[i + 1 * 3] * reaction[1] +
                              + MLocal[i + 2 * 3] * reaction[2] ;

    ls_iter = 0 ;
    success =0;
    rho_k=rho / tau;

    normUT = sqrt(velocity_k[1] * velocity_k[1] + velocity_k[2] * velocity_k[2]);
    while (!success && (ls_iter < ls_itermax))
    {
      rho_k = rho_k * tau ;
      reaction[0] = reaction_k[0] - rho_k * (velocity_k[0] + mu_i * normUT);
      reaction[1] = reaction_k[1] - rho_k * velocity_k[1];
      reaction[2] = reaction_k[2] - rho_k * velocity_k[2];


      projectionOnCone(&reaction[0], mu_i);

      /* velocity <- q  */
      /* cblas_dcopy(nLocal , qLocal , 1 , velocity, 1); */  
      /* velocity <- q + M * reaction  */
      /* cblas_dgemv(CblasColMajor,CblasNoTrans, nLocal, nLocal, 1.0, MLocal, 3, reaction, incx, 1.0, velocity, incy); */


      for (i = 0; i < 3; i++) velocity[i] = MLocal[i + 0 * 3] * reaction[0] + qLocal[i]
                                + MLocal[i + 1 * 3] * reaction[1] +
                                + MLocal[i + 2 * 3] * reaction[2] ;


      
      a1 = sqrt((velocity_k[0] - velocity[0]) * (velocity_k[0] - velocity[0]) +
                (velocity_k[1] - velocity[1]) * (velocity_k[1] - velocity[1]) +
                (velocity_k[2] - velocity[2]) * (velocity_k[2] - velocity[2]));

      a2 = sqrt((reaction_k[0] - reaction[0]) * (reaction_k[0] - reaction[0]) +
                (reaction_k[1] - reaction[1]) * (reaction_k[1] - reaction[1]) +
                (reaction_k[2] - reaction[2]) * (reaction_k[2] - reaction[2]));



      success = (rho_k*a1 <= L * a2)?1:0;

      /* printf("rho_k = %12.8e\t", rho_k); */
      /* printf("a1 = %12.8e\t", a1); */
      /* printf("a2 = %12.8e\t", a2); */
      /* printf("norm reaction = %12.8e\t",sqrt(( reaction[0]) * (reaction[0]) + */
      /*           ( reaction[1]) *  reaction[1]) + */
      /*           ( reaction[2]) * ( reaction[2])); */
      /* printf("success = %i\n", success); */

      ls_iter++;
    }

    /* printf("--  localiter = %i\t, rho= %.10e\t, error = %.10e \n", localiter, rho, localerror); */

    /* compute local error */
    localerror =0.0;
    fc3d_unitary_compute_and_add_error(reaction , velocity, mu_i, &localerror, worktmp);


    /*Update rho*/
      if ((rho_k*a1 < Lmin * a2) && (localerror < localerror_k))
      {
        rho =rho_k*tauinv;
      }
      else
        rho =rho_k;
    if (verbose > 1)
    {
      printf("--  fc3d_projectionOnConeWithLocalIteration_solve localiter = %i\t, rho= %.10e\t, error = %.10e \n", localiter, rho, localerror);

    }


  }
  options->dWork[options->iparam[4]] =rho;
  options->dparam[1] = localerror ;
  DEBUG_PRINTF("final rho  =%e\n", rho);
  
  DEBUG_END("fc3d_projectionOnConeWithLocalIteration_solve(...)\n");
  if (localerror > localtolerance)
    return 1;
  return 0;

}

void fc3d_projectionOnCylinder_initialize(FrictionContactProblem * problem, FrictionContactProblem * localproblem, SolverOptions* options)
{
  assert(localproblem);
  assert(!localproblem->mu);
  localproblem->mu = options->dWork;
}

void fc3d_projectionOnCylinder_update(int contact, FrictionContactProblem* problem, FrictionContactProblem* localproblem, double* reaction, SolverOptions* options)
{
  /* Build a local problem for a specific contact
     reaction corresponds to the global vector (size n) of the global problem.
  */

  /* Call the update function which depends on the storage for MGlobal/MBGlobal */
  /* Build a local problem for a specific contact
   reaction corresponds to the global vector (size n) of the global problem.
  */

  /* The part of MGlobal which corresponds to the current block is copied into MLocal */
  fc3d_local_problem_fill_M(problem, localproblem, contact);

  /****  Computation of qLocal = qBlock + sum over a row of blocks in MGlobal of the products MLocal.reactionBlock,
     excluding the block corresponding to the current contact. ****/
  fc3d_local_problem_compute_q(problem, localproblem, reaction, contact);

}


int fc3d_projectionOnCone_solve(FrictionContactProblem* localproblem, double* reaction, SolverOptions * options)
{


  /*  /\* Builds local problem for the current contact *\/ */
  /*   fc3d_projection_update(contact, reaction); */



  double * MLocal = localproblem->M->matrix0;
  double * qLocal = localproblem->q;
  double mu_i = localproblem->mu[0];
  /* int nLocal = 3; */

  /* this part is critical for the success of the projection */
  /*double an = 1./(MLocal[0]);*/
  /*   double alpha = MLocal[nLocal+1] + MLocal[2*nLocal+2]; */
  /*   double det = MLocal[1*nLocal+1]*MLocal[2*nLocal+2] - MLocal[2*nLocal+1] + MLocal[1*nLocal+2]; */
  /*   double beta = alpha*alpha - 4*det; */
  /*   double at = 2*(alpha - beta)/((alpha + beta)*(alpha + beta)); */

  //double an = 1./(MLocal[0]+mu_i);
  double an = 1. / (MLocal[0]);


  /* int incx = 1, incy = 1; */
  double worktmp[3];
  double normUT;
  /* cblas_dcopy_msan(nLocal , qLocal, incx , worktmp , incy); */
  /* cblas_dgemv(CblasColMajor,CblasNoTrans, nLocal, nLocal, 1.0, MLocal, 3, reaction, incx, 1.0, worktmp, incy); */

  for (int i = 0; i < 3; i++) worktmp[i] = MLocal[i + 0 * 3] * reaction[0] + qLocal[i]
                              + MLocal[i + 1 * 3] * reaction[1] +
                              + MLocal[i + 2 * 3] * reaction[2] ;

  
  normUT = sqrt(worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2]);
  reaction[0] -= an * (worktmp[0] + mu_i * normUT);
  reaction[1] -= an * worktmp[1];
  reaction[2] -= an * worktmp[2];


  projectionOnCone(reaction, mu_i);
  return 0;

}



void fc3d_projection_free(FrictionContactProblem * problem, FrictionContactProblem * localproblem, SolverOptions* localsolver_options )
{
}

void fc3d_projection_with_regularization_free(FrictionContactProblem * problem, FrictionContactProblem * localproblem, SolverOptions* localsolver_options )
{
  free(localproblem->M->matrix0);
  localproblem->M->matrix0 = NULL;
}



int fc3d_projectionOnCone_velocity_solve(FrictionContactProblem* localproblem, double* velocity,  SolverOptions* options)
{
  /* int and double parameters */
  /*     int* iparam = options->iparam; */
  /*     double* dparam = options->dparam; */
  /* Current block position */

  /* Builds local problem for the current contact */
  /*   fc3d_projection_update(contact, velocity); */

  double * MLocal = localproblem->M->matrix0;
  double * qLocal = localproblem->q;
  double mu_i = localproblem->mu[0];
  /* int nLocal = 3; */


  /*double an = 1./(MLocal[0]);*/
  /*   double alpha = MLocal[nLocal+1] + MLocal[2*nLocal+2]; */
  /*   double det = MLocal[1*nLocal+1]*MLocal[2*nLocal+2] - MLocal[2*nLocal+1] + MLocal[1*nLocal+2]; */
  /*   double beta = alpha*alpha - 4*det; */
  /*   double at = 2*(alpha - beta)/((alpha + beta)*(alpha + beta)); */

  double an = 1. / (MLocal[0]);
  int i;
  /* int incx = 1, incy = 1; */
  double worktmp[3];
  double normUT;

  /* cblas_dcopy(nLocal , qLocal, incx , worktmp , incy); */
  /* cblas_dgemv(CblasColMajor,CblasNoTrans, nLocal, nLocal, 1.0, MLocal, 3, velocity, incx, 1.0, worktmp, incy); */
  for (i = 0; i < 3; i++) worktmp[i] = MLocal[i + 0 * 3] * velocity[0] + qLocal[i]
                              + MLocal[i + 1 * 3] * velocity[1] +
                              + MLocal[i + 2 * 3] * velocity[2] ;

  
  
  normUT = sqrt(velocity[1] * velocity[1] + velocity[2] * velocity[2]);
  velocity[0] -=  - mu_i * normUT + an * (worktmp[0]);
  velocity[1] -= an * worktmp[1];
  velocity[2] -= an * worktmp[2];
  double invmui = 1.0 / mu_i;
  projectionOnCone(velocity, invmui);

  normUT = sqrt(velocity[1] * velocity[1] + velocity[2] * velocity[2]);
  velocity[0] -= mu_i * normUT;
  return 0;
}



int fc3d_projectionOnCylinder_solve(FrictionContactProblem *localproblem , double* reaction, SolverOptions* options)
{
  /* int and double parameters */
  /*   int* iparam = options->iparam; */
  /*   double* dparam = options->dparam; */
  double * MLocal = localproblem->M->matrix0;
  double * qLocal = localproblem->q;
  /* int nLocal = 3; */


  /* Builds local problem for the current contact */
  /*   fc3d_projection_update(contact, reaction); */


  /*double an = 1./(MLocal[0]);*/
  /*   double alpha = MLocal[nLocal+1] + MLocal[2*nLocal+2]; */
  /*   double det = MLocal[1*nLocal+1]*MLocal[2*nLocal+2] - MLocal[2*nLocal+1] + MLocal[1*nLocal+2]; */
  /*   double beta = alpha*alpha - 4*det; */
  /*   double at = 2*(alpha - beta)/((alpha + beta)*(alpha + beta)); */

  double an = 1. / (MLocal[0]);
  int i;
  /* int incx = 1, incy = 1; */
  double worktmp[3];

  double R  = localproblem->mu[options->iparam[4]];
  //printf("R=%e\n", R);
  /* cblas_dcopy(nLocal , qLocal, incx , worktmp , incy); */
  /* cblas_dgemv(CblasColMajor,CblasNoTrans, nLocal, nLocal, 1.0, MLocal, 3, reaction, incx, 1.0, worktmp, incy); */
  for (i = 0; i < 3; i++) worktmp[i] = MLocal[i + 0 * 3] * reaction[0] + qLocal[i]
                              + MLocal[i + 1 * 3] * reaction[1] +
                              + MLocal[i + 2 * 3] * reaction[2] ;
  reaction[0] -= an * worktmp[0];
  reaction[1] -= an * worktmp[1];
  reaction[2] -= an * worktmp[2];

  projectionOnCylinder(reaction, R);
  return 0;

}

void fc3d_projectionOnCylinderWithLocalIteration_initialize(
  FrictionContactProblem * problem, FrictionContactProblem * localproblem,
  SolverOptions* options, SolverOptions* localsolver_options )
{
  int nc = problem->numberOfContacts;
  /* printf("fc3d_projectionOnConeWithLocalIteration_initialize. Allocation of dwork\n"); */
  if(localproblem->mu)
  {
    free(localproblem->mu);
  }

  localproblem->mu = options->dWork;

  if (!localsolver_options->dWork)
  {
    localsolver_options->dWork = (double *)malloc(nc * sizeof(double));
    localsolver_options->dWorkSize = nc;
  }
  else
  {
    fprintf(stderr, "Numerics, fc3d_projectionOnCylinderWithLocalIteration_initialize failed. localsolver_options->dWork is different from NULL.\n");
    exit(EXIT_FAILURE);
  }
  for (int i = 0; i < nc; i++)
  {
    localsolver_options->dWork[i]=1.0;
  }
}
void fc3d_projectionOnCylinderWithLocalIteration_free(FrictionContactProblem * problem, FrictionContactProblem * localproblem, SolverOptions* localsolver_options )
{
  localproblem->mu = NULL;
  free(localsolver_options->dWork);
  localsolver_options->dWork=NULL;
}

void fc3d_projectionOnCylinder_free(FrictionContactProblem * problem, FrictionContactProblem * localproblem, SolverOptions* localsolver_options )
{
  localproblem->mu = NULL;
}

int fc3d_projectionOnCylinderWithLocalIteration_solve(FrictionContactProblem* localproblem, double* reaction, SolverOptions* options)
{
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;

  double * MLocal = localproblem->M->matrix0;
  double * qLocal = localproblem->q;
  /* int nLocal = 3; */

  /*   /\* Builds local problem for the current contact *\/ */
  /*   fc3d_projection_update(localproblem, reaction); */

  /*double an = 1./(MLocal[0]);*/
  /*   double alpha = MLocal[nLocal+1] + MLocal[2*nLocal+2]; */
  /*   double det = MLocal[1*nLocal+1]*MLocal[2*nLocal+2] - MLocal[2*nLocal+1] + MLocal[1*nLocal+2]; */
  /*   double beta = alpha*alpha - 4*det; */
  /*   double at = 2*(alpha - beta)/((alpha + beta)*(alpha + beta)); */

  /* double an = 1. / (MLocal[0]); */

  /* double at = 1.0 / (MLocal[4] + mu_i); */
  /* double as = 1.0 / (MLocal[8] + mu_i); */
  /* at = an; */
  /* as = an; */
  double rho=   options->dWork[options->iparam[4]] , rho_k;
  /* printf ("saved rho = %14.7e\n",rho );  */
  /* printf ("options->iparam[4] = %i\n",options->iparam[4] );  */



  /* int incx = 1, incy = 1; */
  int i;
  double velocity[3],velocity_k[3],reaction_k[3], worktmp[3];

  double localerror = 1.0;
  //printf ("localerror = %14.7e\n",localerror );
  int localiter = 0;
  double localtolerance = dparam[0];


  /* Variable for Line_search */
  double a1,a2;
  int success = 0;
  double localerror_k;
  int ls_iter = 0;
  int ls_itermax = 10;
  /* double tau=dparam[4], tauinv=dparam[5], L= dparam[6], Lmin = dparam[7]; would be better */
  double tau=2.0/3.0, tauinv = 3.0/2.0,  L= 0.9, Lmin =0.3;

  double R  = localproblem->mu[options->iparam[4]];

  /* printf ("R = %14.7e\n",R ); */
  while ((localerror > localtolerance) && (localiter < iparam[0]))
  {
    localiter ++;

    /*    printf ("reaction[0] = %14.7e\n",reaction[0]); */
    /*    printf ("reaction[1] = %14.7e\n",reaction[1]); */
    /*    printf ("reaction[2] = %14.7e\n",reaction[2]); */

    /* Store the error */
    localerror_k = localerror;

    /* /\* store the reaction at the beginning of the iteration *\/ */
    /* cblas_dcopy_msan(nLocal , reaction , 1 , reaction_k, 1); */

    /* /\* velocity_k <- q  *\/ */
    /* cblas_dcopy_msan(nLocal , qLocal , 1 , velocity_k, 1); */

    /* /\* velocity_k <- q + M * reaction  *\/ */
    /* cblas_dgemv(CblasColMajor,CblasNoTrans, nLocal, nLocal, 1.0, MLocal, 3, reaction, incx, 1.0, velocity_k, incy); */
    reaction_k[0]=reaction[0];
    reaction_k[1]=reaction[1];
    reaction_k[2]=reaction[2];
    
    /* /\* velocity_k <- q  *\/ */
    /* cblas_dcopy_msan(nLocal , qLocal , 1 , velocity_k, 1); */
    /* /\* velocity_k <- q + M * reaction  *\/ */
    /* cblas_dgemv(CblasColMajor,CblasNoTrans, nLocal, nLocal, 1.0, MLocal, 3, reaction, incx, 1.0, velocity_k, incy); */
    for (i = 0; i < 3; i++) velocity_k[i] = MLocal[i + 0 * 3] * reaction[0] + qLocal[i]
                              + MLocal[i + 1 * 3] * reaction[1] +
                              + MLocal[i + 2 * 3] * reaction[2] ;


    ls_iter = 0 ;
    success =0;
    rho_k=rho / tau;


    while (!success && (ls_iter < ls_itermax))
    {
      rho_k = rho_k * tau ;

      reaction[0] = reaction_k[0] - rho_k * velocity_k[0];
      reaction[1] = reaction_k[1] - rho_k * velocity_k[1];
      reaction[2] = reaction_k[2] - rho_k * velocity_k[2];

      projectionOnCylinder(&reaction[0], R);

      /* /\* velocity <- q  *\/ */
      /* cblas_dcopy(nLocal , qLocal , 1 , velocity, 1); */
      /* /\* velocity <- q + M * reaction  *\/ */
      /* cblas_dgemv(CblasColMajor,CblasNoTrans, nLocal, nLocal, 1.0, MLocal, 3, reaction, incx, 1.0, velocity, incy); */

      for (i = 0; i < 3; i++) velocity[i] = MLocal[i + 0 * 3] * reaction[0] + qLocal[i]
                                + MLocal[i + 1 * 3] * reaction[1] +
                                + MLocal[i + 2 * 3] * reaction[2] ;
      
      a1 = sqrt((velocity_k[0] - velocity[0]) * (velocity_k[0] - velocity[0]) +
                (velocity_k[1] - velocity[1]) * (velocity_k[1] - velocity[1]) +
                (velocity_k[2] - velocity[2]) * (velocity_k[2] - velocity[2]));

      a2 = sqrt((reaction_k[0] - reaction[0]) * (reaction_k[0] - reaction[0]) +
                (reaction_k[1] - reaction[1]) * (reaction_k[1] - reaction[1]) +
                (reaction_k[2] - reaction[2]) * (reaction_k[2] - reaction[2]));



      success = (rho_k*a1 <= L * a2)?1:0;

      /* printf("rho_k = %12.8e\t", rho_k); */
      /* printf("a1 = %12.8e\t", a1); */
      /* printf("a2 = %12.8e\t", a2); */
      /* printf("norm reaction = %12.8e\t",sqrt(( reaction[0]) * (reaction[0]) + */
      /*           ( reaction[1]) *  reaction[1]) + */
      /*           ( reaction[2]) * ( reaction[2])); */
      /* printf("success = %i\n", success); */

      ls_iter++;
    }
    /* if (verbose>2) */
    /*   printf("--  localiter = %i\t, rho= %.10e\t, error = %.10e \n", localiter, rho, localerror);  */

    /* compute local error */
    localerror =0.0;
    fc3d_Tresca_unitary_compute_and_add_error(reaction , velocity, R, &localerror,worktmp);


    /*Update rho*/
      if ((rho_k*a1 < Lmin * a2) && (localerror < localerror_k))
      {
        rho =rho_k*tauinv;
      }
      else
        rho =rho_k;

    if (verbose > 1)
    {
      printf("--  localiter = %i\t, rho= %.10e\t, error = %.10e \n", localiter, rho, localerror);

    }

  }
  options->dWork[options->iparam[4]] =rho;

 
 if (verbose > 1)
    {
      printf("--  localiter = %i\t, rho= %.10e\t, error = %.10e \n", localiter, rho, localerror);

    }

  if (localerror > localtolerance)
    return 1;
  return 0;

}
int fc3d_projectionOnConeWithDiagonalization_setDefaultSolverOptions(SolverOptions* options)
{

  numerics_printf("Set the Default SolverOptions for the ONECONTACT_ProjectionOnConeWithDiagonalization  Solver\n");
  

  options->solverId = SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithDiagonalization;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 0;
  options->dSize = 0;
  options->iSize = 10;
  options->dSize = 10;
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  solver_options_nullify(options);


  return 0;
}
int fc3d_projectionOnConeWithRegularization_setDefaultSolverOptions(SolverOptions* options)
{

  numerics_printf("Set the Default SolverOptions for the ONECONTACT_ProjectionOnConeWithRegularization Solver\n");
  

  options->solverId = SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 0;
  options->dSize = 0;
  options->iSize = 10;
  options->dSize = 10;
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  solver_options_nullify(options);


  return 0;
}

int fc3d_projectionOnCone_setDefaultSolverOptions(SolverOptions* options)
{

  numerics_printf("Set the Default SolverOptions for the ONECONTACT_ProjectionOnCone  Solver\n");
  

  options->solverId = SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 10;
  options->dSize = 10;
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  solver_options_nullify(options);

  return 0;
}
int fc3d_projectionOnCone_velocity_setDefaultSolverOptions(SolverOptions* options)
{

  numerics_printf("Set the Default SolverOptions for the ONECONTACT_ProjectionOnCone_velocity  Solver\n");
  

  options->solverId = SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone_velocity;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 10;
  options->dSize = 10;
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  solver_options_nullify(options);


  return 0;
}
int fc3d_projectionOnConeWithLocalIteration_setDefaultSolverOptions(SolverOptions* options)
{

  numerics_printf("Set the Default SolverOptions for the ONECONTACT_ProjectionOnConeWithLocalIteration  Solver\n");
  

  options->solverId = SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 20;
  options->dSize = 20;
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  solver_options_nullify(options);

  options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
  options->dparam[SICONOS_DPARAM_TOL] = 1e-14;

  return 0;
}
int fc3d_projectionOnCylinder_setDefaultSolverOptions(SolverOptions* options)
{

  numerics_printf("Set the Default SolverOptions for the ONECONTACT_projectionOnCylinder  Solver\n");
  

  options->solverId = SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinder;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 10;
  options->dSize = 10;
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  solver_options_nullify(options);


  return 0;
}
int fc3d_projectionOnCylinderWithLocalIteration_setDefaultSolverOptions(SolverOptions* options)
{

  numerics_printf("Set the Default SolverOptions for the ONECONTACT_ProjectionOnCylinderWithLocalIteration  Solver\n");
  

  options->solverId = SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinderWithLocalIteration;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 20;
  options->dSize = 20;
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  solver_options_nullify(options);

  options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
  options->dparam[SICONOS_DPARAM_TOL] = 1e-14;

  return 0;
}
