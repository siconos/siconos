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


#include "rolling_fc3d_Solvers.h"
#include "projectionOnCone.h"
#include "projectionOnRollingCone.h"

#include "rolling_fc3d_compute_error.h"

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


void rolling_fc3d_projection_update(int contact, FrictionContactProblem* problem, FrictionContactProblem* localproblem, double* reaction, SolverOptions* options)
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





void rolling_fc3d_projectionOnConeWithLocalIteration_initialize(FrictionContactProblem * problem, FrictionContactProblem * localproblem, SolverOptions* localsolver_options )
{
  int nc = problem->numberOfContacts;
  /* printf("rolling_fc3d_projectionOnConeWithLocalIteration_initialize. Allocation of dwork\n"); */
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

void rolling_fc3d_projectionOnConeWithLocalIteration_free(FrictionContactProblem * problem, FrictionContactProblem * localproblem, SolverOptions* localsolver_options )
{
  free(localsolver_options->dWork);
  localsolver_options->dWork=NULL;
}

int rolling_fc3d_projectionOnConeWithLocalIteration_solve(FrictionContactProblem* localproblem, double* reaction, SolverOptions* options)
{
  DEBUG_BEGIN("rolling_fc3d_projectionOnConeWithLocalIteration_solve(...)\n");
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;

  double * MLocal = localproblem->M->matrix0;
  double * qLocal = localproblem->q;
  double mu_i = localproblem->mu[0];
  /* int nLocal = 3; */


  /*   /\* Builds local problem for the current contact *\/ */
  /*   rolling_fc3d_projection_update(localproblem, reaction); */


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
    rolling_fc3d_unitary_compute_and_add_error(reaction , velocity, mu_i, &localerror, worktmp);


    /*Update rho*/
      if ((rho_k*a1 < Lmin * a2) && (localerror < localerror_k))
      {
        rho =rho_k*tauinv;
      }
      else
        rho =rho_k;
    if (verbose > 1)
    {
      printf("--  rolling_fc3d_projectionOnConeWithLocalIteration_solve localiter = %i\t, rho= %.10e\t, error = %.10e \n", localiter, rho, localerror);

    }


  }
  options->dWork[options->iparam[4]] =rho;
  options->dparam[1] = localerror ;
  DEBUG_PRINTF("final rho  =%e\n", rho);

  DEBUG_END("rolling_fc3d_projectionOnConeWithLocalIteration_solve(...)\n");
  if (localerror > localtolerance)
    return 1;
  return 0;

}

int rolling_fc3d_projectionOnConeWithLocalIteration_setDefaultSolverOptions(SolverOptions* options)
{

  numerics_printf("Set the Default SolverOptions for the ONECONTACT_ProjectionOnConeWithLocalIteration  Solver\n");


  options->solverId = SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration;
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
