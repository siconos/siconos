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

#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"
#endif

void rolling_fc3d_projection_update(int contact, RollingFrictionContactProblem* problem,
                                    RollingFrictionContactProblem* localproblem, double* reaction, SolverOptions* options)
{
  /* Build a local problem for a specific contact
     reaction corresponds to the global vector (size n) of the global problem.
  */

  /* Call the update function which depends on the storage for MGlobal/MBGlobal */
  /* Build a local problem for a specific contact
   reaction corresponds to the global vector (size n) of the global problem.
  */

  /* The part of MGlobal which corresponds to the current block is copied into MLocal */
  rolling_fc3d_local_problem_fill_M(problem, localproblem, contact);

  /****  Computation of qLocal = qBlock + sum over a row of blocks in MGlobal of the products MLocal.reactionBlock,
     excluding the block corresponding to the current contact. ****/
  rolling_fc3d_local_problem_compute_q(problem, localproblem, reaction, contact);

  /* Friction coefficient for current block*/
  localproblem->mu[0] = problem->mu[contact];
  /* Rolling Friction coefficient for current block*/
  localproblem->mu_r[0] = problem->mu_r[contact];

}



void rolling_fc3d_projection_initialize(RollingFrictionContactProblem * problem,
                                        RollingFrictionContactProblem * localproblem)
{
}

void rolling_fc3d_projection_free(RollingFrictionContactProblem * problem,
                                  RollingFrictionContactProblem * localproblem,
                                  SolverOptions* localsolver_options )
{
}

int rolling_fc3d_projectionOnCone_solve(
  RollingFrictionContactProblem* localproblem, double* reaction, SolverOptions * options)
{
  DEBUG_BEGIN("rolling_fc3d_projectionOnCone_solve(...)\n");



  double * MLocal = localproblem->M->matrix0;
  double * qLocal = localproblem->q;
  double mu_i = localproblem->mu[0];
  double mu_r_i = localproblem->mu_r[0];
  /* int nLocal = 3; */

  /* this part is critical for the success of the projection */
  /*double an = 1./(MLocal[0]);*/
  /*   double alpha = MLocal[nLocal+1] + MLocal[2*nLocal+2]; */
  /*   double det = MLocal[1*nLocal+1]*MLocal[2*nLocal+2] - MLocal[2*nLocal+1] + MLocal[1*nLocal+2]; */
  /*   double beta = alpha*alpha - 4*det; */
  /*   double at = 2*(alpha - beta)/((alpha + beta)*(alpha + beta)); */

  //double an = 1./(MLocal[0]+mu_i);
  //double an = 1. / (MLocal[0]);

  double an=1.0;
  
  /* int incx = 1, incy = 1; */
  double velocity[5];
  double normUT, normOmegaT;
  /* cblas_dcopy_msan(nLocal , qLocal, incx , velocity , incy); */
  /* cblas_dgemv(CblasColMajor,CblasNoTrans, nLocal, nLocal, 1.0, MLocal, 3, reaction, incx, 1.0, velocity, incy); */
   for (int i = 0; i < 5; i++) velocity[i]
                              = MLocal[i + 0 * 5] * reaction[0] + qLocal[i]
                              + MLocal[i + 1 * 5] * reaction[1]
                              + MLocal[i + 2 * 5] * reaction[2]
                              + MLocal[i + 3 * 5] * reaction[3]
                              + MLocal[i + 4 * 5] * reaction[4];
   DEBUG_EXPR(NV_display(velocity,5););

  
  normUT = hypot(velocity[1], velocity[2]);
  normOmegaT = hypot(velocity[3], velocity[4]);
  reaction[0] -= an * (velocity[0] + mu_i * normUT+ mu_r_i*normOmegaT);
  reaction[1] -= an * velocity[1];
  reaction[2] -= an * velocity[2];
  reaction[3] -= an * velocity[3];
  reaction[4] -= an * velocity[4];

#ifdef DEBUG_MESSAGES
  display_status_rolling_cone(projectionOnRollingCone(reaction, mu_i, mu_r_i));
#else
  projectionOnRollingCone(reaction, mu_i, mu_r_i);
#endif
  
  DEBUG_EXPR(NV_display(reaction,5););
  DEBUG_END("rolling_fc3d_projectionOnCone_solve(...)\n");
  return 0;

}




void rolling_fc3d_projectionOnConeWithLocalIteration_initialize(RollingFrictionContactProblem * problem, RollingFrictionContactProblem * localproblem, SolverOptions* localsolver_options )
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

void rolling_fc3d_projectionOnConeWithLocalIteration_free(RollingFrictionContactProblem * problem, RollingFrictionContactProblem * localproblem, SolverOptions* localsolver_options )
{
  free(localsolver_options->dWork);
  localsolver_options->dWork=NULL;
}

int rolling_fc3d_projectionOnConeWithLocalIteration_solve(RollingFrictionContactProblem* localproblem, double* reaction, SolverOptions* options)
{
  DEBUG_BEGIN("rolling_fc3d_projectionOnConeWithLocalIteration_solve(...)\n");
  //DEBUG_EXPR(rollingFrictionContact_display(localproblem););
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;

  double * MLocal = localproblem->M->matrix0;
  double * qLocal = localproblem->q;
  double mu_i = localproblem->mu[0];
  double mu_r_i = localproblem->mu_r[0];  /* int nLocal = 5; */


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

  double rho=   options->dWork[options->iparam[SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_CONTACTNUMBER]];
  DEBUG_PRINTF ("Contact options->iparam[4] = %i,\t",options->iparam[4] );
  DEBUG_PRINTF("saved rho = %14.7e\n",rho );
  assert(rho >0.0);

  /* int incx = 1, incy = 1; */
  int i ;


  double velocity[5],velocity_k[5],reaction_k[5],worktmp[5];

  double localerror = 1.0;
  //printf ("localerror = %14.7e\n",localerror );
  int localiter = 0;
  double localtolerance = dparam[0];

  /* Variable for Line_search */
  double a1,a2;

  /* double tau=dparam[4], tauinv=dparam[5], L= dparam[6], Lmin = dparam[7]; */
  double tau=2.0/3.0, tauinv = 3.0/2.0,  L= 0.9, Lmin =0.3;

  /*     printf ("localtolerance = %14.7e\n",localtolerance ); */
  while ((localerror > localtolerance) && (localiter < iparam[0]))
  {
    DEBUG_PRINT("\n Local iteration starts \n");
    localiter ++;

    double rho_k;
    double normUT, normOmegaT;
    double localerror_k;
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
    reaction_k[3]=reaction[3];
    reaction_k[4]=reaction[4];
    DEBUG_EXPR(NV_display(reaction_k,5););


    /* /\* velocity_k <- q  *\/ */
    /* cblas_dcopy_msan(nLocal , qLocal , 1 , velocity_k, 1); */
    /* /\* velocity_k <- q + M * reaction  *\/ */
    /* cblas_dgemv(CblasColMajor,CblasNoTrans, nLocal, nLocal, 1.0,
       MLocal, 3, reaction, incx, 1.0, velocity_k, incy); */

    for (i = 0; i < 5; i++) velocity_k[i]
                              = MLocal[i + 0 * 5] * reaction[0] + qLocal[i]
                              + MLocal[i + 1 * 5] * reaction[1]
                              + MLocal[i + 2 * 5] * reaction[2]
                              + MLocal[i + 3 * 5] * reaction[3]
                              + MLocal[i + 4 * 5] * reaction[4];
    DEBUG_EXPR(NV_display(velocity_k,5););

    int ls_iter = 0;
    int ls_itermax = 10;
    int success = 0;
    rho_k=rho / tau;

    /* normUT = sqrt(velocity_k[1] * velocity_k[1] + velocity_k[2] * velocity_k[2]); */
    /* normOmegaT = sqrt(velocity_k[3] * velocity_k[3] + velocity_k[4] * velocity_k[4]); */
    normUT = hypot(velocity_k[1], velocity_k[2]);
    normOmegaT = hypot(velocity_k[3], velocity_k[4]);
    /* ls_itermax=1; */
    DEBUG_PRINTF("ls_itermax =%i\n", ls_itermax);
    while (!success && (ls_iter < ls_itermax))
    {
      rho_k = rho_k * tau ;
      DEBUG_PRINTF("rho_k =%f\n", rho_k);
      reaction[0] = reaction_k[0] - rho_k * (velocity_k[0] + mu_i * normUT + mu_r_i*normOmegaT);
      reaction[1] = reaction_k[1] - rho_k * velocity_k[1];
      reaction[2] = reaction_k[2] - rho_k * velocity_k[2];
      reaction[3] = reaction_k[3] - rho_k * velocity_k[3];
      reaction[4] = reaction_k[4] - rho_k * velocity_k[4];
      DEBUG_PRINT("r-rho tilde v before projection")
      DEBUG_EXPR(NV_display(reaction,5););
/* #ifdef DEBUG_MESSAGES */
/*       display_status_rolling_cone(projectionOnRollingCone(&reaction[0], mu_i, mu_r_i)); */
/* #else */
      projectionOnRollingCone(&reaction[0], mu_i, mu_r_i);
/* #endif */
      DEBUG_EXPR(NV_display(reaction,5););
      /* velocity <- q  */
      /* cblas_dcopy(nLocal , qLocal , 1 , velocity, 1); */
      /* velocity <- q + M * reaction  */
      /* cblas_dgemv(CblasColMajor,CblasNoTrans, nLocal, nLocal, 1.0, MLocal, 3, reaction, incx, 1.0, velocity, incy); */


      for (i = 0; i < 5; i++) velocity[i]
                                = MLocal[i + 0 * 5] * reaction[0] + qLocal[i]
                                + MLocal[i + 1 * 5] * reaction[1]
                                + MLocal[i + 2 * 5] * reaction[2]
                                + MLocal[i + 3 * 5] * reaction[3]
                                + MLocal[i + 4 * 5] * reaction[4];
      DEBUG_EXPR(NV_display(velocity,5););


      a1 = sqrt(
        (velocity_k[0] - velocity[0]) * (velocity_k[0] - velocity[0]) +
        (velocity_k[1] - velocity[1]) * (velocity_k[1] - velocity[1]) +
        (velocity_k[2] - velocity[2]) * (velocity_k[2] - velocity[2]) +
        (velocity_k[3] - velocity[3]) * (velocity_k[3] - velocity[3]) +
        (velocity_k[4] - velocity[4]) * (velocity_k[4] - velocity[4]));

      a2 = sqrt(
        (reaction_k[0] - reaction[0]) * (reaction_k[0] - reaction[0]) +
        (reaction_k[1] - reaction[1]) * (reaction_k[1] - reaction[1]) +
        (reaction_k[2] - reaction[2]) * (reaction_k[2] - reaction[2]) +
        (reaction_k[3] - reaction[3]) * (reaction_k[3] - reaction[3]) +
        (reaction_k[4] - reaction[4]) * (reaction_k[4] - reaction[4]));



      success = (rho_k*a1 <= L * a2)?1:0;

      DEBUG_PRINTF("rho_k = %12.8e\t", rho_k);
      DEBUG_PRINTF("a1 = %12.8e\t", a1);
      DEBUG_PRINTF("a2 = %12.8e\t", a2);
      DEBUG_PRINTF("norm reaction = %12.8e\t",
                   sqrt(reaction[0] * reaction[0] +
                        reaction[1] * reaction[1] +
                        reaction[2] * reaction[2] +
                        reaction[3] * reaction[3] +
                        reaction[4] * reaction[4]
                     ) );
      DEBUG_PRINTF("success = %i\n", success);

      ls_iter++;
    }

    /* printf("--  localiter = %i\t, rho= %.10e\t, error = %.10e \n", localiter, rho, localerror); */

    /* compute local error */
    localerror =0.0;
    rolling_fc3d_unitary_compute_and_add_error(reaction , velocity,
                                               mu_i, mu_r_i,
                                               &localerror, worktmp);
    DEBUG_EXPR(NV_display(reaction,5););
    DEBUG_EXPR(NV_display(velocity,5););



    /*Update rho*/
      if ((rho_k*a1 < Lmin * a2) && (localerror < localerror_k))
      {
        rho =rho_k*tauinv;
      }
      else
        rho =rho_k;

      /* rho_k=1.0; */
      /* rho=1.0 */

        ;
    if (verbose > 1)
    {
      printf("--  rolling_fc3d_projectionOnConeWithLocalIteration_solve localiter = %i\t, rho= %.10e\t, error = %.10e \n", localiter, rho, localerror);

    }
    //getchar();

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
int rolling_fc3d_projectionOnCone_setDefaultSolverOptions(SolverOptions* options)
{

  numerics_printf("Set the Default SolverOptions for the ONECONTACT_ProjectionOnCone  Solver\n");
  

  options->solverId = SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnCone;
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
