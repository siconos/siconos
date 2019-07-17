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
#include "op3x3.h"

#include "SiconosBlas.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "SiconosLapack.h"
#include "sanitizer.h"
#include "numerics_verbose.h"

/* #define DEBUG_NOCOLOR */
/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "debug.h"

#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"
#endif
#include "NumericsVector.h"
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

  normUT = sqrt(velocity[1] * velocity[1] + velocity[2] * velocity[2]);
  normOmegaT = sqrt(velocity[3] * velocity[3] + velocity[4] * velocity[4]);
  /* hypot of libm is sure but really slow */
  /* normUT = hypot(velocity[1], velocity[2]); */
  /* normOmegaT = hypot(velocity[3], velocity[4]); */
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




void rolling_fc3d_projectionOnConeWithLocalIteration_initialize(
  RollingFrictionContactProblem * problem,
  RollingFrictionContactProblem * localproblem,
  SolverOptions* localsolver_options )
{
  int nc = problem->numberOfContacts;
  localsolver_options->iparam[17]=nc;
  /* printf("rolling_fc3d_projectionOnConeWithLocalIteration_initialize. Allocation of dwork\n"); */
  if (!localsolver_options->dWork
      || localsolver_options->dWorkSize < nc)
  {
    if (localsolver_options->iparam[SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_IPARAM_USE_TRIVIAL_SOLUTION] ==
        SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION_TRUE)
    {
      /* nc for rho parameter + 30 * nc for inverse of diagonal matrix */
      unsigned int d_size = nc + 30 *nc;
      localsolver_options->dWork = (double *)realloc(localsolver_options->dWork,
                                                     d_size * sizeof(double));
      localsolver_options->dWorkSize = d_size ;
      unsigned int i_size = 5 *nc;
      localsolver_options->iWork = (int *)realloc(localsolver_options->iWork,
                                                     i_size * sizeof(int));
      localsolver_options->iWorkSize = i_size ;

    }
    else
    {
      /* nc for rho parameter */
      unsigned int size = nc ;
      localsolver_options->dWork = (double *)realloc(localsolver_options->dWork,
                                                    size * sizeof(double));
      localsolver_options->dWorkSize = size ;
    }
  }

  for (int i = 0; i < nc; i++)
  {
    localsolver_options->dWork[i]=1.0;

    if (localsolver_options->iparam[SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_IPARAM_USE_TRIVIAL_SOLUTION] ==
        SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION_TRUE)
    {
      /* The part of MGlobal which corresponds to the current block is copied into MLocal */
      rolling_fc3d_local_problem_fill_M(problem, localproblem, i);
      double * MLocal = localproblem->M->matrix0;
      unsigned int pos = i * 30 + nc;
      double * MLocal_LU = &(localsolver_options->dWork[pos]);
      memcpy(MLocal_LU, MLocal , 25*sizeof(double));
      int * ipiv =  &(localsolver_options->iWork[i*5]);
      int info =0;
      DGETRF( 5, 5, MLocal_LU, 5, ipiv, &info );
      assert(!info);
    }
  }
};


void rolling_fc3d_projectionOnConeWithLocalIteration_free(RollingFrictionContactProblem * problem, RollingFrictionContactProblem * localproblem, SolverOptions* localsolver_options )
{
  free(localsolver_options->dWork);
  localsolver_options->dWork=NULL;
}


static int  rolling_fc3d_check_trivial_solution(
  unsigned int contact,
  unsigned int nc,
  double * q,
  double mu, double mur,
  double * reaction,
  SolverOptions *localsolver_options)
{

  if (q[0] >0)
  {
    reaction[0]=0.0;
    reaction[1]=0.0;
    reaction[2]=0.0;
    reaction[3]=0.0;
    reaction[4]=0.0;
    return 1;
  }

  unsigned int pos = contact * 30 + nc;
  double * MLocal_LU = &(localsolver_options->dWork[pos]);
  int * ipiv = &(localsolver_options->iWork[5*contact]);
  /* double * tmp = localsolver_options->dWork[pos+25]; */
  memcpy(reaction, q , 5*sizeof(double));
  int info =0;
  DGETRS( LA_NOTRANS, 5, 1, MLocal_LU, 5, ipiv, reaction, 5, &info );
  assert(!info);

  /* DEBUG_EXPR(NM_dense_display(M, 5,5,5);); */
  /* double * A =  (double * ) malloc(25*sizeof(double)); */

  /* memcpy(A, M , 25*sizeof(double)); */
  /* DEBUG_EXPR(NM_dense_display(A, 5,5,5);); */
  /* double tmp[5];  */

  /* memcpy(tmp, q , 5*sizeof(double)); */
  /* DEBUG_EXPR(NV_display(q,5);); */

  /* DEBUG_EXPR(NV_display(tmp,5);); */

  /* DEBUG_PRINT("solve\n"); */
  /* /\* solve_nxn_gepp(5, A, tmp, reaction); *\/ */

  /* int ipiv[25]; */
  /* int info=0; */
  /* DGESV(5, 1, A, 5, ipiv, tmp, 5, &info); */
  /* memcpy(reaction, tmp , 5*sizeof(double)); */


  /* for (int i = 0; i < 5; i++) tmp[i] */
  /*                           = M[i + 0 * 5] * reaction[0] - q[i] */
  /*                           + M[i + 1 * 5] * reaction[1] */
  /*                           + M[i + 2 * 5] * reaction[2] */
  /*                           + M[i + 3 * 5] * reaction[3] */
  /*                           + M[i + 4 * 5] * reaction[4]; */
  /* /\* NV_display(tmp,5) *\/ */
  /* DEBUG_EXPR(NV_display(reaction,5);); */
  /* if (cblas_dnrm2(5, tmp, 1) >= 1e-10) */
  /* { */
  /*   NM_dense_display(M, 5,5,5); */
  /*   NM_dense_display(MLocal_LU, 5,5,5); */
  /* } */

  /* assert(cblas_dnrm2(5, tmp, 1) < 1e-10); */
  double normT = sqrt(reaction[1] * reaction[1] + reaction[2] * reaction[2]);
  double normMT = sqrt(reaction[3] * reaction[3] + reaction[4] * reaction[4]);

  if ((normT <= mu * -reaction [0]) && (normMT <= mur * -reaction[0]))
  {
    reaction[0]=-reaction[0];
    reaction[1]=-reaction[1];
    reaction[2]=-reaction[2];
    reaction[3]=-reaction[3];
    reaction[4]=-reaction[4];
    return 2;
  }
  return 0;

}


int rolling_fc3d_projectionOnConeWithLocalIteration_solve(
  RollingFrictionContactProblem* localproblem,
  double* reaction,
  SolverOptions* options)
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

  DEBUG_EXPR(NM_dense_display(MLocal, 5,5,5););

  double rho=   options->dWork[options->iparam[SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_CONTACTNUMBER]];
  DEBUG_PRINTF ("Contact number = %i,\t",options->iparam[4] );
  DEBUG_PRINTF("saved rho = %14.7e\n",rho );
  assert(rho >0.0);

  /* int incx = 1, incy = 1; */
  int i ;

  double velocity[5], velocity_k[5], reaction_k[5];

  //double trivial_error=0.0;
  int trivial = rolling_fc3d_check_trivial_solution(options->iparam[4],
                                                    options->iparam[17],
                                                    qLocal,
                                                    mu_i, mu_r_i,
                                                    reaction_k,
                                                    options);
  if (trivial)
  {
    /* rolling_fc3d_unitary_compute_and_add_error(reaction_k , velocity_k, */
    /*                                            mu_i, mu_r_i, */
    /*                                            &trivial_error, worktmp); */
    /* assert(trivial_error < 1e-14); */
    numerics_printf_verbose(2, "found trivial solution = %i\t error = %e", trivial );
    /* printf( "found trivial solution = %i\t error = %e\n", trivial, trivial_error );  */
    options->dparam[1] = 0.0 ;
    memcpy(reaction, reaction_k , 5*sizeof(double));
    /* NV_display(reaction,5); */
    /* NV_display(velocity_k,5); */
    DEBUG_END("rolling_fc3d_projectionOnConeWithLocalIteration_solve(...)\n");
    return 0 ;
  }
  else
    numerics_printf_verbose(2, "trivial solution not found = %i", trivial );

  double localerror = 1.0;
  //printf ("localerror = %14.7e\n",localerror );
  int localiter = 0;
  double localtolerance = dparam[0];

  /* Variable for Line_search */
  double a1,a2;

  /* double tau=dparam[4], tauinv=dparam[5], L= dparam[6], Lmin = dparam[7]; */
  double tau=2.0/3.0, tauinv = 3.0/2.0,  L= 0.9, Lmin =0.3;

  int status = -1;
  numerics_printf_verbose(2,"--  rolling_fc3d_projectionOnConeWithLocalIteration_solve contact = %i", options->iparam[4] );
  numerics_printf_verbose(2,"--  rolling_fc3d_projectionOnConeWithLocalIteration_solve | localiter \t| rho \t\t\t| error\t\t\t| status\t|");
  numerics_printf_verbose(2,"--                                                        | %i \t\t| %.10e\t| %.10e\t|%i \t\t|", localiter, rho, localerror, status);

  /*     printf ("localtolerance = %14.7e\n",localtolerance ); */
  while ((localerror > localtolerance) && (localiter < iparam[0]))
  {
    DEBUG_PRINT("\n Local iteration starts \n");
    localiter ++;

    double rho_k;
    double normUT, normOmegaT;
    double localerror_k;

    /* Store the error */
    localerror_k = localerror;

    /* store the reaction at the beginning of the iteration */
    /* cblas_dcopy(nLocal , reaction , 1 , reaction_k, 1); */
    memcpy(reaction_k, reaction , 5*sizeof(double));
    /* reaction_k[0]=reaction[0]; */
    /* reaction_k[1]=reaction[1]; */
    /* reaction_k[2]=reaction[2]; */
    /* reaction_k[3]=reaction[3]; */
    /* reaction_k[4]=reaction[4]; */
     /* DEBUG_EXPR(NV_display(reaction_k,5);); */


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

    /* memcpy(velocity_k, qLocal , 5*sizeof(double)); */
    /* mvp5x5(MLocal, reaction, velocity_k); */
    /* DEBUG_EXPR(NV_display(velocity_k,5);); */

    int ls_iter = 0;
    int ls_itermax = 100;
    int success = 0;

    rho_k=rho / tau;

    normUT = sqrt(velocity_k[1] * velocity_k[1] + velocity_k[2] * velocity_k[2]);
    normOmegaT = sqrt(velocity_k[3] * velocity_k[3] + velocity_k[4] * velocity_k[4]);
    /* hypot of libm is sure but really slow */
    /* normUT = hypot(velocity_k[1], velocity_k[2]); */
    /* normOmegaT = hypot(velocity_k[3], velocity_k[4]); */
    /* ls_itermax=1; */
    /* DEBUG_PRINTF("ls_itermax =%i\n", ls_itermax); */
    while (!success && (ls_iter < ls_itermax))
    {
      rho_k = rho_k * tau ;
      reaction[0] = reaction_k[0] - rho_k * (velocity_k[0] + mu_i * normUT + mu_r_i*normOmegaT);
      reaction[1] = reaction_k[1] - rho_k * velocity_k[1];
      reaction[2] = reaction_k[2] - rho_k * velocity_k[2];
      reaction[3] = reaction_k[3] - rho_k * velocity_k[3];
      reaction[4] = reaction_k[4] - rho_k * velocity_k[4];
      /* DEBUG_PRINT("r-rho tilde v before projection") */
      /* DEBUG_EXPR(NV_display(reaction,5);); */
/* #ifdef DEBUG_MESSAGES */
/*       display_status_rolling_cone(projectionOnRollingCone(&reaction[0], mu_i, mu_r_i)); */
/* #else */
      status = projectionOnRollingCone(&reaction[0], mu_i, mu_r_i);
/* #endif */
      /* DEBUG_EXPR(NV_display(reaction,5);); */
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



      /* /\* DEBUG_EXPR(NV_display(velocity,5);); *\/ */

      /* The evaluation of Lipschitz constant takes into account the nonlinear part
       *  M r + q +  g(Wr+q) */
      /* double normUT_1 = sqrt(velocity[1] * velocity[1] + velocity[2] * velocity[2]); */
      /* double normOmegaT_1 = sqrt(velocity[3] * velocity[3] + velocity[4] * velocity[4]); */


      /* a1 = sqrt( */
      /*   (velocity_k[0] + mu_i * normUT + mu_r_i*normOmegaT - velocity[0] - mu_i * normUT_1 - mu_r_i*normOmegaT_1) * */
      /*   (velocity_k[0] + mu_i * normUT + mu_r_i*normOmegaT - velocity[0] - mu_i * normUT_1 - mu_r_i*normOmegaT_1) + */
      /*   (velocity_k[1] - velocity[1]) * (velocity_k[1] - velocity[1]) + */
      /*   (velocity_k[2] - velocity[2]) * (velocity_k[2] - velocity[2]) + */
      /*   (velocity_k[3] - velocity[3]) * (velocity_k[3] - velocity[3]) + */
      /*   (velocity_k[4] - velocity[4]) * (velocity_k[4] - velocity[4])); */

      /* The evaluation of Lipschitz constant is only made with the linear part
       * M r + q */
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
      DEBUG_PRINTF("a1/a2 = %12.8e\t", a1/a2);
      DEBUG_PRINTF("norm reaction = %12.8e\t", hypot5(reaction));
      DEBUG_PRINTF("norm velocity = %12.8e\n", hypot5(velocity));
      DEBUG_PRINTF("success = %i\n", success);

      ls_iter++;
    } /* end of the line search loop  */

    /* printf("--  localiter = %i\t, rho= %.10e\t, error = %.10e \n", localiter, rho, localerror); */

    /* compute local error */
    localerror =a2;
    /* rolling_fc3d_unitary_compute_and_add_error(reaction , velocity, */
    /*                                            mu_i, mu_r_i, */
    /*                                            &localerror, worktmp); */
    /* DEBUG_EXPR(NV_display(reaction,5);); */
    /* DEBUG_EXPR(NV_display(velocity,5);); */

    /*Update rho*/
      if ((rho_k*a1 < Lmin * a2) && (localerror < localerror_k))
      {
        rho =rho_k*tauinv;
      }
      else
        rho =rho_k;

      /* rho_k=1e-4; */
      /* rho=1e-4; */
      /* numerics_printf_verbose(3,"--                                                        | %i \t\t| %.10e\t| %.10e\t|%i \t\t|", localiter, rho, localerror, status); */


  } /* end of the global while loop  */

  numerics_printf_verbose(
    2,
    "--                                                        | %i \t\t| %.10e\t| %.10e\t|%i \t\t|", localiter, rho, localerror, status);
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
  options->dparam[SICONOS_DPARAM_TOL] = 1e-12;
  options->iparam[SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_IPARAM_USE_TRIVIAL_SOLUTION] =
    SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION_TRUE;

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
