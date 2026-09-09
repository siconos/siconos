/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


#include "CohesiveFrictionContactProblem.h"
#include "CohesiveFrictionContact_options.h"
#include "SolverOptions.h"
#include "NumericsMatrix.h"
#include "cohesive_friction_3d_projection.h"
#include "cohesive_friction_3d_compute_error.h"

#include "numerics_verbose.h"
#include "numerics_errors.h"
#include "projectionOnDisk.h"

#include "solver_registry.h"

/* #define DEBUG_NOCOLOR */
/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "siconos_debug.h"  // for DEBUG_PRINTF, DEBUG_EXPR, DEBU...
#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"
#endif
#include "NumericsVector.h"

void cohesive_friction_3d_projection_initialize(CohesiveFrictionContactProblem* main_problem,
						SolverOptions* localsolver_options) {
  size_t n_coh = main_problem->numberOfCohesivePoints;  
  /* printf("fc3d_projectionOnConeWithLocalIteration_initialize. Allocation of dwork\n"); */
  if (!localsolver_options->dWork || localsolver_options->dWorkSize < n_coh) {
    localsolver_options->dWork =
      (double*)realloc(localsolver_options->dWork, n_coh * sizeof(double));
    localsolver_options->dWorkSize = n_coh;
  }
  for (size_t i = 0; i < n_coh; i++) {
    localsolver_options->dWork[i] = 1.0;
  }  
}

void cohesive_friction_3d_projection_free(CohesiveFrictionContactProblem* main_problem, SolverOptions * options) {}


int cohesive_friction_3d_projection_solve(CohesiveFrictionContactProblem* localproblem, double* reaction,
					  SolverOptions* options) {



  DEBUG_BEGIN("cohesive_friction_3d_projection_solve(...)\n");

  DEBUG_EXPR(cohesiveFrictionContact_display(localproblem););
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;

  double* MLocal = localproblem->M->matrix0;
  double* qLocal = localproblem->q;
 
  double c_n = localproblem->c_n[0];
  double c_t = localproblem->c_t[0];

  
  reaction[0] = - c_n;
  /* reaction[1] = 0.0; */
  /* reaction[2] = 0.0; */
  /* int nLocal = 3; */


  int cohesive_idx = options->iparam[SICONOS_COHESIVE_FRICTION_IPARAM_CURRENT_CONTACT_NUMBER];                     

  
  double rho = options->dWork[cohesive_idx],
         rho_k;
  DEBUG_PRINTF(" Contact options->iparam[SICONOS_COHESIVE_FRICTION_IPARAM_CURRENT_CONTACT_NUMBER] = %i\n",
              cohesive_idx);
  DEBUG_PRINTF("saved rho = %14.7e\n", rho);
  assert(rho > 0);

  /* int incx = 1, incy = 1; */
  /* int i; */

  double velocity[3], velocity_k[3], reaction_k[3], worktmp[3];
  double localerror = 1.0;
  // printf ("localerror = %14.7e\n",localerror );
  int localiter = 0;
  double localtolerance = dparam[SICONOS_DPARAM_TOL];

  /* Variable for Line_search */
  double a1, a2;
  int success = 0;
  double localerror_k;
  int ls_iter = 0;
  int ls_itermax = 10;

  double tau = 2.0 / 3.0, tauinv = 3.0 / 2.0, L = 0.9, Lmin = 0.3;

  numerics_printf_verbose(2, "--  cohesive_friction_3d_projectionOnConeWithLocalIteration_solve contact = %i",
                          cohesive_idx);
  numerics_printf_verbose(2,
                          "--  cohesive_friction_3d_projectionOnConeWithLocalIteration_solve | localiter \t| "
                          "rho \t\t\t| error\t\t\t|");
  numerics_printf_verbose(
      2,                  "--                                                                | %i \t\t| %.10e\t| %.10e\t|",
      localiter, rho, localerror);

  /* printf ("localtolerance = %14.7e\n",localtolerance ); */
  /* printf("iparam[SICONOS_IPARAM_MAX_ITER] %i \n", iparam[SICONOS_IPARAM_MAX_ITER]); */
  while ((localerror > localtolerance) && (localiter < iparam[SICONOS_IPARAM_MAX_ITER] )) {
    DEBUG_PRINTF("\n Local iteration starts % i \n", localiter);
    localiter++;

    /*    printf ("reaction[0] = %14.7e\n",reaction[0]); */
    /*    printf ("reaction[1] = %14.7e\n",reaction[1]); */
    /*    printf ("reaction[2] = %14.7e\n",reaction[2]); */

    /* Store the error */
    localerror_k = localerror;

    /* store the reaction at the beginning of the iteration */
    /* cblas_dcopy(nLocal , reaction , 1 , reaction_k, 1); */

    reaction_k[0] = reaction[0];
    reaction_k[1] = reaction[1];
    reaction_k[2] = reaction[2];
    DEBUG_EXPR(NV_display(reaction_k, 3););
    /* /\* velocity_k <- q  *\/ */
    /* cblas_dcopy_msan(nLocal , qLocal , 1 , velocity_k, 1); */
    /* /\* velocity_k <- q + M * reaction  *\/ */
    /* cblas_dgemv(CblasColMajor,CblasNoTrans, nLocal, nLocal, 1.0, MLocal, 3, reaction,
     * incx, 1.0, velocity_k, incy); */
    for (int i = 0; i < 3; i++)
      velocity_k[i] = MLocal[i + 0 * 3] * reaction[0] + qLocal[i] +
                      MLocal[i + 1 * 3] * reaction[1] + MLocal[i + 2 * 3] * reaction[2];

    ls_iter = 0;
    success = 0;
    rho_k = rho / tau;

    while (!success && (ls_iter < ls_itermax)) {
      rho_k = rho_k * tau;
      DEBUG_PRINTF("rho_k =%f\n", rho_k);
      DEBUG_EXPR(NV_display(velocity_k, 3););      
      reaction[0] = reaction_k[0] - rho_k * velocity_k[0];
      reaction[1] = reaction_k[1] - rho_k * velocity_k[1];
      reaction[2] = reaction_k[2] - rho_k * velocity_k[2];
      DEBUG_PRINT("r-rho tilde v before projection:  \n");
      DEBUG_EXPR(NV_display(reaction, 3););

      reaction[0] = -c_n;
      //reaction[1] = -qLocal[1];
      //reaction[2] = 0;
      projectionOnDisk(&reaction[1], c_t);
      
      DEBUG_PRINT("reaction after projection:  \n");
      DEBUG_EXPR(NV_display(reaction, 3););

      /* velocity <- q  */
      /* cblas_dcopy(nLocal , qLocal , 1 , velocity, 1); */
      /* velocity <- q + M * reaction  */
      /* cblas_dgemv(CblasColMajor,CblasNoTrans, nLocal, nLocal, 1.0, MLocal, 3, reaction,
       * incx, 1.0, velocity, incy); */

      for (int i = 0; i < 3; i++)
        velocity[i] = MLocal[i + 0 * 3] * reaction[0] + qLocal[i]
	  + MLocal[i + 1 * 3] * reaction[1]
	  + MLocal[i + 2 * 3] * reaction[2];

      a1 = sqrt((velocity_k[0] - velocity[0]) * (velocity_k[0] - velocity[0]) +
                (velocity_k[1] - velocity[1]) * (velocity_k[1] - velocity[1]) +
                (velocity_k[2] - velocity[2]) * (velocity_k[2] - velocity[2]));

      a2 = sqrt((reaction_k[0] - reaction[0]) * (reaction_k[0] - reaction[0]) +
                (reaction_k[1] - reaction[1]) * (reaction_k[1] - reaction[1]) +
                (reaction_k[2] - reaction[2]) * (reaction_k[2] - reaction[2]));

      success = (rho_k * a1 <= L * a2) ? 1 : 0;

      DEBUG_PRINTF("rho_k = %12.8e\t", rho_k);
      DEBUG_PRINTF("a1 = %12.8e\t", a1);
      DEBUG_PRINTF("a2 = %12.8e\t", a2);
      DEBUG_PRINTF("norm reaction = %12.8e\t",
                   sqrt(reaction[0] * reaction[0] + reaction[1] * reaction[1] +
                        reaction[2] * reaction[2]));
      DEBUG_PRINTF("success = %i\n", success);

      ls_iter++;
    }

    /* printf("--  localiter = %i\t, rho= %.10e\t, error = %.10e \n", localiter, rho,
     * localerror); */

    /* compute local error */
    localerror = 0.0;
    cohesive_friction_3d_unitary_compute_and_add_error(reaction, velocity, c_t, &localerror,
                                                       worktmp);
    localerror = sqrt(localerror);
    DEBUG_PRINTF("localerror = %e\n", localerror)
    ;
    /*Update rho*/
    if ((rho_k * a1 < Lmin * a2) && (localerror < localerror_k)) {
      rho = rho_k * tauinv;
    } else
      rho = rho_k;

    numerics_printf_verbose(
      2,                  "--                                                                | %i \t\t| %.10e\t| %.10e\t|",
      localiter, rho, localerror);
  }
  options->dWork[cohesive_idx] = rho;
  options->dparam[SICONOS_DPARAM_RESIDU] = localerror;
  DEBUG_PRINTF("final rho  =%e\n", rho);
  /* printf("velocity:");NV_display(velocity,3); */
  /* printf("reaction:");NV_display(reaction,3); */
  /* printf("norm reaction: %e\n", sqrt(reaction[1]*reaction[1] + reaction[2]* reaction[2])); */
   
  DEBUG_END("cohesive_friction_3d_projectionOnConeWithLocalIteration_solve(...)\n");
  if (localerror > localtolerance) return 1;


  
 
  return 0;
}

/* Projection on Cone */
static void cohesive_friction_3d_proj_set_default(SolverOptions* options) { /* No specific defaults */ }

static int cohesive_friction_3d_proj_init_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  (void)options;
  return NUMERICS_OK;
}

static int cohesive_friction_3d_proj_solve_wrap(void* problem, double* reaction, double* velocity,
                                     SolverOptions* options) {
  (void)velocity;
  return cohesive_friction_3d_projection_solve((CohesiveFrictionContactProblem*)problem, reaction, options);
}

REGISTER_SOLVER(SICONOS_COHESIVE_FRICTION_3D_PROJECTION,
                "SICONOS_COHESIVE_FRICTION_3D_PROJECTION", "Projection solver (local solver) for cohesion",
                cohesive_friction_3d_proj_init_wrap,
                cohesive_friction_3d_proj_solve_wrap, NULL, NULL,
                cohesive_friction_3d_proj_set_default,
                100, 1e-14, 1)
