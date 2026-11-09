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

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "FrictionContactProblem.h"
#include "FrictionContact_options.h"
#include "GlobalFrictionContactProblem.h"
#include "NonSmoothDrivers.h"
#include "SolverOptions.h"
#include "numerics_verbose.h"

int main(int argc, char** argv) {
  (void)argc;
  (void)argv;

  verbose = 0;

  const char* filename = "./data/gfc2d_granularflowonwall_00230_00206.dat";

  GlobalFrictionContactProblem* problem = globalFrictionContact_new_from_filename(filename);
  if (!problem) {
    fprintf(stderr, "Failed to load %s\n", filename);
    return 1;
  }

  const int dim = problem->dimension;
  const int nc = problem->numberOfContacts;
  const int m = dim * nc;
  const int n_global = problem->M->size0;

  printf("Loaded %s\n", filename);
  printf("  dimension = %d, contacts = %d, local dofs = %d, global dofs = %d\n", dim, nc, m,
         n_global);

  /* Build an initial guess by solving the reduced problem with FC2D_NSGS. */
  FrictionContactProblem* rproblem = globalFrictionContact_reduce_to_FrictionContact(problem);
  if (!rproblem) {
    fprintf(stderr, "Failed to reduce global problem to FrictionContactProblem\n");
    return 1;
  }

  double* reaction0 = (double*)calloc((size_t)m, sizeof(double));
  double* velocity0 = (double*)calloc((size_t)m, sizeof(double));

  SolverOptions* options_nsgs = solver_options_create(SICONOS_FRICTION_2D_NSGS);
  options_nsgs->dparam[SICONOS_DPARAM_TOL] = 1e-4;
  options_nsgs->iparam[SICONOS_IPARAM_MAX_ITER] = 100;

  int info_nsgs = fc2d_driver(rproblem, reaction0, velocity0, options_nsgs);
  printf("\nFC2D_NSGS initial guess: info = %d, iterations = %d, residual = %e\n", info_nsgs,
         options_nsgs->iparam[SICONOS_IPARAM_ITER_DONE],
         options_nsgs->dparam[SICONOS_DPARAM_RESIDU]);

  if (info_nsgs != 0) {
    fprintf(stderr, "Initial guess computation failed\n");
    return 1;
  }

  double* globalvelocity0 = (double*)calloc((size_t)n_global, sizeof(double));
  globalFrictionContact_computeGlobalVelocity(problem, reaction0, globalvelocity0);

  /* Solve with constant rho (warm-start NO) using the NSGS initial guess. */
  GlobalFrictionContactProblem* problem_no = globalFrictionContact_copy(problem);
  double* reaction_no = (double*)calloc((size_t)m, sizeof(double));
  double* velocity_no = (double*)calloc((size_t)m, sizeof(double));
  double* globalvelocity_no = (double*)calloc((size_t)n_global, sizeof(double));
  memcpy(reaction_no, reaction0, (size_t)m * sizeof(double));
  memcpy(velocity_no, velocity0, (size_t)m * sizeof(double));
  memcpy(globalvelocity_no, globalvelocity0, (size_t)n_global * sizeof(double));

  SolverOptions* options_no = solver_options_create(SICONOS_GLOBAL_FRICTION_2D_NSN_AC);
  options_no->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  options_no->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
  options_no->iparam[SICONOS_FRICTION_2D_NSN_WARM_START_STRATEGY] =
      SICONOS_FRICTION_2D_NSN_WARM_START_NO;
  options_no->iparam[SICONOS_FRICTION_2D_NSN_REGUL] =
      SICONOS_FRICTION_2D_NSN_REGULARIZATION_YES;
  options_no->dparam[SICONOS_FRICTION_2D_NSN_REGUL] = 1e-6;

  int info_no =
      gfc2d_driver(problem_no, reaction_no, velocity_no, globalvelocity_no, options_no);
  printf("\nGFC2D_NSN_AC (constant rho): info = %d, iterations = %d, residual = %e\n", info_no,
         options_no->iparam[SICONOS_IPARAM_ITER_DONE],
         options_no->dparam[SICONOS_DPARAM_RESIDU]);

  /* Solve with rho from the diagonal approximation using the same initial guess. */
  GlobalFrictionContactProblem* problem_diag = globalFrictionContact_copy(problem);
  double* reaction_diag = (double*)calloc((size_t)m, sizeof(double));
  double* velocity_diag = (double*)calloc((size_t)m, sizeof(double));
  double* globalvelocity_diag = (double*)calloc((size_t)n_global, sizeof(double));
  memcpy(reaction_diag, reaction0, (size_t)m * sizeof(double));
  memcpy(velocity_diag, velocity0, (size_t)m * sizeof(double));
  memcpy(globalvelocity_diag, globalvelocity0, (size_t)n_global * sizeof(double));

  SolverOptions* options_diag = solver_options_create(SICONOS_GLOBAL_FRICTION_2D_NSN_AC);
  options_diag->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  options_diag->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
  options_diag->iparam[SICONOS_FRICTION_2D_NSN_WARM_START_STRATEGY] =
      SICONOS_FRICTION_2D_NSN_WARM_START_RHO_FROM_DIAGONAL;
  options_diag->iparam[SICONOS_FRICTION_2D_NSN_REGUL] =
      SICONOS_FRICTION_2D_NSN_REGULARIZATION_YES;
  options_diag->dparam[SICONOS_FRICTION_2D_NSN_REGUL] = 1e-6;

  int info_diag = gfc2d_driver(problem_diag, reaction_diag, velocity_diag, globalvelocity_diag,
                               options_diag);
  printf("GFC2D_NSN_AC (rho-from-diagonal): info = %d, iterations = %d, residual = %e\n",
         info_diag, options_diag->iparam[SICONOS_IPARAM_ITER_DONE],
         options_diag->dparam[SICONOS_DPARAM_RESIDU]);

  int failures = 0;
  if (info_no != 0) {
    fprintf(stderr, "NSN with constant rho failed to converge (info = %d)\n", info_no);
    failures++;
  }
  if (info_diag != 0) {
    fprintf(stderr, "NSN with rho-from-diagonal failed to converge (info = %d)\n", info_diag);
    failures++;
  }

  if (info_no == 0 && info_diag == 0) {
    const double abs_tol = 1e-6;
    const double rel_tol = 1e-5;
    double max_diff = 0.0;
    int worst = -1;

    for (int i = 0; i < m; ++i) {
      double diff = fabs(velocity_no[i] - velocity_diag[i]);
      double scale = fmax(fabs(velocity_no[i]), fabs(velocity_diag[i]));
      double tol = fmax(abs_tol, rel_tol * scale);
      if (diff > max_diff) {
        max_diff = diff;
        worst = i;
      }
      if (diff > tol) {
        if (failures < 10) {
          printf(
              "  velocity mismatch at %d: constant = %e, diagonal = %e, diff = %e (tol = "
              "%e)\n",
              i, velocity_no[i], velocity_diag[i], diff, tol);
        }
        failures++;
      }
    }
    printf("Max velocity difference = %e at component %d\n", max_diff, worst);
  }

  solver_options_delete(options_nsgs);
  solver_options_delete(options_no);
  solver_options_delete(options_diag);
  free(reaction0);
  free(velocity0);
  free(globalvelocity0);
  free(reaction_no);
  free(velocity_no);
  free(globalvelocity_no);
  free(reaction_diag);
  free(velocity_diag);
  free(globalvelocity_diag);
  frictionContactProblem_free(rproblem);
  globalFrictionContactProblem_free(problem_no);
  globalFrictionContactProblem_free(problem_diag);
  globalFrictionContactProblem_free(problem);

  return (failures == 0) ? 0 : 1;
}
