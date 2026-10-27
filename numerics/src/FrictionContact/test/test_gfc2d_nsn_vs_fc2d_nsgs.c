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

  /* Use the largest available gfc2d granular-flow problem. */
  const char* filename = "./data/gfc2d_granularflowonwall_00230_00206.dat";

  GlobalFrictionContactProblem* gproblem = globalFrictionContact_new_from_filename(filename);
  if (!gproblem) {
    fprintf(stderr, "Failed to load %s\n", filename);
    return 1;
  }

  const int dim = gproblem->dimension;
  const int nc = gproblem->numberOfContacts;
  const int m = dim * nc;
  const int n_global = gproblem->M->size0;

  printf("Loaded %s\n", filename);
  printf("  dimension = %d, contacts = %d, local dofs = %d, global dofs = %d\n", dim, nc, m,
         n_global);
  printf("  M storage = %d, H storage = %d\n", gproblem->M->storageType,
         gproblem->H->storageType);

  /* Solve the global problem with GFC2D_NSN_AC. */
  double* reaction_g = (double*)calloc((size_t)m, sizeof(double));
  double* velocity_g = (double*)calloc((size_t)m, sizeof(double));
  double* globalvelocity_g = (double*)calloc((size_t)n_global, sizeof(double));

  SolverOptions* options_nsn = solver_options_create(SICONOS_GLOBAL_FRICTION_2D_NSN_AC);
  options_nsn->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  options_nsn->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;

  int info_nsn = gfc2d_driver(gproblem, reaction_g, velocity_g, globalvelocity_g, options_nsn);
  printf("\nGFC2D_NSN_AC: info = %d, iterations = %d, residual = %e\n", info_nsn,
         options_nsn->iparam[SICONOS_IPARAM_ITER_DONE],
         options_nsn->dparam[SICONOS_DPARAM_RESIDU]);

  /* Reduce the same global problem to a local FrictionContactProblem. */
  FrictionContactProblem* rproblem = globalFrictionContact_reduce_to_FrictionContact(gproblem);
  if (!rproblem) {
    fprintf(stderr, "Failed to reduce global problem to FrictionContactProblem\n");
    return 1;
  }
  printf("\nReduced problem: local dofs = %d, M storage = %d\n",
         rproblem->dimension * rproblem->numberOfContacts, rproblem->M->storageType);

  /* Solve the reduced problem with fc2d_nsgs. */
  double* reaction_r = (double*)calloc((size_t)m, sizeof(double));
  double* velocity_r = (double*)calloc((size_t)m, sizeof(double));

  SolverOptions* options_nsgs = solver_options_create(SICONOS_FRICTION_2D_NSGS);
  options_nsgs->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  options_nsgs->iparam[SICONOS_IPARAM_MAX_ITER] = 100000;

  int info_nsgs = fc2d_driver(rproblem, reaction_r, velocity_r, options_nsgs);
  printf("FC2D_NSGS:    info = %d, iterations = %d, residual = %e\n", info_nsgs,
         options_nsgs->iparam[SICONOS_IPARAM_ITER_DONE],
         options_nsgs->dparam[SICONOS_DPARAM_RESIDU]);

  /* Compare local velocities. */
  const double abs_tol = 1e-6;
  const double rel_tol = 1e-5;
  int failures = 0;
  double max_diff = 0.0;
  double sum_diff = 0.0;
  int worst = -1;

  for (int i = 0; i < m; ++i) {
    double diff = fabs(velocity_g[i] - velocity_r[i]);
    double scale = fmax(fabs(velocity_g[i]), fabs(velocity_r[i]));
    double tol = fmax(abs_tol, rel_tol * scale);
    sum_diff += diff;
    if (diff > max_diff) {
      max_diff = diff;
      worst = i;
    }
    if (diff > tol) {
      if (failures < 10) {
        printf("  velocity mismatch at %d: gfc2d = %e, fc2d = %e, diff = %e (tol = %e)\n", i,
               velocity_g[i], velocity_r[i], diff, tol);
      }
      failures++;
    }
  }

  printf("\nVelocity comparison:\n");
  printf("  max absolute difference = %e (at index %d)\n", max_diff, worst);
  printf("  average absolute difference = %e\n", sum_diff / (double)m);
  printf("  components above threshold = %d / %d\n", failures, m);

  solver_options_delete(options_nsn);
  solver_options_delete(options_nsgs);
  free(reaction_g);
  free(velocity_g);
  free(globalvelocity_g);
  free(reaction_r);
  free(velocity_r);
  frictionContactProblem_free(rproblem);
  globalFrictionContactProblem_free(gproblem);

  if (info_nsn != 0 || info_nsgs != 0) {
    printf("\nFAIL: at least one solver did not converge.\n");
    return 1;
  }
  if (failures > 0) {
    printf("\nFAIL: velocity differences exceed threshold.\n");
    return 1;
  }

  printf("\nPASS: velocities agree within threshold.\n");
  return 0;
}
