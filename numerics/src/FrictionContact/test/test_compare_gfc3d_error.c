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

#include <stdio.h>   // for printf
#include <stdlib.h>  // for calloc, free

#include "FrictionContact_options.h"
#include "GlobalFrictionContactProblem.h"  // for GlobalFrictionContactProblem
#include "NonSmoothDrivers.h"              // for gfc3d_driver
#include "SiconosBlas.h"                   // for cblas_dnrm2
#include "SolverOptions.h"                 // for solver_options_create, ...
#include "gfc3d_compute_error.h"           // for gfc3d_compute_error, ...

typedef struct {
  const char* label;
  int max_iter;
  double tol;
} TestCaseConfig;

typedef struct {
  int norm;
  int scope;
  int reduction;
  const char* label;
} ScalingCase;

/**
 * @brief Compare the GFC3D complementarity error measured with the default
 *        error computation and with each available scaling option on a
 *        selection of problems solved with NSGS.
 *
 *        The repository only contains a handful of plain-text GFC3D data
 *        files, so the test cases are obtained by running a few NSGS presets
 *        on each available problem.
 */
int main(void) {
  const char* filenames[] = {"./data/GFC3D_Example0.dat",
                             "./data/GFC3D_Example00.dat",
                             "./data/GFC3D_Example00_badly_scaled.dat",
                             "./data/GFC3D_Example1.dat",
                             "./data/GFC3D_OneContact.dat",
                             "./data/GFC3D_TwoRods1.dat",
                             "./data/GFC3D_Example0_SBM.dat"};

  const int n_files = (int)(sizeof(filenames) / sizeof(filenames[0]));

  /* A few NSGS presets to give around 20 test cases. */
  TestCaseConfig presets[] = {
      {"default", 100, 1e-3}, {"tight", 200, 1e-5}, {"fast", 50, 1e-2}};
  const int n_presets = (int)(sizeof(presets) / sizeof(presets[0]));

  ScalingCase scalings[] = {
      {SICONOS_FRICTION_3D_SCALING_NO, SICONOS_FRICTION_3D_SCALING_SCOPE_PER_CONE,
       SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX, "no"},
      {SICONOS_FRICTION_3D_SCALING_EUCLIDEAN, SICONOS_FRICTION_3D_SCALING_SCOPE_PER_CONE,
       SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX, "euc_pc_max"},
      {SICONOS_FRICTION_3D_SCALING_EUCLIDEAN, SICONOS_FRICTION_3D_SCALING_SCOPE_PER_VECTOR,
       SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX, "euc_pv_max"},
      {SICONOS_FRICTION_3D_SCALING_CONE, SICONOS_FRICTION_3D_SCALING_SCOPE_PER_CONE,
       SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX, "cone_pc_max"},
      {SICONOS_FRICTION_3D_SCALING_CONE_MU, SICONOS_FRICTION_3D_SCALING_SCOPE_PER_CONE,
       SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX, "cmu_pc_max"},
      {SICONOS_FRICTION_3D_SCALING_CONE_MU, SICONOS_FRICTION_3D_SCALING_SCOPE_PER_VECTOR,
       SICONOS_FRICTION_3D_SCALING_REDUCTION_MEAN, "cmu_pv_mean"},
      {SICONOS_FRICTION_3D_SCALING_CONE_WEIGHTED, SICONOS_FRICTION_3D_SCALING_SCOPE_PER_CONE,
       SICONOS_FRICTION_3D_SCALING_REDUCTION_MIXED, "wgt_pc_mix"},
      {SICONOS_FRICTION_3D_SCALING_CONE_WEIGHTED, SICONOS_FRICTION_3D_SCALING_SCOPE_PER_VECTOR,
       SICONOS_FRICTION_3D_SCALING_REDUCTION_MAX, "wgt_pv_max"}};
  const int n_scalings = (int)(sizeof(scalings) / sizeof(scalings[0])); /* used for loops */
#define N_SCALINGS_LOCAL ((int)(sizeof(scalings) / sizeof(scalings[0])))

  printf(
      "---------------------------------------------------------------------------------------"
      "-----------"
      "------------------------------------------------\n");
  printf("| %-46s | info | iter |   standard   ", "problem (preset)");
  for (int s = 0; s < n_scalings; ++s) {
    printf("| %10s   ", scalings[s].label);
  }
  printf("|\n");
  printf(
      "---------------------------------------------------------------------------------------"
      "-----------"
      "------------------------------------------------\n");

  for (int i = 0; i < n_files; ++i) {
    const char* filename = filenames[i];
    GlobalFrictionContactProblem* problem = globalFrictionContact_new_from_filename(filename);

    if (!problem || !problem->M || !problem->H || !problem->q || !problem->b || !problem->mu) {
      printf("| %-46s |     N/A     |  read failed  |", filename);
      for (int s = 0; s < n_scalings; ++s) printf("               |");
      printf("\n");
      if (problem) globalFrictionContact_free(problem);
      continue;
    }

    const int dim = problem->dimension;
    if (dim != 3) {
      printf("| %-46s |     N/A     |  dim != 3     |", filename);
      for (int s = 0; s < n_scalings; ++s) printf("               |");
      printf("\n");
      globalFrictionContact_free(problem);
      continue;
    }

    const int nc = problem->numberOfContacts;
    const int n = (int)problem->M->size0;
    const int m = dim * nc;

    for (int p = 0; p < n_presets; ++p) {
      double* reaction = (double*)calloc(m, sizeof(double));
      double* velocity = (double*)calloc(m, sizeof(double));
      double* globalvelocity = (double*)calloc(n, sizeof(double));

      SolverOptions* options = solver_options_create(SICONOS_GLOBAL_FRICTION_3D_NSGS);
      options->dparam[SICONOS_DPARAM_TOL] = presets[p].tol;
      options->iparam[SICONOS_IPARAM_MAX_ITER] = presets[p].max_iter;

      const int info = gfc3d_driver(problem, reaction, velocity, globalvelocity, options);
      const int iter = options->iparam[SICONOS_IPARAM_ITER_DONE];

      const double norm_q = cblas_dnrm2(n, problem->q, 1);
      const double norm_b = cblas_dnrm2(m, problem->b, 1);

      /* Standard error: use the default scaling (the option is left untouched). */
      double error_standard = 0.0;
      gfc3d_compute_error(problem, reaction, velocity, globalvelocity,
                          /* tolerance = */ 1.0, options, norm_q, norm_b, &error_standard);

      /* Error with each explicit scaling choice. */
      double errors[N_SCALINGS_LOCAL];
      for (int s = 0; s < n_scalings; ++s) {
        options->iparam[SICONOS_FRICTION_3D_SCALING_NORM] = scalings[s].norm;
        options->iparam[SICONOS_FRICTION_3D_SCALING_SCOPE] = scalings[s].scope;
        options->iparam[SICONOS_FRICTION_3D_SCALING_REDUCTION] = scalings[s].reduction;
        options->dparam[SICONOS_FRICTION_3D_SCALING_ALPHA] = 1.0;
        errors[s] = 0.0;
        gfc3d_compute_error(problem, reaction, velocity, globalvelocity,
                            /* tolerance = */ 1.0, options, norm_q, norm_b, &errors[s]);
      }

      char label[128];
      snprintf(label, sizeof(label), "%s (%s)", filename, presets[p].label);

      printf("| %-46s |  %3d |  %4d | %12.4e ", label, info, iter, error_standard);
      for (int s = 0; s < n_scalings; ++s) {
        printf("| %12.4e ", errors[s]);
      }
      printf("|\n");
      fflush(stdout);

      solver_options_delete(options);
      free(reaction);
      free(velocity);
      free(globalvelocity);
    }

    globalFrictionContact_free(problem);
  }

  printf(
      "---------------------------------------------------------------------------------------"
      "-----------"
      "------------------------------------------------\n");

  return 0;
}
