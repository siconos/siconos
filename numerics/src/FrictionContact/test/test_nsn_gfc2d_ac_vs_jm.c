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

/*!\file test_nsn_gfc2d_ac_vs_jm.c
 * \brief Compare the Alart--Curnier and Jean--Moreau formulations of
 *        GFC2D_NSN_AC on the four possible contact statuses.
 *
 * Four single-contact problems are built so that the unique solution is in
 * one of the following states:
 *   - no-contact:   r_N = r_T = 0,            u_N > 0
 *   - stick:        r_N > 0, |r_T| < mu r_N, u_N = u_T = 0
 *   - slip +:       r_N > 0, r_T = +mu r_N,  u_N = 0, u_T > 0
 *   - slip -:       r_N > 0, r_T = -mu r_N,  u_N = 0, u_T < 0
 *
 * For each problem both formulations are solved and the number of iterations
 * and final relative errors are reported.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "FrictionContact_options.h"
#include "GlobalFrictionContactProblem.h"
#include "NonSmoothDrivers.h"
#include "NumericsMatrix.h"
#include "SiconosBlas.h"
#include "SolverOptions.h"
#include "gfc2d_compute_error.h"
#include "numerics_verbose.h"

/** Test scenario descriptor. */
typedef struct {
  const char* name;
  double v[2]; /**< desired global velocity */
  double r[2]; /**< desired reaction (rN, rT) */
  double mu;   /**< friction coefficient */
} Scenario;

static GlobalFrictionContactProblem* build_problem(const Scenario* s) {
  const int n = 2;   /* global DOFs */
  const int nc = 1;  /* contacts */
  const int dim = 2; /* 2D contact space */

  GlobalFrictionContactProblem* problem = globalFrictionContactProblem_new();
  problem->dimension = dim;
  problem->numberOfContacts = nc;

  /* M = I_2, column-major. */
  double* M_data = (double*)malloc((size_t)n * n * sizeof(double));
  M_data[0] = 1.0;
  M_data[1] = 0.0;
  M_data[2] = 0.0;
  M_data[3] = 1.0;
  problem->M = NM_create_with_data(NM_DENSE, n, n, M_data);

  /* H = [normal tangent] with normal=(0,1), tangent=(-1,0).
   * Local unknowns are (rN, rT).  Hence u = H^T v = (v_y, -v_x). */
  double* H_data = (double*)malloc((size_t)n * dim * sizeof(double));
  H_data[0] = 0.0;
  H_data[1] = 1.0; /* col 0 = normal */
  H_data[2] = -1.0;
  H_data[3] = 0.0; /* col 1 = tangent */
  problem->H = NM_create_with_data(NM_DENSE, n, dim, H_data);

  /* b = 0. */
  problem->b = (double*)calloc((size_t)dim * nc, sizeof(double));

  /* Reaction. */
  problem->mu = (double*)malloc((size_t)nc * sizeof(double));
  problem->mu[0] = s->mu;

  /* q = v - H*r. */
  double Hr[2];
  Hr[0] = -s->r[1]; /* (H*r)_x = -r_T */
  Hr[1] = s->r[0];  /* (H*r)_y =  r_N */
  problem->q = (double*)malloc((size_t)n * sizeof(double));
  problem->q[0] = s->v[0] - Hr[0];
  problem->q[1] = s->v[1] - Hr[1];

  problem->norm_q = cblas_dnrm2(n, problem->q, 1);
  problem->norm_b = cblas_dnrm2(dim * nc, problem->b, 1);

  return problem;
}

static int solve_scenario(const Scenario* s, int formulation, double* out_error,
                          int* out_iterations) {
  GlobalFrictionContactProblem* problem = build_problem(s);

  double* reaction = (double*)calloc((size_t)2, sizeof(double));
  double* velocity = (double*)calloc((size_t)2, sizeof(double));
  double* globalvelocity = (double*)calloc((size_t)2, sizeof(double));

  SolverOptions* options = solver_options_create(SICONOS_GLOBAL_FRICTION_2D_NSN_AC);
  if (!options) {
    fprintf(stderr, "Failed to create solver options\n");
    return 1;
  }
  options->dparam[SICONOS_DPARAM_TOL] = 1e-10;
  options->iparam[SICONOS_IPARAM_MAX_ITER] = 200;
  options->iparam[SICONOS_FRICTION_2D_NSN_FORMULATION] = formulation;
  options->iparam[SICONOS_FRICTION_2D_NSN_WARM_START_STRATEGY] =
      SICONOS_FRICTION_2D_NSN_WARM_START_YES;
  options->isSet = 1;

  int info = gfc2d_driver(problem, reaction, velocity, globalvelocity, options);

  double error = 0.0;
  int err_info = gfc2d_compute_error(problem, reaction, velocity, globalvelocity,
                                     options->dparam[SICONOS_DPARAM_TOL], options,
                                     problem->norm_q, problem->norm_b, &error);

  *out_error = error;
  *out_iterations = options->iparam[SICONOS_IPARAM_ITER_DONE];

  int ok = (info == 0) && (err_info == 0) && (error < 1e-8);

  solver_options_delete(options);
  free(reaction);
  free(velocity);
  free(globalvelocity);
  globalFrictionContact_free(problem);

  return ok ? 0 : 1;
}

int main(int argc, char** argv) {
  (void)argc;
  (void)argv;

  verbose = 0;
  setvbuf(stdout, NULL, _IONBF, 0);

  const double mu = 0.5;
  const double g = 1.0;

  Scenario scenarios[] = {
      {"no-contact", {0.0, 1.0}, {0.0, 0.0}, mu},
      {"stick", {0.0, 0.0}, {g, 0.0}, mu},
      {"slip+", {-0.7, 0.0}, {g, +mu * g}, mu},
      {"slip-", {+0.7, 0.0}, {g, -mu * g}, mu},
  };
  const int nscenarios = sizeof(scenarios) / sizeof(scenarios[0]);

  printf("%-12s %-20s %10s %12s %6s\n", "Scenario", "Formulation", "Info", "Error", "Iter");
  printf("------------------------------------------------------------------\n");

  int all_ok = 1;
  for (int i = 0; i < nscenarios; ++i) {
    for (int formulation = SICONOS_FRICTION_2D_NSN_FORMULATION_ALART_CURNIER;
         formulation <= SICONOS_FRICTION_2D_NSN_FORMULATION_JEAN_MOREAU; ++formulation) {
      double error = 0.0;
      int iterations = 0;
      int fail = solve_scenario(&scenarios[i], formulation, &error, &iterations);
      all_ok = all_ok && !fail;

      const char* name = (formulation == SICONOS_FRICTION_2D_NSN_FORMULATION_ALART_CURNIER)
                             ? "Alart-Curnier"
                             : "Jean-Moreau";
      printf("%-12s %-20s %10s %12.4e %6d\n", scenarios[i].name, name, fail ? "FAIL" : "OK",
             error, iterations);
    }
  }

  printf("------------------------------------------------------------------\n");
  if (all_ok) {
    printf("All scenarios converged with both formulations.\n");
    return 0;
  } else {
    printf("At least one scenario failed to converge.\n");
    return 1;
  }
}
