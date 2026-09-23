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

/*!\file test_vi_fpp_gfc2d_single_particle.c
 * \brief Minimal sanity test for GFC2D_VI_FPP.
 *
 * A single 2D particle with identity mass matrix is pressed against a
 * horizontal wall by a vertical gravity force.  There is one contact, with
 * normal pointing upward and tangent pointing rightward.  The expected
 * solution is zero global/local velocity and a normal reaction equal to the
 * gravity magnitude.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "FrictionContact_options.h"
#include "GlobalFrictionContactProblem.h"
#include "NonSmoothDrivers.h"
#include "NumericsMatrix.h"
#include "SiconosBlas.h"
#include "SolverOptions.h"
#include "gfc2d_compute_error.h"
#include "numerics_verbose.h"

static int nearly_equal(double a, double b, double tol) { return fabs(a - b) <= tol; }

int main(int argc, char** argv) {
  (void)argc;
  (void)argv;

  /* Problem sizes */
  const int n = 2;    /* global degrees of freedom */
  const int nc = 1;   /* number of contacts */
  const int dim = 2;  /* 2D contact space */

  /* Physical data */
  const double g = 1.0;
  const double mu = 1.0;

  verbose = 1;
  setvbuf(stdout, NULL, _IONBF, 0);

  GlobalFrictionContactProblem* problem = globalFrictionContactProblem_new();
  problem->dimension = dim;
  problem->numberOfContacts = nc;

  /* M = I_2, stored column-major */
  double* M_data = (double*)malloc((size_t)n * n * sizeof(double));
  M_data[0] = 1.0;
  M_data[1] = 0.0; /* column 0 */
  M_data[2] = 0.0;
  M_data[3] = 1.0; /* column 1 */
  problem->M = NM_create_from_data(NM_DENSE, n, n, M_data);

  /* H = [n t] with n = (0,1) (upward normal) and t = (-1,0) (leftward tangent) */
  double* H_data = (double*)malloc((size_t)n * dim * sizeof(double));
  H_data[0] = 0.0;
  H_data[1] = 1.0; /* column 0 = normal */
  H_data[2] = -1.0;
  H_data[3] = 0.0; /* column 1 = tangent */
  problem->H = NM_create_from_data(NM_DENSE, n, dim, H_data);

  /* q = (0.8, -g): gravity pulls the particle downward */
  problem->q = (double*)malloc((size_t)n * sizeof(double));
  problem->q[0] = 0.8;
  problem->q[1] = -g;

  /* b = 0: zero initial gap/relative velocity */
  problem->b = (double*)calloc((size_t)dim * nc, sizeof(double));

  problem->mu = (double*)malloc((size_t)nc * sizeof(double));
  problem->mu[0] = mu;

  problem->norm_q = cblas_dnrm2(n, problem->q, 1);
  problem->norm_b = cblas_dnrm2(dim * nc, problem->b, 1);

  /* Solution vectors */
  double* reaction = (double*)calloc((size_t)dim * nc, sizeof(double));
  double* velocity = (double*)calloc((size_t)dim * nc, sizeof(double));
  double* globalvelocity = (double*)calloc((size_t)n, sizeof(double));

  /* Solver options */
  SolverOptions* options = solver_options_create(SICONOS_GLOBAL_FRICTION_2D_VI_FPP);
  if (!options) {
    fprintf(stderr, "Failed to create solver options\n");
    return 1;
  }
  options->dparam[SICONOS_DPARAM_TOL] = 1e-8;
  options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
  options->isSet = 1;
  
  /* Solve */
  int info = gfc2d_driver(problem, reaction, velocity, globalvelocity, options);

  /* Compute the relative error */
  double error = 0.0;
  int err_info = gfc2d_compute_error(problem, reaction, velocity, globalvelocity,
                                     options->dparam[SICONOS_DPARAM_TOL], options,
                                     problem->norm_q, problem->norm_b, &error);

  printf("Solver info          = %d\n", info);
  printf("Compute error info   = %d\n", err_info);
  printf("Relative error       = %e\n", error);
  printf("Global velocity      = [% e, % e]\n", globalvelocity[0], globalvelocity[1]);
  printf("Reaction             = [% e, % e]\n", reaction[0], reaction[1]);
  printf("Local velocity       = [% e, % e]\n", velocity[0], velocity[1]);
  printf("Iterations done      = %d\n", options->iparam[SICONOS_IPARAM_ITER_DONE]);

  int test_ok = 1;
  if (info != 0) test_ok = 0;
  if (err_info != 0) test_ok = 0;
  if (error > 1e-6) test_ok = 0;
  if (!nearly_equal(globalvelocity[0], 0.0, 1e-6)) test_ok = 0;
  if (!nearly_equal(globalvelocity[1], 0.0, 1e-6)) test_ok = 0;
  if (!nearly_equal(reaction[0], g, 1e-6)) test_ok = 0;
  if (!nearly_equal(reaction[1], 0.8, 1e-6)) test_ok = 0;

  if (test_ok) {
    printf("test: success\n");
  } else {
    printf("test: failure\n");
  }

  solver_options_delete(options);
  free(reaction);
  free(velocity);
  free(globalvelocity);
  globalFrictionContact_free(problem);

  return test_ok ? 0 : 1;
}
