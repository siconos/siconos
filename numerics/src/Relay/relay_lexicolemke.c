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
#include <stdlib.h>  // for malloc, free

#include "LCP_Solvers.h"                   // for lcp_compute_error
#include "LinearComplementarityProblem.h"  // for LinearComplementarityProblem
#include "NonSmoothDrivers.h"              // for linearComplementarity_driver
#include "NumericsFwd.h"                   // for LinearComplementarityProblem
#include "RelayProblem.h"                  // for RelayProblem
#include "Relay_Solvers.h"                 // for relay_to_lcp, relay_lexico...
#include "relay_cst.h"
#include "SolverOptions.h"                 // for SolverOptions, SICONOS_DPA...
#include "lcp_cst.h"                       // for SICONOS_LCP_LEMKE
#include "relay_cst.h"                     // for SICONOS_LCP_LEMKE

/* Solver registration system */
#include "utils/solver_registry.h"
#include "utils/numerics_errors.h"
void relay_lexicolemke(RelayProblem *problem, double *z, double *w, int *info,
                       SolverOptions *options) {
  int i;
  // conversion into LCP
  LinearComplementarityProblem *lcp_problem =
      (LinearComplementarityProblem *)malloc(sizeof(LinearComplementarityProblem));

  /* Relay_display(problem); */

  relay_to_lcp(problem, lcp_problem);

  /* linearComplementarity_display(lcp_problem); */

  double *zlcp = (double *)calloc(lcp_problem->size, sizeof(double));
  double *wlcp = (double *)calloc(lcp_problem->size, sizeof(double));

  /*  FILE * fcheck = fopen("lcp_relay.dat","w"); */
  /*  info = linearComplementarity_printInFile(lcp_problem,fcheck); */

  // Call the lcp_solver
  options->solverId = SICONOS_LCP_LEMKE;
  *info = linearComplementarity_driver(lcp_problem, zlcp, wlcp, options);
  if (options->filterOn > 0)
    lcp_compute_error(lcp_problem, zlcp, wlcp, options->dparam[SICONOS_DPARAM_TOL],
                      &(options->dparam[SICONOS_DPARAM_RESIDU]));

  // Conversion of result
  for (i = 0; i < problem->size; i++) {
    /* z[i] = 1.0/2.0*(zlcp[i]-wlcp[i+problem->size]); works only for ub=1 and lb=-1 */
    z[i] = zlcp[i] + problem->lb[i];
    /* printf("problem->ub[ %i]=%12.10e\n", i, problem->ub[i]); */
    /* printf("zlcp[ %i]=%12.10e\n", i, zlcp[i]); */
    /* printf("z[ %i]=%12.10e\n", i, z[i]); */

    w[i] = wlcp[i] - zlcp[i + problem->size];
    // printf("w[ %i]=%12.10e\n", i, w[i]);
  }

  /*        for (i=0; i< lcp_problem->size; i++){  */
  /*     printf("zlcp[ %i]=%12.10e,\t wlcp[ %i]=%12.10e \n", i, zlcp[i],i, wlcp[i]); */
  /*        } */
  /*        printf("\n"); */

  /*        for (i=0; i< problem->size; i++){  */
  /*     printf("z[ %i]=%12.10e,\t w[ %i]=%12.10e\n", i, z[i],i, w[i]); */
  /*        } */

  options->solverId = SICONOS_RELAY_LEMKE;
  /*        printf("\n"); */
  free(zlcp);
  free(wlcp);
  freeLinearComplementarityProblem(lcp_problem);
}

/* ===========================================================================
 * Solver Registration
 * ===========================================================================
 * This registers SICONOS_RELAY_LEMKE in the global solver registry.
 */

static int relay_lexicolemke_init_wrap(void* problem, SolverOptions* options) {
  SOLVER_MAX_ITER(options) = 1000;
  SOLVER_TOL(options) = 1e-6;
  return NUMERICS_OK;
}

static int relay_lexicolemke_solve_wrap(void* problem, double* z, double* w, SolverOptions* options) {
  int info = NUMERICS_OK;
  relay_lexicolemke((RelayProblem*)problem, z, w, &info, options);
  return info;
}

static void relay_lexicolemke_free_wrap(void* problem, SolverOptions* options) {
  /* Cleanup if needed */
  (void)problem;
  (void)options;
}

REGISTER_SOLVER(SICONOS_RELAY_LEMKE, "RELAY_LEMKE",
                "Lexicographic Lemke solver for Relay problems",
                relay_lexicolemke_init_wrap,
                relay_lexicolemke_solve_wrap,
                relay_lexicolemke_free_wrap,
                NULL,  /* error function */
                1000,  /* default_max_iter */
                1e-6,  /* default_tol */
                0      /* is_local_solver */);
