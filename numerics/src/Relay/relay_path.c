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
#include "NumericsFwd.h"    // for SolverOptions, RelayProblem
#include "Relay_Solvers.h"  // for relay_path, relay_path_setDefaultSolverOp...
#include "Relay_options.h"
#include "naming_conventions.h"
#include "numerics_errors.h"
#include "solver_registry.h"

#ifdef HAVE_PATHFERRIS
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "NumericsMatrix.h"
#include "Relay_options.h"
#include "SimpleLCP.h"
#include "numerics_verbose.h"
#endif /*HAVE_PATHFERRIS*/

void relay_path(RelayProblem *problem, double *z, double *w, int *info,
                SolverOptions *options) {
  *info = 1;
#ifdef HAVE_PATHFERRIS
  /* matrix M/vector q of the relay */
  if (!problem->M->matrix0) {
    printf("relay_path only implemented for dense storage\n");
    //      return info;
  }

  double *M = problem->M->matrix0;
  double *q = problem->q;

  int nnz, i, j;
  /* size of the RELAY */
  int n = problem->size;

  //  double tol = options->dparam[SICONOS_DPARAM_TOL];
  MCP_Termination termination;

  nnz = nbNonNulElems(n, M, 1.0e-18);
  int *m_i = (int *)calloc(nnz + 1, sizeof(int));
  int *m_j = (int *)calloc(nnz + 1, sizeof(int));
  double *m_ij = (double *)calloc(nnz + 1, sizeof(double));
  double *lb = (double *)calloc(n + 1, sizeof(double));
  double *ub = (double *)calloc(n + 1, sizeof(double));
  double err = 1e24, val;

  FortranToPathSparse(n, M, 1.0e-18, m_i, m_j, m_ij);
  for (i = 0; i < n; i++) {
    lb[i] = problem->lb[i];
    ub[i] = problem->ub[i];
  }
  /*   for (i=0;i<n;i++){ */
  /*       printf(" problem->lb[%i] = %e",i, problem->lb[i]); */
  /*       printf(" problem->ub[%i] = %e",i, problem->ub[i]); */
  /*   } */
  SimpleLCP(n, nnz, m_i, m_j, m_ij, q, lb, ub, &termination, z);
  /*   printLCP(n, nnz, m_i, m_j, m_ij, q, lb, ub); */

  if (termination == MCP_Error) {
    *info = 1;
    if (verbose > 0) printf("PATH : Error in the solution.\n");
  } else if (termination == MCP_Solved) {
    for (i = 0; i < n; i++) {
      val = q[i];
      for (j = 0; j < n; j++) {
        val += M[i + j * n] * z[j];
      }
      w[i] = val;
    }
    *info = 0;
    /* **** Criterium convergence **** */
    // relay_compute_error(problem,z,w,tol,&err);

    if (verbose > 0) printf("PATH : RELAY Solved, error %10.7f.\n", err);
  } else {
    if (verbose > 0) printf("PATH : Other error: %d\n", termination);
  }
  free(m_i);
  free(m_j);
  free(m_ij);
  free(lb);
  free(ub);

#endif /*HAVE_PATHFERRIS*/

  return;
}

/* ===========================================================================
 * Solver Registration
 * ===========================================================================
 * This registers SICONOS_RELAY_PATH in the global solver registry.
 */

static void relay_path_set_default(SolverOptions *options) {
  SOLVER_MAX_ITER(options) = 1000;
  SOLVER_TOL(options) = 1e-6;
}

static int relay_path_init_wrap(void *problem, SolverOptions *options) {
  (void)problem;
  return NUMERICS_OK;
}

static int relay_path_solve_wrap(void *problem, double *z, double *w, SolverOptions *options) {
  int info = NUMERICS_OK;
  relay_path((RelayProblem *)problem, z, w, &info, options);
  return info;
}

static void relay_path_free_wrap(void *problem, SolverOptions *options) {
  /* Cleanup if needed */
  (void)problem;
  (void)options;
}

REGISTER_SOLVER(SICONOS_RELAY_PATH, "RELAY_PATH",
                "PATH solver for Relay problems (MCP formulation)", relay_path_init_wrap,
                relay_path_solve_wrap, relay_path_free_wrap, NULL, /* error function */
                relay_path_set_default, 1000,                      /* default_max_iter */
                1e-6,                                              /* default_tol */
                0 /* is_local_solver */);
