// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include <assert.h>  // for assert
#include <stdio.h>   // for printf
#include <stdlib.h>  // for free, malloc

#include "NumericsFwd.h"        // for RelayProblem
#include "NumericsMatrix.h"     // for NM_gemv
#include "Relay_Solvers.h"      // for relay_compute_error
#include "SiconosBlas.h"        // for cblas_daxpy, cblas_dnrm2, cblas_dcopy
#include "convex_set.h"         // for ConvexSet
#include "natural_map_error.h"  // for natural_map_error
#include "numerics_errors.h"
#include "numerics_verbose.h"  // for verbose

int relay_compute_error(RelayProblem* problem, const double* restrict z, double* restrict w,
                        double tolerance, double* restrict error) {
  /* Checks inputs */
  CHECK_NULL(problem);
  CHECK_NULL(z);
  CHECK_NULL(w);

  /* Computes w = Mz + q */
  int n = problem->size;
  cblas_dcopy(n, problem->q, 1, w, 1);  // w <-q
  NM_gemv(1.0, problem->M, z, 1.0, w);
  /* natural map z - P_{[lb, ub]}(z - w) */
  ConvexSet K = convex_set(CONVEX_SET_BOX, 1, n, problem->lb, problem->ub);
  const NaturalMapErrorScaling no_scaling = {NATURAL_MAP_ERROR_NORM_NO, 0, 0, 0.0};
  *error = natural_map_error(&no_scaling, &K, z, w);

  /* Computes error */
  double norm_q = cblas_dnrm2(n, problem->q, 1);
  *error = *error / (norm_q + 1.0);
  if (*error > tolerance) {
    if (verbose > 0)
      printf(" Numerics - relay_compute_error: error = %g > tolerance = %g.\n", *error,
             tolerance);
    return 1;
  } else
    return 0;
}
