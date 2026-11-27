// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include <float.h>   // for DBL_EPSILON
#include <math.h>    // for fabs, fmax, pow, sqrt
#include <stddef.h>  // for NULL

#include "LCP_Solvers.h"                   // for lcp_compute_error, lcp_com...
#include "LinearComplementarityProblem.h"  // for LinearComplementarityProblem
#include "NumericsFwd.h"                   // for LinearComplementarityProblem
#include "NumericsMatrix.h"                // for NM_gemv
#include "SiconosBlas.h"                   // for cblas_dcopy, cblas_dnrm2
#include "convex_set.h"                    // for ConvexSet
#include "natural_map_error.h"             // for natural_map_error
#include "numerics_errors.h"
#include "numerics_verbose.h"  // for numerics_error, numerics_p...

void lcp_compute_error_only(unsigned int n, const double* restrict z, const double* restrict w,
                            double* restrict error) {
  /* natural map z - P_{R^n_+}(z - w) = z - max(0, z - w) */
  ConvexSet K = convex_set(CONVEX_SET_ORTHANT, 1, (int)n, NULL, NULL);
  const NaturalMapErrorScaling no_scaling = {NATURAL_MAP_ERROR_NORM_NO, 0, 0, 0.0};
  *error = natural_map_error(&no_scaling, &K, z, w);
}

int lcp_compute_error(LinearComplementarityProblem* problem, double* z, double* w,
                      double tolerance, double* error) {
  /* Input validation using standardized macros */
  CHECK_NULL(problem);
  CHECK_NULL(z);
  CHECK_NULL(w);
  CHECK_NULL(error);
  CHECK_MATRIX(problem->M);
  CHECK_NULL(problem->q);

  /* Computes w = Mz + q */
  int incx = 1, incy = 1;
  unsigned int n = problem->size;
  cblas_dcopy(n, problem->q, incx, w, incy);  // w <-q
  NM_gemv(1.0, problem->M, z, 1.0, w);
  double norm_q = cblas_dnrm2(n, problem->q, incx);
  lcp_compute_error_only(n, z, w, error);

  if (fabs(norm_q) > DBL_EPSILON) *error /= norm_q;

  if (*error > tolerance) {
    numerics_printf("lcp_compute_error : error = %g > tolerance = %g.", *error, tolerance);
    return 1;
  } else
    return 0;
}
