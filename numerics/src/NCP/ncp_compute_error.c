/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#include "LCP_Solvers.h"  // for lcp_compute_error_only
#include "NCP_Solvers.h"  // for ncp_compute_error
#include "numerics_errors.h"
#include "siconos_debug.h"  // for DEBUG_END, DEBUG_BEGIN, DEBUG_PRINTF

int ncp_compute_error(int n, double* z, double* F, double tol, double* err) {
  DEBUG_BEGIN("ncp_compute_error(int n, double* z, double * F, double tol, double* err)\n")
  lcp_compute_error_only(n, z, F, err);

  DEBUG_PRINTF("ncp_compute_error err= %g, tol =%g \n", *err, tol);

  if (*err >= tol) {
    DEBUG_END("ncp_compute_error(int n, double* z, double * F, double tol, double* err)\n");
    return 1;
  } else {
    DEBUG_END("ncp_compute_error(int n, double* z, double * F, double tol, double* err)\n");
    return 0;
  }
}
