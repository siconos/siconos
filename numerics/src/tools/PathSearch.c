/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include "PathSearch.h"

#include <assert.h>

#include "SolverOptions.h"  // for SolverOptions
#include "lcp_cst.h"
#include "line_search.h"  // for ARCSEARCH, LINESEARCH, NM_LS_MEAN
#include "numerics_errors.h"
void pathsearch_set_default(SolverOptions* options) {
  options->iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS] = NM_LS_MEAN;
  options->iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS_M] = 10;
  options->iparam[SICONOS_IPARAM_PATHSEARCH_STACKSIZE] = 5;
  options->iparam[SICONOS_IPARAM_NMS_WATCHDOG_TYPE] = LINESEARCH;
  options->iparam[SICONOS_IPARAM_NMS_PROJECTED_GRADIENT_TYPE] = ARCSEARCH;
  options->iparam[SICONOS_IPARAM_NMS_N_MAX] = 10;

  options->dparam[SICONOS_DPARAM_NMS_DELTA] = 20;
  options->dparam[SICONOS_DPARAM_NMS_DELTA_VAR] = .8;
  options->dparam[SICONOS_DPARAM_NMS_SIGMA] = .01;
  options->dparam[SICONOS_DPARAM_NMS_ALPHA_MIN_WATCHDOG] = 1e-12;
  options->dparam[SICONOS_DPARAM_NMS_ALPHA_MIN_PGRAD] = 1e-12;
  options->dparam[SICONOS_DPARAM_NMS_MERIT_INCR] = 1.1; /* XXX 1.1 ?*/

  assert(options->numberOfInternalSolvers == 1);
  options->internalSolvers[0] = solver_options_create(SICONOS_LCP_PIVOT);

  SolverOptions* lcp_options = options->internalSolvers[0];

  /* We always allocate once and for all since we are supposed to solve
   * many LCPs */
  lcp_options->iparam[SICONOS_IPARAM_PREALLOC] = 1;
  /* set the right pivot rule */
  lcp_options->iparam[SICONOS_LCP_IPARAM_PIVOTING_METHOD_TYPE] = SICONOS_LCP_PIVOT_PATHSEARCH;
  /* set the right stacksize */
  lcp_options->iparam[SICONOS_IPARAM_PATHSEARCH_STACKSIZE] =
      options->iparam[SICONOS_IPARAM_PATHSEARCH_STACKSIZE];
}
