/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

#include <assert.h>
#include "PathSearch.h"
#include "SolverOptions.h"  // for SolverOptions
#include "line_search.h"    // for ARCSEARCH, LINESEARCH, NM_LS_MEAN
#include "lcp_cst.h"

void pathsearch_set_default(SolverOptions* options)
{
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

  SolverOptions * lcp_options = options->internalSolvers[0];

  /* We always allocate once and for all since we are supposed to solve
   * many LCPs */
  lcp_options->iparam[SICONOS_IPARAM_PREALLOC] = 1;
  /* set the right pivot rule */
  lcp_options->iparam[SICONOS_LCP_IPARAM_PIVOTING_METHOD_TYPE] = SICONOS_LCP_PIVOT_PATHSEARCH;
  /* set the right stacksize */
  lcp_options->iparam[SICONOS_IPARAM_PATHSEARCH_STACKSIZE] = options->iparam[SICONOS_IPARAM_PATHSEARCH_STACKSIZE];

  
}


