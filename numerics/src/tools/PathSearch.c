/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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

#include "PathSearch.h"

#include <assert.h>

#include "SolverOptions.h"

void free_solverData_PathSearch(void* solverData)
{
  assert(solverData);
  pathsearch_data* solverData_PathSearch = (pathsearch_data*) solverData;
  free_NMS_data(solverData_PathSearch->data_NMS);
  free(solverData_PathSearch->lsa_functions);
}

void pathsearch_default_SolverOption(SolverOptions* options)
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
}


