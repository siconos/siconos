/* Siconos-Numerics, Copyright INRIA 2005-2014
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
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


