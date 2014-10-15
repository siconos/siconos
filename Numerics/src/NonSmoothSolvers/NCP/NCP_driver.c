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

#include "SolverOptions.h"
#include <assert.h>

#include "NCP_Solvers.h"
#include "NCP_cst.h"

char SICONOS_NCP_NEWTON_FBLSA_STR[] = "NCP Newton FBLSA";
char SICONOS_NCP_NEWTON_MINFBLSA_STR[] = "NCP Newton minFBLSA";
char SICONOS_NCP_PATHSEARCH_STR[] = "NCP Path search";

int ncp_driver(NCP_struct* problem, double *z , double *F, SolverOptions* options,  NumericsOptions* global_options)
{
  assert(options && "ncp_driver null input for solver options.\n");

  /* Set global options */
  if (global_options)
    setNumericsOptions(global_options);

  /* Checks inputs */
  assert(problem && z && F && "ncp_driver null input for MixedComplementarityProblem and/or unknowns (z,w)");

  /* Output info. : 0: ok -  >0: error (which depends on the chosen solver) */
  int info = -1;

  switch (options->solverId)
  {
  case SICONOS_NCP_NEWTON_FBLSA: // Fischer-Burmeister + Newton w/ LS
    ncp_newton_FBLSA(problem, z, F, &info, options);
    break;
  case SICONOS_NCP_NEWTON_MINFBLSA: // min (+ FB as backup) + Newton w/ LS
    ncp_newton_minFBLSA(problem, z, F, &info, options);
    break;
  case SICONOS_NCP_PATHSEARCH: // pathsearch method
    ncp_pathsearch(problem, z, F, &info, options);
    break;
  default:
    fprintf(stderr, "ncp_driver error: unknown solver id: %d\n", options->solverId);
    exit(EXIT_FAILURE);
  }

  /* check the conditions 0 <= z _|_ F(z) >= 0 */
  if (options->filterOn > 0)
  {
    int info_ = ncp_compute_error(problem->n, z, F, options->dparam[0], &(options->dparam[1]));
    if (info <= 0) /* info was not set or the solver was happy */
      info = info_;
  }

  return info;
}
