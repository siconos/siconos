/* Siconos-Numerics, Copyright INRIA 2005-2015
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

#include "NumericsOptions.h"
#include "NonSmoothDrivers.h"

#include "AVI_Solvers.h"

char *  SICONOS_AVI_CAOFERRIS_STR = "AVI from Cao & Ferris";

int avi_driver(AffineVariationalInequalities* problem, double *z , double *w, SolverOptions* options,  NumericsOptions* global_options)
{
  assert(options && "avi_driver : null input for solver options");

  /* Set global options */
  if (global_options)
  {
    setNumericsOptions(global_options);
  }

  /* Checks inputs */
  assert(problem && z && w &&
      "avi_driver : input for LinearComplementarityProblem and/or unknowns (z,w)");

  assert(problem->M->storageType == 0 &&
      "avi_driver_DenseMatrix : forbidden type of storage for the matrix M of the AVI");

  /* If the options for solver have not been set, read default values in .opt file */
  if (options->isSet == 0)
  {
    readSolverOptions(0, options);
    options->filterOn = 1;
  }

  if (verbose > 0)
  {
    printSolverOptions(options);
  }

  /* Output info. : 0: ok -  >0: problem (depends on solver) */
  int info = -1;

  if (verbose == 1)
    printf(" ========================== Call %s solver for AVI ==========================\n", idToName(options->solverId));

  int id = options->solverId;
  switch (id)
  {
  case SICONOS_AVI_CAOFERRIS:
    info = avi_caoferris(problem, z, w, options);
    break;
  /*error */
  default:
  {
    fprintf(stderr, "avi_driver error: unknown solver name: %s\n", idToName(options->solverId));
    exit(EXIT_FAILURE);
  }
  }
  /*************************************************
   *  3 - Computes w = Mz + q and checks validity
   *************************************************/
  if ((options->filterOn > 0) && (info <= 0))
  {
    /* info was not set or the solver was happy */
    /* TODO implement avi_compute_error, for instance evaluate the normal map*/
    /* /info = avi_compute_error(problem, z, w, options->dparam[0], &(options->dparam[1]));*/
  }
  return info;

}
