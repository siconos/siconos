/* Siconos-Numerics version 3.0.0, Copyright INRIA 2005-2008.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef MEXFLAG
#include "NonSmoothDrivers.h"
#endif
#include <time.h>

int relay_driver(Relay_Problem* problem, double *z , double *w, Solver_Options* options, Numerics_Options* global_options)
{
  if (options == NULL || global_options == NULL)
    numericsError("Relay_driver", "null input for solver and/or global options");

  /* Set global options */
  setNumericsOptions(global_options);

  /* Checks inputs */
  if (problem == NULL || z == NULL || w == NULL)
    numericsError("Relay_driver", "null input for Relay_Problem and/or unknowns (z,w)");

  /* Output info. : 0: ok -  >0: problem (depends on solver) */
  int info = -1;

  /* Switch to DenseMatrix or SparseBlockMatrix solver according to the type of storage for M */
  /* Storage type for the matrix M of the LCP */
  int storageType = problem->M->storageType;

  /* Sparse Block Storage */
  if (storageType == 1)
  {
    numericsError("Relay_driver", "not yet implemented for sparse storage.");
  }
  // else

  /*************************************************
   *  2 - Call specific solver (if no trivial sol.)
   *************************************************/

  /* Solver name */
  char * name = options->solverName;

  if (verbose == 1)
    printf(" ========================== Call %s solver for Relayproblem ==========================\n", name);

  /****** NLGS algorithm ******/
  if (strcmp(name , "PGS") == 0)
    relay_pgs(problem, z , w , &info , options);

  else if (strcmp(name , "NLGS") == 0)
    fprintf(stderr, "Relay_driver error: NLGS solver obsolete:\n");

  /****** Latin algorithm ******/
  else if (strcmp(name , "Latin") == 0)
    pr_latin(problem, z , w , &info , options);

  else if (strcmp(name , "path") == 0)
    relay_path(problem, z , w , &info , options);

  /*error */
  else
  {
    fprintf(stderr, "Relay_driver error: unknown solver name: %s\n", name);
    exit(EXIT_FAILURE);
  }
  if (options[0].filterOn > 0)
    info = relay_compute_error(problem, z, w, options[0].dparam[0], &(options[0].dparam[1]));

  return info;
}

