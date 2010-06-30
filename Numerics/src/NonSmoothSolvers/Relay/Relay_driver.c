/* Siconos-Numerics, Copyright INRIA 2005-2010.
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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef MEXFLAG
#include "NonSmoothDrivers.h"
#endif
#include <time.h>

char  SICONOS_RELAY_PGS_STR[] = "RELAY_PGS";
char  SICONOS_RELAY_PATH_STR[] = "RELAY_PATH";
char  SICONOS_RELAY_ENUM_STR[] = "RELAY_ENUM";
char  SICONOS_RELAY_NLGS_STR[] = "RELAY_NLGS";
char  SICONOS_RELAY_LEMKE_STR[] = "RELAY_LEMKE";
char  SICONOS_RELAY_LATIN_STR[] = "RELAY_LATIN";

int relay_driver(RelayProblem* problem, double *z , double *w,
                 SolverOptions* options, NumericsOptions* global_options)
{


  //Relay_display(problem);

  if (options == NULL || global_options == NULL)
    numericsError("Relay_driver", "null input for solver and/or global options");

  /* Set global options */
  setNumericsOptions(global_options);

  /* Checks inputs */
  if (problem == NULL || z == NULL || w == NULL)
    numericsError("Relay_driver", "null input for RelayProblem and/or unknowns (z,w)");

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
  if (verbose > 0)
    printSolverOptions(options);

  /* Solver name */
  //char * name = options->solverName;

  if (verbose == 1)
    printf(" ========================== Call %s solver for Relayproblem ==========================\n", idToName(options->solverId));

  switch (options->solverId)
  {
  case SICONOS_RELAY_PGS:
  {
    relay_pgs(problem, z , w , &info , options);
    break;
  }
  case SICONOS_RELAY_NLGS:
  {
    fprintf(stderr, "Relay_driver error: NLGS solver obsolete use PGS:\n");
    break;
  }
  case SICONOS_RELAY_LEMKE:
  {

    char filename[20] = "relay_simple.dat";

    FILE *FP = fopen(filename, "w");
    relay_printInFile(problem, FP);
    fclose(FP);
    relay_lexicolemke(problem, z , w , &info , options, global_options);
    break;
  }
  case SICONOS_RELAY_ENUM:
  {
    relay_enum(problem, z , w , &info , options, global_options);
    break;
  }
  case SICONOS_RELAY_PATH:
  {
    relay_path(problem, z , w , &info , options);
    break;
  }
  /*error */
  default:
  {
    fprintf(stderr, "Relay_driver error: unknown solver name: %s\n", idToName(options->solverId));
    exit(EXIT_FAILURE);
  }
  }
  if (options[0].filterOn > 0)
    info = relay_compute_error(problem, z, w, options[0].dparam[0], &(options[0].dparam[1]));

  return info;
}


