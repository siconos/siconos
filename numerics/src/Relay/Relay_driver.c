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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef MEXFLAG
#include "NonSmoothDrivers.h"
#endif
#include "relay_cst.h"
#include <time.h>
#include "misc.h"

char*  SICONOS_RELAY_PGS_STR = "RELAY_PGS";
char*  SICONOS_RELAY_PATH_STR = "RELAY_PATH";
char*  SICONOS_RELAY_ENUM_STR = "RELAY_ENUM";
char*  SICONOS_RELAY_NLGS_STR = "RELAY_NLGS";
char*  SICONOS_RELAY_LEMKE_STR = "RELAY_LEMKE";
char*  SICONOS_RELAY_LATIN_STR = "RELAY_LATIN";
char*  SICONOS_RELAY_AVI_CAOFERRIS_STR = "RELAY_AVI_CAOFERRIS";
char*  SICONOS_RELAY_AVI_CAOFERRIS_TEST_STR = "test version of the solver by Cao & Ferris; DO NOT USE!";

int relay_driver(RelayProblem* problem, double *z , double *w,
                 SolverOptions* options)
{


  //Relay_display(problem);

  if (options == NULL)
    numericsError("Relay_driver", "null input for solver and/or global options");

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
    solver_options_print(options);

  /* Solver name */
  //char * name = options->solverName;

  if (verbose == 1)
    printf(" ========================== Call %s solver for Relayproblem ==========================\n", solver_options_id_to_name(options->solverId));

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

#ifdef DEBUG_RELAY
    char filename[20] = "relay_simple.dat";

    FILE *FP = fopen(filename, "w");
    relay_printInFile(problem, FP);
    fclose(FP);
#endif
    relay_lexicolemke(problem, z , w , &info , options);
    break;
  }
  case SICONOS_RELAY_ENUM:
  {
    relay_enum(problem, z , w , &info , options);
    break;
  }
  case SICONOS_RELAY_PATH:
  {
    relay_path(problem, z , w , &info , options);
    break;
  }
  case SICONOS_RELAY_AVI_CAOFERRIS:
  {
    relay_avi_caoferris(problem, z , w , &info , options);
    break;
  }
  case SICONOS_RELAY_AVI_CAOFERRIS_TEST:
  {
    relay_avi_caoferris_test(problem, z , w , &info , options);
    break;
  }
  /*error */
  default:
  {
    fprintf(stderr, "Relay_driver error: unknown solver name: %s\n", solver_options_id_to_name(options->solverId));
    exit(EXIT_FAILURE);
  }
  }
  if (options[0].filterOn > 0)
    info = relay_compute_error(problem, z, w, options[0].dparam[0], &(options[0].dparam[1]));

  return info;
}


