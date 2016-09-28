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
#include <time.h>
#ifndef MEXFLAG
#include "NonSmoothDrivers.h"
#endif
#include "relay_cst.h"

int dr_driver(RelayProblem* problem, double *z , double *w, SolverOptions* options)
{
  if (options == NULL)
    numerics_error("dr_driver", "null input for solver options");

  /* Checks inputs */
  if (problem == NULL || z == NULL || w == NULL)
    numerics_error("dr_driver", "null input for LinearComplementarityProblem and/or unknowns (z,w)");

  /* Output info. : 0: ok -  >0: problem (depends on solver) */
  int info = -1;

  /* Switch to DenseMatrix or SparseBlockMatrix solver according to the type of storage for M */
  /* Storage type for the matrix M of the LCP */
  int storageType = problem->M->storageType;

  /* Sparse Block Storage */
  if (storageType == 1)
  {
    numerics_error("dr_driver", "not yet implemented for sparse storage.");
  }
  // else

  /*************************************************
   *  2 - Call specific solver (if no trivial sol.)
   *************************************************/

  /* Solver name
  char * name = options->solverName;*/

  if (verbose == 1)
    printf(" ========================== Call %s solver for Relayproblem ==========================\n", solver_options_id_to_name(options->solverId));

  /****** NLGS algorithm ******/
  if (options->solverId == SICONOS_RELAY_NLGS)
    dr_nlgs(problem, z , w , &info , options);

  /****** Latin algorithm ******/
  else if (options->solverId == SICONOS_RELAY_LATIN)
    dr_latin(problem, z , w , &info , options);

  /*error */
  else
  {
    fprintf(stderr, "dr_driver error: unknown solver named: %s\n", solver_options_id_to_name(options->solverId));
    exit(EXIT_FAILURE);
  }

  return info;
}
