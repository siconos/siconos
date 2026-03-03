/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

#include <stdio.h>   // for fprintf, printf, stderr
#include <stdlib.h>  // for exit, EXIT_FAILURE

#include "AVI_Solvers.h"                    // for avi_caoferris, avi_pathavi
#include "AVI_cst.h"                        // for SICONOS_AVI_CAOFERRIS
#include "AffineVariationalInequalities.h"  // for AffineVariationalInequali...
#include "NonSmoothDrivers.h"               // for avi_driver
#include "NumericsFwd.h"                    // for SolverOptions, AffineVari...
#include "NumericsMatrix.h"                 // for NM_DENSE, NumericsMatrix
#include "SolverOptions.h"                  // for solver_options_id_to_name
#include "assert.h"                         // for assert
#include "numerics_verbose.h"
#include "utils/numerics_errors.h"

int avi_driver(AffineVariationalInequalities* problem, double* z, double* w,
               SolverOptions* options) {
  /* Input validation using standardized macros */
  CHECK_NULL(problem);
  CHECK_NULL(z);
  CHECK_NULL(w);
  CHECK_OPTIONS(options);
  CHECK_MATRIX(problem->M);
  CHECK_NULL(problem->q);

  assert(options->isSet);
  
  /* Check storage type */
  if (problem->M->storageType != NM_DENSE) {
    fprintf(stderr, "avi_driver: forbidden type of storage for the matrix M of the AVI\n");
    return NUMERICS_ERR_INVALID_ARGUMENT;
  }

  if (verbose > 0) {
    solver_options_print(options);
  }

  /* Output info. : 0: ok -  >0: problem (depends on solver) */
  int info = -1;

  if (verbose == 1)
    printf(" ========================== Call %s solver for AVI ==========================\n",
           solver_options_id_to_name(options->solverId));

  int id = options->solverId;
  switch (id) {
    case SICONOS_AVI_CAOFERRIS: {
      info = avi_caoferris(problem, z, w, options);
      break;
    }
    case SICONOS_AVI_PATHAVI: {
      info = avi_pathavi(problem, z, w, options);
      break;
    }
    /*error */
    default: {
      fprintf(stderr, "avi_driver error: unknown solver name: %s\n",
              solver_options_id_to_name(options->solverId));
      exit(EXIT_FAILURE);
    }
  }
  /*************************************************
   *  3 - Computes w = Mz + q and checks validity
   *************************************************/
  if ((options->filterOn > 0) && (info <= 0)) {
    /* info was not set or the solver was happy */
    /* TODO implement avi_compute_error, for instance evaluate the normal map*/
    /* /info = avi_compute_error(problem, z, w, options->dparam[0], &(options->dparam[1]));*/
  }
  return info;
}
