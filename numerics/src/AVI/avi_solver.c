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

#include "NonSmoothDrivers.h"
#include "AVI_Solvers.h"
#include "AVI_cst.h"
#include "numerics_verbose.h"
#include "assert.h"
#include "NumericsMatrix.h"

const char* const SICONOS_AVI_CAOFERRIS_STR = "AVI from Cao & Ferris";
const char* const SICONOS_AVI_PATHAVI_STR = "PATHVI";

int avi_driver(AffineVariationalInequalities* problem, double *z , double *w, SolverOptions* options)
{

  assert(options && "avi_driver : null input for solver options");
  /* Checks inputs */
  assert(problem && z && w &&
      "avi_driver : input for LinearComplementarityProblem and/or unknowns (z,w)");

  assert(problem->M->storageType == NM_DENSE &&
      "avi_driver_DenseMatrix : forbidden type of storage for the matrix M of the AVI");

  assert(options->isSet);

  if (!problem || !problem->M || !problem->q || !options)
  {
    numerics_error_nonfatal("avi_driver", "Problem data is incomplete: you need to set at least M, q and poly");
    return -1;
  }

  if (verbose > 0)
  {
    solver_options_print(options);
  }

  /* Output info. : 0: ok -  >0: problem (depends on solver) */
  int info = -1;

  if (verbose == 1)
    printf(" ========================== Call %s solver for AVI ==========================\n", solver_options_id_to_name(options->solverId));

  int id = options->solverId;
  switch (id)
  {
  case SICONOS_AVI_CAOFERRIS:
  {
    info = avi_caoferris(problem, z, w, options);
    break;
  }
  case SICONOS_AVI_PATHAVI:
  {
    info = avi_pathavi(problem, z, w, options);
    break;
  }
  /*error */
  default:
  {
    fprintf(stderr, "avi_driver error: unknown solver name: %s\n", solver_options_id_to_name(options->solverId));
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
