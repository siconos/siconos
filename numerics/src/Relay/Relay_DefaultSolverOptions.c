/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.

 * Copyright 2016 INRIA.

 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at

 * http://www.apache.org/licenses/LICENSE-2.0

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
#include <float.h>
#include "NumericsOptions.h"
#include "Relay_Solvers.h"
#include "NonSmoothDrivers.h"
#include "relay_cst.h"

int relay_setDefaultSolverOptions(RelayProblem* problem, SolverOptions* options, int solverId)
{

  int info = -1;
  switch (solverId)
  {
  case SICONOS_RELAY_PGS:
  {
    info =    relay_pgs_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_RELAY_LEMKE:
  {
    info =    relay_lexicolemke_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_RELAY_ENUM:
  {
    info =    relay_enum_setDefaultSolverOptions(problem, options);
    break;
  }
  case SICONOS_RELAY_PATH:
  {
    info =    relay_path_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_RELAY_AVI_CAOFERRIS:
  {
    info =    relay_avi_caoferris_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_RELAY_AVI_CAOFERRIS_TEST:
  {
    info =    relay_avi_caoferris_test_setDefaultSolverOptions(options);
    break;
  }
/* XXX: to implement ?
   case SICONOS_RELAY_LATIN:
  {
    info =    relay_latin_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_RELAY_NLGS:
  {
    info =    relay_nlgs_setDefaultSolverOptions(options);
    break;
  }
  */
  default:
  {
    numericsError("Relay_setDefaultSolverOptions", " Unknown Solver");
  }
  }

  return info;
}
