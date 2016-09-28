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
#include <float.h>
#include "VariationalInequality_Solvers.h"
#include "NonSmoothDrivers.h"
#include "numerics_verbose.h"

int variationalInequality_setDefaultSolverOptions(SolverOptions* options, int solverId)
{
  solver_options_nullify(options);

  int info = -1;
  switch (solverId)
  {
  case SICONOS_VI_EG:
  {
    info =    variationalInequality_ExtraGradient_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_VI_FPP:
  {
    info =    variationalInequality_FixedPointProjection_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_VI_HP:
  {
    info =    variationalInequality_HyperplaneProjection_setDefaultSolverOptions(options);
    break;
  }
  case SICONOS_VI_BOX_QI:
  {
    info =    variationalInequality_common_setDefaultSolverOptions(options, solverId);
    break;
  }
  default:
  {
    numerics_error("variationalInequality_setDefaultSolverOptions", "Unknown Solver");

  }
  }

  return info;
}

int variationalInequality_common_setDefaultSolverOptions(SolverOptions* options, int solverId)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for a VI Solver\n");
  }

  options->solverId = solverId;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 8;
  options->dSize = 8;
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  solver_options_nullify(options);
  for (i = 0; i < 8; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->iparam[0] = 1000;
  options->dparam[0] = 1e-10;

  options->internalSolvers = NULL;

  return 0;
}
