/* Siconos-Numerics, Copyright INRIA 2005-2012.
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
#include <time.h>
#include <float.h>

#include "NumericsOptions.h"
#include "VariationalInequality_Solvers.h"
#include "NonSmoothDrivers.h"

int variationalInequality_setDefaultSolverOptions(SolverOptions* options, int solverId)
{
  options->iparam = NULL;
  options->dparam = NULL;
  options->dWork = NULL;
  options->iWork = NULL;
  options->internalSolvers = NULL;
  options->numericsOptions = NULL;
  options->callback = NULL;

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
    numericsError("variationalInequality_setDefaultSolverOptions", "Unknown Solver");

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
  options->dWork = NULL;
  null_SolverOptions(options);
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
