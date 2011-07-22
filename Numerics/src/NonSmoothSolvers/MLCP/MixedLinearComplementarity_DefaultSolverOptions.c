/* Siconos-Numerics, Copyright INRIA 2005-2011.
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
#include "LA.h"
#include "NumericsOptions.h"
#include "MLCP_Solvers.h"
#include "NonSmoothDrivers.h"



void  mixedLinearComplementarity_default_setDefaultSolverOptions(MixedLinearComplementarityProblem* problem, SolverOptions* pOptions)
{
  pOptions->isSet = 0;
  pOptions->iSize = 10;
  pOptions->iparam = 0;
  pOptions->dSize = 10;
  pOptions->dparam = 0;
  pOptions->filterOn = 0;
  pOptions->dWork = 0;
  pOptions->iWork = 0;
  pOptions->iparam = (int*)malloc(10 * sizeof(int));
  pOptions->dparam = (double*)malloc(10 * sizeof(double));
  pOptions->numberOfInternalSolvers = 1;


  pOptions->dparam[0] = 10 - 7;
  /*default number of it*/
  pOptions->iparam[0] = 1000;
  /*enum case : do not use dgels*/
  pOptions->iparam[4] = 0;
  pOptions->iparam[5] = 3; /*Number of registered configurations*/
  pOptions->iparam[8] = 0; /*Prb nedd a update*/
  pOptions->dparam[5] = 1e-12; /*tol used by direct solver to check complementarity*/
  pOptions->dparam[6] = 1e-12; /*tol for direct solver to determinate if a value is positive*/

  int sizeOfIwork = mlcp_driver_get_iwork(problem, pOptions);
  if (sizeOfIwork)
    pOptions->iWork = (int*)malloc(sizeOfIwork * sizeof(int));
  int sizeOfDwork = mlcp_driver_get_dwork(problem, pOptions);
  if (sizeOfDwork)
    pOptions->dWork = (double*)malloc(sizeOfDwork * sizeof(double));



}

void  mixedLinearComplementarity_deleteDefaultSolverOptions(MixedLinearComplementarityProblem* problem, SolverOptions* pOptions)
{
  if (pOptions->iparam)
    free(pOptions->iparam);
  if (pOptions->dparam)
    free(pOptions->dparam);
  if (pOptions->iWork)
    free(pOptions->iWork);
  if (pOptions->dWork)
    free(pOptions->dWork);
  pOptions->iparam = NULL;
  pOptions->dparam = NULL;
  pOptions->iWork = NULL;
  pOptions->dWork = NULL;

}

int mixedLinearComplementarity_setDefaultSolverOptions(MixedLinearComplementarityProblem* problem, SolverOptions* pOptions)
{
  int info = -1;

  switch (pOptions->solverId)
  {
  case SICONOS_MLCP_DIRECT_ENUM:
  {
    info =    mixedLinearComplementarity_directEnum_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  case SICONOS_MLCP_PATH_ENUM:
  {
    info =    mixedLinearComplementarity_pathEnum_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  case  SICONOS_MLCP_DIRECT_PATH_ENUM:
  {
    info =    mixedLinearComplementarity_directPathEnum_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  case SICONOS_MLCP_DIRECT_SIMPLEX:
  {
    info =    mixedLinearComplementarity_directSimplex_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  case SICONOS_MLCP_DIRECT_PATH:
  {
    info =    mixedLinearComplementarity_directPath_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  case SICONOS_MLCP_DIRECT_FB:
  {
    info =    mixedLinearComplementarity_directFB_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  case SICONOS_MLCP_SIMPLEX:
  {
    info =    mixedLinearComplementarity_simplex_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  case SICONOS_MLCP_PGS:
  {
    info =    mixedLinearComplementarity_pgs_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  case SICONOS_MLCP_RPGS:
  {
    info =    mixedLinearComplementarity_rpgs_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  case SICONOS_MLCP_RPSOR:
  {
    info =    mixedLinearComplementarity_rpsor_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  case SICONOS_MLCP_PATH:
  {
    info =    mixedLinearComplementarity_path_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  case SICONOS_MLCP_ENUM:
  {
    info =    mixedLinearComplementarity_enum_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  case SICONOS_MLCP_FB:
  {
    info =    mixedLinearComplementarity_fb_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  case SICONOS_MLCP_PSOR:
  {
    info = mixedLinearComplementarity_psor_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  default:
  {
    numericsError("mixedLinearComplementarity_setDefaultSolverOptions", "Unknown Solver");

  }
  }
  return info;
}

