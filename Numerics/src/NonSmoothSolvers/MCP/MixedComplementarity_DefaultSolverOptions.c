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
#include "MCP_Solvers.h"
#include "MCP_cst.h"
#include "NonSmoothDrivers.h"

void  mixedComplementarity_default_setDefaultSolverOptions(MixedComplementarityProblem* problem, SolverOptions* pOptions)
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
  pOptions->numberOfInternalSolvers = 0;

  for (int i = 0; i < 10; i++)
  {
    pOptions-> iparam[i] = 0;
    pOptions->dparam[i] = 0.0;
  }


  /*default tolerance of it*/
  pOptions->dparam[0] = 10e-7;
  /*default number of it*/
  pOptions->iparam[0] = 10;




  /* int sizeOfIwork = mcp_driver_get_iwork(problem, pOptions); */
  /* if(sizeOfIwork) */
  /*   pOptions->iWork = (int*)malloc(sizeOfIwork*sizeof(int)); */
  /* int sizeOfDwork = mcp_driver_get_dwork(problem, pOptions); */
  /* if(sizeOfDwork) */
  /*   pOptions->dWork = (double*)malloc(sizeOfDwork*sizeof(double)); */
}



int mixedComplementarity_setDefaultSolverOptions(MixedComplementarityProblem* problem, SolverOptions* pOptions)
{
  int info = -1;

  switch (pOptions->solverId)
  {
  case SICONOS_MCP_FB:
  {
    info =    mixedComplementarity_FB_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  default:
  {
    numericsError("mixedLinearComplementarity_setDefaultSolverOptions", "Unknown Solver");
  }
  }
  return info;
}

