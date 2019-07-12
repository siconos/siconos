/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#include "MCP_Solvers.h"
#include "MCP_cst.h"
#include "NonSmoothDrivers.h"
#include "numerics_verbose.h"

void  mcp_old_default_setDefaultSolverOptions(MixedComplementarityProblem_old* problem, SolverOptions* pOptions)
{
  pOptions->isSet = 0;
  pOptions->iSize = 10;
  pOptions->iparam = 0;
  pOptions->dSize = 10;
  pOptions->dparam = 0;
  pOptions->filterOn = 0;
  pOptions->dWork = 0;
  pOptions->iWork = 0;
  pOptions->iparam = (int*)calloc(10, sizeof(int));
  pOptions->dparam = (double*)calloc(10, sizeof(double));
  pOptions->numberOfInternalSolvers = 0;
  solver_options_nullify(pOptions);


  /*default tolerance of it*/
  pOptions->dparam[0] = 10e-7;
  /*default number of it*/
  pOptions->iparam[0] = 10;




  /* int sizeOfIwork = mcp_old_driver_get_iwork(problem, pOptions); */
  /* if(sizeOfIwork) */
  /*   pOptions->iWork = (int*)malloc(sizeOfIwork*sizeof(int)); */
  /* int sizeOfDwork = mcp_old_driver_get_dwork(problem, pOptions); */
  /* if(sizeOfDwork) */
  /*   pOptions->dWork = (double*)malloc(sizeOfDwork*sizeof(double)); */
}



int mcp_old_setDefaultSolverOptions(MixedComplementarityProblem_old* problem, SolverOptions* pOptions)
{
  int info = -1;

  switch (pOptions->solverId)
  {
  case SICONOS_MCP_FB:
  {
    info =    mcp_old_FB_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  default:
  {
    numerics_error("mixedLinearComplementarity_setDefaultSolverOptions", "Unknown Solver");
  }
  }
  return info;
}

