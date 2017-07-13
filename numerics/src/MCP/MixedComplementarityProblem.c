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

#ifndef MCP_PROBLEM_C
#define MCP_PROBLEM_C


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "MixedComplementarityProblem.h"
#include "NumericsMatrix.h"

void freeMixedComplementarityProblem(MixedComplementarityProblem* problem)
{
//  if (problem->Fmcp) free(problem->Fmcp);
//  if (problem->nablaFmcp) free(problem->nablaFmcp);
  free(problem);
}

void freeMCP(MixedComplementarityProblem2* mcp)
{
  if (mcp->nabla_Fmcp)
  {
    NM_free(mcp->nabla_Fmcp);
    free(mcp->nabla_Fmcp);
    mcp->nabla_Fmcp = NULL;
  }

  free(mcp);
}

MixedComplementarityProblem2* newMCP(void)
{
  MixedComplementarityProblem2* mcp = (MixedComplementarityProblem2*) malloc(sizeof(MixedComplementarityProblem2));

  mcp->n1 = 0;
  mcp->n2 = 0;
  mcp->compute_Fmcp = NULL;
  mcp->compute_nabla_Fmcp = NULL;
  mcp->nabla_Fmcp = NULL;
  mcp->env = NULL;

  return mcp;
}


#endif
