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

#ifndef MCP_PROBLEM_C
#define MCP_PROBLEM_C
#include "MixedComplementarityProblem.h"

#include <stdio.h>   // for NULL
#include <stdlib.h>  // for free, malloc

#include "NumericsMatrix.h"  // for NM_clear

void mixedComplementarityProblem_old_free(MixedComplementarityProblem_old* problem) {
  //  if (problem->Fmcp) free(problem->Fmcp);
  //  if (problem->nablaFmcp) free(problem->nablaFmcp);
  free(problem);
}

void mixedComplementarityProblem_free(MixedComplementarityProblem* mcp) {
  if (mcp->nabla_Fmcp) {
    NM_clear(mcp->nabla_Fmcp);
    free(mcp->nabla_Fmcp);
    mcp->nabla_Fmcp = NULL;
  }

  free(mcp);
}

MixedComplementarityProblem* mixedComplementarityProblem_new(void) {
  MixedComplementarityProblem* mcp =
      (MixedComplementarityProblem*)malloc(sizeof(MixedComplementarityProblem));

  mcp->n1 = 0;
  mcp->n2 = 0;
  mcp->compute_Fmcp = NULL;
  mcp->compute_nabla_Fmcp = NULL;
  mcp->nabla_Fmcp = NULL;
  mcp->env = NULL;

  return mcp;
}

#endif
