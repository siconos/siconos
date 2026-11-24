// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef MCP_PROBLEM_C
#define MCP_PROBLEM_C
#include "MixedComplementarityProblem.h"

#include <stdio.h>   // for NULL
#include <stdlib.h>  // for free, malloc

#include "NumericsMatrix.h"  // for NM_clear
#include "numerics_errors.h"

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
