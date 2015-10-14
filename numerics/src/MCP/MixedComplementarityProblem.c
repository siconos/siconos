#ifndef MCP_PROBLEM_C
#define MCP_PROBLEM_C


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "MixedComplementarityProblem.h"

void freeMixedComplementarityProblem(MixedComplementarityProblem* problem)
{
//  if (problem->Fmcp) free(problem->Fmcp);
//  if (problem->nablaFmcp) free(problem->nablaFmcp);
  free(problem);
}
#endif
