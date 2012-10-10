#ifndef MCP_PROBLEM_C
#define MCP_PROBLEM_C


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "MixedComplementarityProblem.h"

void freeMixedComplementarityProblem(MixedComplementarityProblem* problem)
{
  free(problem->Fmcp);
  free(problem->nablaFmcp);
  free(problem);
  problem = NULL;
}
#endif
