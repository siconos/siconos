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
#include "MCP_Solvers.h"       // for mcp_newton_FB_FBLSA_setDefaultSolverOp...
#include "MCP_cst.h"           // for SICONOS_MCP_NEWTON_FB_FBLSA, SICONOS_M...
#include "NumericsFwd.h"       // for SolverOptions, MixedComplementarityPro...
#include "SolverOptions.h"     // for SolverOptions
#include "numerics_verbose.h"  // for numerics_error

int mcp_setDefaultSolverOptions(
  MixedComplementarityProblem* problem,
  SolverOptions* options)
{
  int info = -1;

  switch (options->solverId)
  {
  case SICONOS_MCP_NEWTON_FB_FBLSA:
  {
    info =     mcp_newton_FB_FBLSA_setDefaultSolverOptions(problem, options);
    break;
  }
  case SICONOS_MCP_NEWTON_MIN_FBLSA:
  {
    info =     mcp_newton_min_FBLSA_setDefaultSolverOptions(problem, options);
    break;
  }
  default:
  {
    numerics_error("mcp_setDefaultSolverOptions", "Unknown Solver");
  }
  }
  return info;
}


int mcp_old_setDefaultSolverOptions(MixedComplementarityProblem_old* problem, SolverOptions* pOptions)
{
  int info = -1;

  switch (pOptions->solverId)
  {
  case SICONOS_MCP_OLD_FB:
  {
    info =    mcp_old_FB_setDefaultSolverOptions(problem, pOptions);
    break;
  }
  default:
  {
    numerics_error("mcp_old_setDefaultSolverOptions", "Unknown Solver");
  }
  }
  return info;
}

