/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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
#include <assert.h>            // for assert
#include <stdio.h>             // for fprintf, NULL, stderr
#include <stdlib.h>            // for exit, EXIT_FAILURE
#include "MCP_Solvers.h"       // for mcp_newton_FB_FBLSA, mcp_newton_min_FBLSA
#include "MCP_cst.h"           // for SICONOS_MCP_OLD_FB, SICONOS_MCP_NEWTON...
#include "NonSmoothDrivers.h"  // for mcp_driver, mcp_old_driver
#include "NumericsFwd.h"       // for SolverOptions, MixedComplementarityPro...
#include "SolverOptions.h"     // for SolverOptions
#include "numerics_verbose.h"  // for numerics_error

const char* const SICONOS_MCP_OLD_FB_STR = "NewtonFB";
const char* const SICONOS_MCP_NEWTON_FB_FBLSA_STR = "MCP Newton FBLSA";
const char* const SICONOS_MCP_NEWTON_MIN_FBLSA_STR = "MCP Newton minFBLSA";

int mcp_driver(MixedComplementarityProblem* problem, double *z, double *Fmcp, SolverOptions* options)
{
  assert(options != NULL);
  /* Checks inputs */
  assert(problem != NULL);
  assert(z != NULL);
  assert(Fmcp != NULL);
  /* Output info. : 0: ok -  >0: error (which depends on the chosen solver) */
  int info = -1;

  switch(options->solverId)
  {
  case SICONOS_MCP_NEWTON_FB_FBLSA: // Fischer-Burmeister/Newton -- new version
    mcp_newton_FB_FBLSA(problem, z, Fmcp, &info, options);
    break;

  case SICONOS_MCP_NEWTON_MIN_FBLSA: // Fischer-Burmeister/Newton + min descent direction
    mcp_newton_min_FBLSA(problem, z, Fmcp, &info, options);
    break;

  default:
    fprintf(stderr, "mcp_driver error: unknown solver id: %d\n", options->solverId);
    exit(EXIT_FAILURE);
  }

  return info;
}

int mcp_old_driver(MixedComplementarityProblem_old* problem, double *z, double *w, SolverOptions* options)
{
  if(options == NULL)
    numerics_error("mcp_old_driver ", "null input for solver options.\n");

  /* Checks inputs */
  if(problem == NULL || z == NULL || w == NULL)
    numerics_error("mcp_old_driver", "null input for MixedComplementarityProblem_old and/or unknowns (z,w)");
  /* Output info. : 0: ok -  >0: error (which depends on the chosen solver) */
  int info = -1;

  switch(options->solverId)
  {
  case SICONOS_MCP_OLD_FB: // Fischer-Burmeister/Newton
    mcp_old_FischerBurmeister(problem, z, w, &info, options);
    break;

  default:
    fprintf(stderr, "mcp_old_driver error: unknown solver id: %d\n", options->solverId);
    exit(EXIT_FAILURE);

  }

  return info;
}

void mcp_old_driver_init(MixedComplementarityProblem_old* problem, SolverOptions* options)
{
  switch(options->solverId)
  {
  case SICONOS_MCP_OLD_FB :
    mcp_old_FischerBurmeister_init(problem, options) ;
    break ;
  default :
    fprintf(stderr, "mcp_old_driver_init error: unknown solver id: %d\n", options->solverId);
    exit(EXIT_FAILURE);
  }

}

void mcp_old_driver_reset(MixedComplementarityProblem_old* problem, SolverOptions* options)
{
  switch(options->solverId)
  {
  case SICONOS_MCP_OLD_FB :
    mcp_old_FischerBurmeister_reset(problem, options) ;
    break ;
  default :
    fprintf(stderr, "mcp_old_driver_init error: unknown solver id: %d\n", options->solverId);
    exit(EXIT_FAILURE);
  }
}
