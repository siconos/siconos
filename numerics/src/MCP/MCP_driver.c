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

//#include "SolverOptions.h"
//#include "MixedComplementarityProblem.h"
#include <assert.h>

#include "MCP_Solvers.h"
#include "MCP_cst.h"
#include "MCP_FischerBurmeister.h"

#include "NonSmoothDrivers.h"
#include "numerics_verbose.h"

const char* const SICONOS_MCP_FB_STR = "NewtonFB";
const char* const SICONOS_MCP_NEWTON_FBLSA_STR = "Newton FBLSA";
const char* const SICONOS_MCP_NEWTON_MINFBLSA_STR = "Newton minFBLSA";

int mcp_driver2(MixedComplementarityProblem2* problem, double *z , double *Fmcp, SolverOptions* options)
{
  assert(options != NULL);
  /* Checks inputs */
  assert(problem != NULL);
  assert(z != NULL);
  assert(Fmcp != NULL);
  /* Output info. : 0: ok -  >0: error (which depends on the chosen solver) */
  int info = -1;

  switch (options->solverId)
  {
  case SICONOS_MCP_NEWTON_FBLSA: // Fischer-Burmeister/Newton -- new version
    mcp_newton_FBLSA(problem, z, Fmcp, &info, options);
    break;

  case SICONOS_MCP_NEWTON_MINFBLSA: // Fischer-Burmeister/Newton + min descent direction
    mcp_newton_minFBLSA(problem, z, Fmcp, &info, options);
    break;

  default:
    fprintf(stderr, "mcp_driver error: unknown solver id: %d\n", options->solverId);
    exit(EXIT_FAILURE);
  }

  return info;
}

int mcp_driver(MixedComplementarityProblem* problem, double *z , double *w, SolverOptions* options)
{
  if (options == NULL)
    numerics_error("mcp_driver ", "null input for solver options.\n");

  /* Checks inputs */
  if (problem == NULL || z == NULL || w == NULL)
    numerics_error("mcp_driver", "null input for MixedComplementarityProblem and/or unknowns (z,w)");
  /* Output info. : 0: ok -  >0: error (which depends on the chosen solver) */
  int info = -1;

  switch (options->solverId)
  {
  case SICONOS_MCP_FB: // Fischer-Burmeister/Newton
    mcp_FischerBurmeister(problem, z, w, &info, options);
    break;

  default:
    fprintf(stderr, "mcp_driver error: unknown solver id: %d\n", options->solverId);
    exit(EXIT_FAILURE);

  }

  return info;
}

void mcp_driver_init(MixedComplementarityProblem* problem, SolverOptions* options)
{
  switch (options->solverId)
  {
  case SICONOS_MCP_FB :
    mcp_FischerBurmeister_init(problem, options) ;
    break ;
  default :
    fprintf(stderr, "mcp_driver_init error: unknown solver id: %d\n", options->solverId);
    exit(EXIT_FAILURE);
  }

}

void mcp_driver_reset(MixedComplementarityProblem* problem, SolverOptions* options)
{
  switch (options->solverId)
  {
  case SICONOS_MCP_FB :
    mcp_FischerBurmeister_reset(problem, options) ;
    break ;
  default :
    fprintf(stderr, "mcp_driver_init error: unknown solver id: %d\n", options->solverId);
    exit(EXIT_FAILURE);
  }

}
