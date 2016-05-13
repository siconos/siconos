/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.

 * Copyright 2016 INRIA.

 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at

 * http://www.apache.org/licenses/LICENSE-2.0

 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

#ifndef MCPFB_H
#define MCPFB_H

#include "FischerBurmeister.h"

/*!\file MCP_FischerBurmeister.h

  routines required for the MCP solver based on Fischer-Burmeister functions and Semi-Smooth Newton solver.

*/

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif


  /** Initialisation of the MCP Fischer solver (set problem, allocate working memory and so on. This routine must be called before any attempt to run the mcp_driver.
   *   \param[in] problem  description of the MCP
   *   \param[in] options for the solver
   */
  void mcp_FischerBurmeister_init(MixedComplementarityProblem * problem, SolverOptions* options);

  /** Reset of the MCP Fischer solver (free local variable)
   *  \param[in] problem  description of the MCP
   *  \param[in] options for the solver
   */
  void mcp_FischerBurmeister_reset(MixedComplementarityProblem * problem, SolverOptions* options);


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
