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

#ifndef MCP_CST_H
#define MCP_CST_H
/*!\file MCP_cst.h
  \brief Constants to define the list of available MCP solvers. See the solver list \ref mcpSolversList
*/
/**\enum MCP_SOLVER
   Each SICONOS_MCP_XXX refers to number of the solver XXX for MCP. See the solver list \ref mcpSolversList
 */
enum MCP_SOLVER
{
  SICONOS_MCP_FB = 700,
  SICONOS_MCP_NEWTON_FBLSA = 701,
  SICONOS_MCP_NEWTON_MINFBLSA = 702
};


extern const char* const SICONOS_MCP_FB_STR;
extern const char* const SICONOS_MCP_NEWTON_FBLSA_STR;
extern const char* const SICONOS_MCP_NEWTON_MINFBLSA_STR;

#endif
