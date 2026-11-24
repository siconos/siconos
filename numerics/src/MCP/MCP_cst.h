// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef MCP_CST_H
#define MCP_CST_H
/*!\file MCP_cst.h
  \brief Constants to define the list of available MCP solvers.

  \rst
  See the detailed documentation in :ref:`mcp_solvers`
  \endrst
*/

/**\enum MCP_SOLVER

   Each SICONOS_MCP_XXX refers to number of the solver XXX for MCP.
 */
enum MCP_SOLVER {

  SICONOS_MCP_NEWTON_FB_FBLSA = 700,
  SICONOS_MCP_NEWTON_MIN_FBLSA = 701,

  /* old interface for MCP */
  SICONOS_MCP_OLD_FB = 710,
};

#endif
