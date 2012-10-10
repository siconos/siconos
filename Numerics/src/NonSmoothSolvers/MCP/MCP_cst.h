/* Siconos-Numerics, Copyright INRIA 2005-2011.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

#ifndef MCP_CST_H
#define MCP_CST_H
/*!\file mcp_cst.h
  \brief Constants to define the list of available MCP solvers. See the solver list \ref mcpSolversList
*/
/**\enum MCP_SOLVER
   \brief Each SICONOS_MCP_XXX refers to number of the solver XXX for MCP. See the solver list \ref mcpSolversList
 */
enum MCP_SOLVER
{
  SICONOS_MCP_FB = 700
};


extern char SICONOS_MCP_FB_STR[];

#endif
