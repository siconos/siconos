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
      \param[in] : the description of the MCP
      \param[in] : options for the solver
  */
  void mcp_FB_init(MixedComplementarityProblem * problem, SolverOptions* options);


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
