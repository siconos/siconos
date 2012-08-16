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
#ifndef MCP_SOLVERS_H
#define MCP_SOLVERS_H

/*!\file MCP_Solvers.h
  \brief List of all the available solvers for the resolution of Mixed Complementarity Problems.\n

  \author siconos-team@lists.gforge.inria.fr
*/

/*! \page MCPSolvers MixedComplementarity Problems Solvers

  \section mcp_FB  semi-smooth Newton/Fisher-Burmeister solver.
  a nonsmooth Newton method based based on the Fischer-Bursmeister convex function

  function: mcp_FB() \n
  parameters:
  - iparam[0] (in): maximum number of iterations allowed
  - iparam[1] (out): number of iterations processed
  - dparam[0] (in): tolerance
  - dparam[1] (out): resulting error

 */

#include "MixedComplementarityProblem.h"
#include "SolverOptions.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Fischer Burmeister solver
      \param[in] problem a structure which represents the MCP
      \param[in,out] z a m+n-vector, initial solution + returns the solution of the problem.
      \param[out] w a m+n-vector, solution of the problem.
      \param[out] info termination value:   0 if success else >0.
      \param[in,out] options structure used to define the solver and its parameters.
      \author Franck Pérignon
  */
  void mcp_FB(MixedComplementarityProblem* problem, double* z, double* w, int* info, SolverOptions* options);


  /** Initialisation of the MCP solver (set problem, allocate working memory and so on. This routine must be called before any attempt to run the mcp_driver.
      \param[in] : the description of the MCP
      \param[in] : options for the solver
  */
  void mcp_driver_init(MixedComplementarityProblem * problem, SolverOptions* options);


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif


#endif
