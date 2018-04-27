/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#ifndef MCP_SOLVERS_H
#define MCP_SOLVERS_H

/*!\file MCP_Solvers.h
  \brief List of all the available solvers for the resolution of Mixed Complementarity Problems.\n

  \author siconos-team@lists.gforge.inria.fr
*/

#include "MixedComplementarityProblem.h"
#include "SolverOptions.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif
  /** Initialisation of the MCP solver (set problem, allocate working memory and so on. This routine must be called before any attempt to run the mcp_driver.
      \param[in] problem the description of the MCP
      \param[in] options for the solver
  */
  void mcp_driver_init(MixedComplementarityProblem * problem, SolverOptions* options);

  /** Reset of the MCP solver
     \param[in] problem the description of the MCP
     \param[in] options for the solver
  */
  void mcp_driver_reset(MixedComplementarityProblem * problem, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for MixedLinearComplementarity
      \param problem  the pointer to the array of options to set.
      \param pOptions the pointer to the SolverOptions stucture.
  */
  int mixedComplementarity_setDefaultSolverOptions(MixedComplementarityProblem* problem, SolverOptions* pOptions);

  /** set the default solver parameters and perform memory allocation for MixedLinearComplementarity
      \param problem  the pointer to the array of options to set.
      \param pOptions the pointer to the SolverOptions stucture.
  */
  void  mixedComplementarity_default_setDefaultSolverOptions(MixedComplementarityProblem* problem, SolverOptions* pOptions);

  /** Fischer Burmeister solver
      \param[in] problem a structure which represents the MCP
      \param[in,out] z a m+n-vector, initial solution + returns the solution of the problem.
      \param[out] w not used
      \param[out] info termination value:   0 if success else >0.
      \param[in,out] options structure used to define the solver and its parameters.
      \author Franck Pérignon
  */
  void mcp_FischerBurmeister(MixedComplementarityProblem* problem, double* z, double* w, int* info, SolverOptions* options);

  /** Solver based on Fischer-Burmeister reformulation and line search (VFBLSA
   * in Facchinei--Pang 2003 p. 834)
      \param[in] problem a structure which represents the MCP
      \param[in,out] z a n1+n2-vector, initial solution + returns the solution of the problem.
      \param[out] Fmcp n1+n2-vector which contains the value of Fmcp(z) = (G(z), H(z))
      \param[out] info termination value:   0 if success else >0.
      \param[in,out] options structure used to define the solver and its parameters.
      \author Olivier Huber
  */
  void mcp_newton_FBLSA(MixedComplementarityProblem2* problem, double* z, double* Fmcp, int* info, SolverOptions* options);

  /** Solver based on Fischer-Burmeister reformulation and line search. The
   * descent direction is found using a min reformulation (minFBLSA in
   * Facchinei--Pang 2003 p. 855)
      \param[in] problem a structure which represents the MCP
      \param[in,out] z a n1+n2-vector, initial solution + returns the solution of the problem.
      \param[out] Fmcp n1+n2-vector which contains the value of Fmcp(z) = (G(z), H(z))
      \param[out] info termination value:   0 if success else >0.
      \param[in,out] options structure used to define the solver and its parameters.
      \author Olivier Huber
  */
  void mcp_newton_minFBLSA(MixedComplementarityProblem2* problem, double* z, double* Fmcp, int* info, SolverOptions* options);

  /** set the default solver parameters and perform memory allocation for MixedLinearComplementarity
      \param problem  the pointer to the array of options to set.
      \param pSolver the pointer to the SolverOptions stucture.
  */
  int mixedComplementarity_FB_setDefaultSolverOptions(MixedComplementarityProblem* problem, SolverOptions* pSolver);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif


#endif
