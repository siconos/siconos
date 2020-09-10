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
#ifndef MCP_SOLVERS_H
#define MCP_SOLVERS_H

/*!\file MCP_Solvers.h
  \brief List of all the available solvers for the resolution of Mixed Complementarity Problems.

  \rst
  See the detailed documentation in :ref:`mcp_solvers`
  \endrst

*/

#include "MixedComplementarityProblem.h"
#include "SolverOptions.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Solver based on Fischer-Burmeister reformulation and line search (VFBLSA
   * in Facchinei--Pang 2003 p. 834)
      \param[in] problem a structure which represents the MCP
      \param[in,out] z a n1+n2-vector, initial solution + returns the solution of the problem.
      \param[out] Fmcp n1+n2-vector which contains the value of Fmcp(z) = (G(z), H(z))
      \param[out] info termination value:   0 if success else >0.
      \param[in,out] options structure used to define the solver and its parameters.
  */
  void mcp_newton_FB_FBLSA(MixedComplementarityProblem* problem, double* z, double* Fmcp, int* info, SolverOptions* options);
  /** Solver based on Fischer-Burmeister reformulation and line search. The
   * descent direction is found using a min reformulation (minFBLSA in
   * Facchinei--Pang 2003 p. 855)
      \param[in] problem a structure which represents the MCP
      \param[in,out] z a n1+n2-vector, initial solution + returns the solution of the problem.
      \param[out] Fmcp n1+n2-vector which contains the value of Fmcp(z) = (G(z), H(z))
      \param[out] info termination value:   0 if success else >0.
      \param[in,out] options structure used to define the solver and its parameters.
  */
  void mcp_newton_min_FBLSA(MixedComplementarityProblem* problem, double* z, double* Fmcp, int* info, SolverOptions* options);

  int mcp_compute_error(MixedComplementarityProblem* problem, double *z , double *w,  double * error);
  /** Initialisation of the MCP solver (set problem, allocate working memory and so on. This routine must be called before any attempt to run the mcp_old_driver.
      \param[in] problem the description of the MCP
      \param[in] options for the solver
  */
  
  void mcp_old_driver_init(MixedComplementarityProblem_old * problem, SolverOptions* options);

  /** Reset of the MCP solver
     \param[in] problem the description of the MCP
     \param[in] options for the solver
  */
  void mcp_old_driver_reset(MixedComplementarityProblem_old * problem, SolverOptions* options);

  /** Fischer Burmeister solver
      \param[in] problem a structure which represents the MCP
      \param[in,out] z a m+n-vector, initial solution + returns the solution of the problem.
      \param[out] w not used
      \param[out] info termination value:   0 if success else >0.
      \param[in,out] options structure used to define the solver and its parameters.
  */
  void mcp_old_FischerBurmeister(MixedComplementarityProblem_old* problem, double* z, double* w, int* info, SolverOptions* options);


  /** Initialisation of the MCP Fischer solver (set problem, allocate working memory and so on. This routine must be called before any attempt to run the mcp_old_driver.
   *   \param[in] problem  description of the MCP
   *   \param[in] options for the solver
   */
  void mcp_old_FischerBurmeister_init(MixedComplementarityProblem_old * problem, SolverOptions* options);

  /** Reset of the MCP Fischer solver (free local variable)
   *  \param[in] problem  description of the MCP
   *  \param[in] options for the solver
   */
  void mcp_old_FischerBurmeister_reset(MixedComplementarityProblem_old * problem, SolverOptions* options);
  
  int mcp_old_compute_error(MixedComplementarityProblem_old* problem, double *z , double *w,  double * error);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif


#endif
