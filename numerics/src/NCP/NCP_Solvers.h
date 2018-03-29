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

#ifndef NCP_H
#define NCP_H

/*!\file NCP_Solvers.h
  \brief Functions related to NCP formulation and solvers.
  \author Franck Perignon, Olivier Huber
*/

#include "SparseBlockMatrix.h"
#include "NCP_FixedP.h"

#include "SiconosConfig.h"
#include "SolverOptions.h"
#include "NonlinearComplementarityProblem.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /**
   * This function compute the complementarity error of the NCP: \f$ 0 \leq z \perp F(z) \geq 0\f$.
   * \param n size of the vectors
   * \param[in] z solution
   * \param[in] F value of the function at the solution
   * \param[in] tol tolerance for the error
   * \param[out] err value of the error
   * \return 0 if the solution is good enough, 1 otherwise
   * \author Olivier Huber
  */
  int ncp_compute_error(int n, double* z, double* F, double tol, double *err);

  /** NCP Solver using the FB merit function and a Newton-based method with
   * line-search
   * \param problem the formalization of the NCP problem
   * \param[in,out] z on input, initial guess; on output the solution
   * \param F the value of the function at the solution
   * \param info 0 if everything worked
   * \param options struct used to specify the solver parameters
   */
   void ncp_newton_FBLSA(NonlinearComplementarityProblem* problem, double *z, double* F, int *info, SolverOptions* options);

   /** NCP Solver using the min merit function (+ the FB as backup) and a Newton-based method with
   * line-search
   * \param problem the formalization of the NCP problem
   * \param[in,out] z on input, initial guess; on output the solution
   * \param F the value of the function at the solution
   * \param info 0 if everything worked
   * \param options struct used to specify the solver parameters
   */
   void ncp_newton_minFBLSA(NonlinearComplementarityProblem* problem, double *z, double* F, int *info, SolverOptions* options);


  /** NCP Solver using a path search algorithm, following the work of D. Ralph.
   * M. Ferris, and many other collaborators of the latter.
   * \param problem the formalization of the NCP problem
   * \param[in,out] z on input, initial guess; on output the solution
   * \param F the value of the function at the solution
   * \param info 0 if everything worked
   * \param options struct used to specify the solver parameters
   */
   void ncp_pathsearch(NonlinearComplementarityProblem* problem, double* z, double* F, int *info , SolverOptions* options);

  /** NCP Solver using the PATH solver
   * \param problem the formalization of the NCP problem
   * \param[in,out] z on input, initial guess; on output the solution
   * \param F the value of the function at the solution
   * \param info 0 if everything worked
   * \param options struct used to specify the solver parameters
   */
   void ncp_path(NonlinearComplementarityProblem* problem, double* z, double* F, int *info , SolverOptions* options);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
