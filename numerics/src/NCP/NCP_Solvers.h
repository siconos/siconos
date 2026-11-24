/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef NCP_H
#define NCP_H

/*!\file NCP_Solvers.h
  \brief Functions related to NCP formulation and solvers.
*/

#include "NonlinearComplementarityProblem.h"  // IWYU pragma: keep

#if defined(__cplusplus)
extern "C" {
#endif

/** This function compute the complementarity error of the NCP
 *  \f$ 0 \leq z \perp F(z) \geq 0 \f$.
 *
 *  \param n size of the vectors
 *  \param[in] z solution
 *  \param[in] F value of the function at the solution
 *  \param[in] tol tolerance for the error
 *  \param[out] err value of the error
 *  \return 0 if the solution is good enough, 1 otherwise
 */
int ncp_compute_error(int n, double* z, double* F, double tol, double* err);

/** NCP Solver using the FB merit function and a Newton-based method with
 *  line-search
 *
 *  \param problem the formalization of the NCP problem
 *  \param[in,out] z on input, initial guess; on output the solution
 *  \param F the value of the function at the solution
 *  \param info 0 if everything worked
 *  \param options struct used to specify the solver parameters
 */
void ncp_newton_FBLSA(NonlinearComplementarityProblem* problem, double* z, double* F,
                      int* info, SolverOptions* options);

/** NCP Solver using the min merit function (+ the FB as backup) and a Newton-based method with
 *  line-search
 *
 *  \param problem the formalization of the NCP problem
 *  \param[in,out] z on input, initial guess; on output the solution
 *  \param F the value of the function at the solution
 *  \param info 0 if everything worked
 *  \param options struct used to specify the solver parameters
 */
void ncp_newton_minFBLSA(NonlinearComplementarityProblem* problem, double* z, double* F,
                         int* info, SolverOptions* options);

/** NCP Solver using a path search algorithm, following the work of D. Ralph.
 *  M. Ferris, and many other collaborators of the latter.
 *
 *  \param problem the formalization of the NCP problem
 *  \param[in,out] z on input, initial guess; on output the solution
 *  \param F the value of the function at the solution
 *  \param info 0 if everything worked
 *  \param options struct used to specify the solver parameters
 */
void ncp_pathsearch(NonlinearComplementarityProblem* problem, double* z, double* F, int* info,
                    SolverOptions* options);

/** NCP Solver using the PATH solver
 *
 *  \param problem the formalization of the NCP problem
 *  \param[in,out] z on input, initial guess; on output the solution
 *  \param F the value of the function at the solution
 *  \param info 0 if everything worked
 *  \param options struct used to specify the solver parameters
 */
void ncp_path(NonlinearComplementarityProblem* problem, double* z, double* F, int* info,
              SolverOptions* options);

#if defined(__cplusplus)
}
#endif

#endif
