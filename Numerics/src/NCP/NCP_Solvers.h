/* Siconos-Numerics, Copyright INRIA 2005-2012.
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

#ifndef NCP_H
#define NCP_H

/*! \page NCProblem Nonlinear Complementarity Problems (NCP)

  \section ncpIntro The problem
  Find \f$z \in \mathcal{R}^n_+\f$ such that:
  \f{equation*}{
  0 \le z \perp F(z) \ge 0
  \f}

  \section ncpSolvers Available solvers/formulations:
   - ncp_newton_FBLSA() with the FB merit function and a Newton with line-search
   - ncp_newton_minFBLSA() with the min merit function (with the FB as backup) and a Newton with line-search
   - ncp_pathsearch() solver using a path search
   - NCP_Path() Interface to Path (Ferris)
*/

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
   void ncp_pathsearch(NCP_struct* problem, double* z, double* F, int *info , SolverOptions* options);

  /** NCP Solver using the PATH solver
   * \param problem the formalization of the NCP problem
   * \param[in,out] z on input, initial guess; on output the solution
   * \param F the value of the function at the solution
   * \param info 0 if everything worked
   * \param options struct used to specify the solver parameters
   */
   void ncp_path(NCP_struct* problem, double* z, double* F, int *info , SolverOptions* options);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
