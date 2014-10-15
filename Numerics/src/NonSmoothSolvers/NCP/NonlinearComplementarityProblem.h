/* Siconos-Numerics, Copyright INRIA 2005-2014.
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
#ifndef NCP_PROBLEM_H
#define NCP_PROBLEM_H

#include "NumericsMatrix.h"

/*!\file NonlinearComplementarityProblem.h
 * \brief data structure to formalize a Nonlinear Complementarity Problem (NCP)
 *
 * \author Olivier Huber
*/

/*! \page NCProblem Nonlinear Complementarity problems (NCP)
  \section ncpIntro  Problem Statement
  Given a sufficiently smooth function \f${F}\colon {{\mathrm{I\!R}}}^{n}  \to {{\mathrm{I\!R}}}^{n}\f$
  The Nonlinear Complementarity Problem (NCP) is to find two vectors \f$(z,w \in {{\mathrm{I\!R}}}^{n})\f$ such that:
  \f{align*}{
  w &= F(z) \\
  0 &\le w \perp z \ge 0
  \f}

  \section ncpSolversList Available solvers:
  - ncp_FBLSA(), nonsmooth Newton method based on Fisher-Burmeister function with a line search.
  - ncp_pathsearch(), a solver based on a path search method
*/

/** type for user defined function used to compute F and its jacobian.
 */
typedef void (*ptrFunctionNCP)(void* env, int n, double* z, double* F);
typedef void (*ptrFunctionJacNCP)(void* env, int n, double* z, NumericsMatrix* jacF);

/** \struct  NonlinearComplementarityProblem NonlinearComplementarityProblem.h
 * The structure that defines a Nonlinear Complementarity Problem (NCP) : Find two vectors \f$(z,w \in {{\mathrm{I\!R}}}^{n})\f$ such that:\n
  \f{align*}{
  w &= F(z) \\
  0 &\le w \perp z \ge 0
  \f}
 */
typedef struct
{
  int n; /**< size of the problem */
  ptrFunctionNCP compute_F; /**< pointer to the function used to compute \f$F(z)\f$ */
  ptrFunctionJacNCP compute_nabla_F; /**< pointer to the function used to compute \f$\nabla_z F(z)\f$ */
  NumericsMatrix* nabla_F; /**< storage for \f$\nabla_z F\f$*/
  void* env; /**< environment for the compute_Fmcp and compute_nabla_F function.
               When called from Python, it contains an object with compute_F and compute_nabla_F as methods.
               When called from C, it can reference a data struct containing variables needed for the computations.*/
} NonlinearComplementarityProblem;

typedef NonlinearComplementarityProblem NCP_struct;

#include "NumericsConfig.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** compute the value of the Fischer-Burmeister function with an NCP
   * \param data_opaque the strcture formalizing the NCP problem
   * \param z current iterate
   * \param F current value of the nonlinear function
   * \param[out] F_FB value of the FB function
   */
  void ncp_FB(void* data_opaque, double* z, double* F, double* F_FB);
void FB_compute_F_ncp(void* data_opaque, double* z, double* Fmcp);
void FB_compute_H_ncp(void* data_opaque, double* z, double* Fmcp, double* workV1, double* workV2, double* H);
void FB_compute_error_ncp(void* data_opaque, double* z, double* w, double* Jac_F_merit, double tol, double* err);


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
