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
#ifndef GLOBALFRICTIONCONTACTPROBLEM_H
#define GLOBALFRICTIONCONTACTPROBLEM_H

/*! \page globalFcProblem Global-Friction-contact problems (2 or 3-dimensional)
 * \section pfcIntro Problem statement.
 *  Given
 * <ul>
 *   <li> a symmetric positive semi--definite  matrix \f${M} \in {{\mathrm{I\!R}}}^{n \times n} \f$ </li>
 *   <li> a matrix \f${H} \in {{\mathrm{I\!R}}}^{n \times {d\, n_c}} \f$ </li>
 *   <li> a vector \f$ {q} \in {{\mathrm{I\!R}}}^n\f$</li>
 *   <li> a vector \f$ {b} \in {{\mathrm{I\!R}}}^{d\, n_c}\f$</li>
 *   <li> a vector of coefficients of friction \f$\mu \in{{\mathrm{I\!R}}}^{n_c}\f$</li>
 *</ul>
 * the (global or global) frictional contact problem  is to find three vectors \f$v\in{{\mathrm{I\!R}}}^n\f$,
 * the (global) velocity \f$u\in{{\mathrm{I\!R}}}^{d\,n_c}\f$,
 * the relative local velocity and \f$r\in {{\mathrm{I\!R}}}^{d,n_c}\f$,
 * the contact forces denoted by \f$\mathrm{PFC}(M,H,q,b,\mu)\f$  such that
 * \f{eqnarray*}{
 * \begin{cases}
 *  M v =  q +  H r \\
 *   u = H^\top v + b \\
 *    \hat u = u +\left[
 *      \left[\begin{array}{c}
 *          \mu^\alpha \|u^\alpha_{T}\|\\
 *         0 \\
 *         0
 *        \end{array}\right]^T, \alpha = 1 \ldots n_c
 *    \right]^T \\ \                                \
 *    C^\star_{\mu} \ni {\hat u} \perp r \in C_{\mu}
 * \end{cases}
 * \f}
 * and the set \f$C^{\alpha,\star}_{\mu^\alpha}\f$ is its dual.
 * The modified local velocity \f$\widehat u \f$ is not considered as an unknown since it can obtained uniquely
 * from the local velocity \f$u\f$.
 * Coulomb's friction law with Signorini's condition for the unilateral contact written in terms
 * of second order complementarity condition
 * \f{eqnarray}{
 *    C^\star_{\mu} \ni {\hat u} \perp r \in C_{\mu}
 * \f}
 * can be interpreted in a more usual form
 *
 * \f{eqnarray}{
 * \begin{cases}
 *  0 \leq u_{N} \perp r_N \geq 0  \quad\quad\text{ Signorini condition}\\
 *  u_T = 0 \Rightarrow \|r_T\| \leq \mu |r_n|  \quad\quad\text{ Sticking mode} \\
 *  u_T \neq 0 \Rightarrow r_T = - \mu |r_n| \frac{u_T }{\|u_T\|}  \quad\quad\text{ Sliding mode}
 * \end{cases}
 * \f}
 *
 * This problem models any instance of discretized frictional contact problem obtained from
 * <ul>
 * <li>the time-discretization  of dynamics contact problems with event-capturing of event-tracking schemes, </li>
 * <li>the time-discretization  of quasi-static contact problems, </li>
 * <li>the modeling  of static contact problems. In this last case, \f$u\f$ plays the role of the relative displacement at contact  </li>
 * </ul>
 *
 * The problem is stored and given to the solver in Siconos/Numerics thanks to
 *  a C structure GlobalFrictionContactProblem .
 *
 *
 * \section pfc3DSolversList Available solvers for Friction Contact 3D
 * Use the generic function globalFrictionContact3D_driver() to call one the the specific solvers listed below:
 * <ul>
 *  <li> globalfrictionContact3D_nsgs() : non-smooth Gauss-Seidel solver </li>
 * </ul>
 * (see the functions/solvers list in GlobalFrictionContact3D_Solvers.h)
 *
 * \section pfc3DParam Required and optional parameters
 * GlobalFrictionContact3D problems needs some specific parameters, given to the GlobalFrictionContact3D_driver() function thanks to a SolverOptions structure. \n
 *
 *
 */

/*!\file GlobalFrictionContactProblem.h
  \brief Definition of a structure to handle with friction-contact (2D or 3D) problems.
  \author Vincent Acary.
*/

#include "NumericsMatrix.h"

/** \struct GlobalFrictionContactProblem GlobalFrictionContactProblem.h
 * The structure that defines a Friction-Contact (3D or 2D)problem \f$\mathrm{PFC}(M,H,q,b,\mu)\f$  such that
 * \f{eqnarray*}{
 * \begin{cases}
 *  M v =  q +  H r \\
 *  u = H^\top v + b \\
 *    \hat u = u +\left[
 *      \left[\begin{array}{c}
 *          \mu^\alpha \|u^\alpha_{T}\|\\
 *         0 \\
 *         0
 *        \end{array}\right]^T, \alpha = 1 \ldots n_c
 *    \right]^T \\ \                                \
 *    C^\star_{\mu} \ni {\hat u} \perp r \in C_{\mu}
 * \end{cases}
 * \f}
 * and the set \f$C^{\alpha,\star}_{\mu^\alpha}\f$ is its dual.

*/
typedef struct
{
  /** dimension \f$d=2\f$ or \f$d=3\f$ of the contact space (3D or 2D ) */
  int dimension;
  /** the number of contacts \f$ n_c \f$ */
  int numberOfContacts;
  /** \f${M} \in {{\mathrm{I\!R}}}^{n \times n} \f$,
      a matrix with \f$ n\f$ stored in NumericsMatrix structure */
  NumericsMatrix* M;
  /**  \f${H} \in {{\mathrm{I\!R}}}^{n \times m} \f$,
      a matrix with \f$ m = d  n_c\f$ stored in NumericsMatrix structure */
  NumericsMatrix* H;
  /** \f${q} \in {{\mathrm{I\!R}}}^{n} \f$ */
  double* q;
  /** \f${b} \in {{\mathrm{I\!R}}}^{m} \f$ */
  double* b;
  /** mu \f${\mu} \in {{\mathrm{I\!R}}}^{n_c} \f$, vector of friction coefficients
      (\f$ n_c =\f$ numberOfContacts) */
  double* mu;
} GlobalFrictionContactProblem;

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif
  int globalFrictionContact_printInFile(GlobalFrictionContactProblem*  problem, FILE* file);

  int globalFrictionContact_newFromFile(GlobalFrictionContactProblem*  problem, FILE* file);

  void freeGlobalFrictionContact_problem(GlobalFrictionContactProblem* problem);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif
#endif
