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
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#ifndef FRICTIONCONTACTPROBLEM_H
#define FRICTIONCONTACTPROBLEM_H

/*! \page fcProblem Friction-contact problems (2 or 3-dimensional)
 *
 * \section fcIntro The problem
 *  Given
 * <ul>
 *   <li> a symmetric positive semi--definite  matrix \f${M} \in {{\mathrm{I\!R}}}^{m \times m} \f$ </li>
 *   <li> a vector \f$ {q} \in {{\mathrm{I\!R}}}^m\f$</li>
 *   <li> a vector of coefficients of friction \f$\mu \in{{\mathrm{I\!R}}}^{n_c}\f$</li>
 *</ul>
 * the (reduced or dual) frictional contact problem  is to find two vectors \f$u\in{{\mathrm{I\!R}}}^m\f$,
 * the relative local velocity and \f$r\in {{\mathrm{I\!R}}}^m\f$,
 * the contact forces denoted by \f$\mathrm{FC}(W,q,\mu)\f$  such that
 * \f{eqnarray*}{
 * \begin{cases}
 *   u = M r + q \\
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
 *
 * The Coulomb's friction law with the Signorini condition for the unilateral contact written in terms of second order complementarity condition
 * \f{eqnarray}{
 *    C^\star_{\mu} \ni {\hat u} \perp r \in C_{\mu}
 * \f}
 * can be interpret usually as
 *
 * \f{eqnarray}{
 * \begin{cases}
 *  0 \leq u_{N} \perp r_N \geq 0  \quad\quad\text{ Signorini condition}\\
 *  u_T = 0 \Rightarrow \|r_T\| \leq \mu |r_n|  \quad\quad\text{ Sticking mode} \\
 *  u_T \neq 0 \Rightarrow r_T = - \mu |r_n| \frac{u_T }{\|u_T\|}  \quad\quad\text{ Sliding mode}
 * \end{cases}
 * \f}


  \section fc3DSolversList Available solvers for Friction Contact 3D
  Use the generic function frictionContact3D_driver() to call one the the specific solvers listed below:

  <ul>

  <li> frictionContact3D_nsgs() : non-smooth Gauss-Seidel solver: SolverId : </li>
  <li> frictionContact3D_nsgs_velocity() : non-smooth Gauss-Seidel solver  based on velocity updates</li>
  <li> frictionContact3D_proximal() : Proximal point solver for friction-contact 3D problem</li>
  <li> frictionContact3D_TrescaFixedPoint() </li>
  <li> frictionContact3D_ProjectedGradientOnCylinder()</li>
  <li> frictionContact3D_DeSaxceFixedPoint()</li>
  <li> frictionContact3D_ExtraGradient()</li>
  <li> frictionContact3D_HyperplaneProjection()</li>
  </li>  frictionContact3D_AlartCurnierNewton()</li>
  </ul>
  (see the functions/solvers list in FrictionContact3D_Solvers.h)

  \section fc3DParam Required and optional parameters
  FrictionContact3D problems needs some specific parameters, given to the FrictionContact3D_driver() function thanks to a SolverOptions structure. \n
  They are:\n
     - the name of the solver (ex: NSGS), used to switch to the right solver function
     - iparam[0]: max. number of iterations allowed
     - iparam[1]:
     - dparam[0]: tolerance
     - isStorageSparse: 1 if a SparseBlockStructuredMatrix is used for M, else 0 (double* storage)

  \section fc2DSolversList Available solvers for Friction Contact 2D

  - pfc_2D_latin(), latin solver
  - pfc_2D_nlgs(), Non Linear Gauss Seidel solver
  - pfc_2D_cpg(), conjugated projected gradient solver
  - dfc_2D_latin(), latin solver

*/

/*!\file FrictionContactProblem.h
  \brief Definition of a structure to handle with friction-contact (2D or 3D) problems.
  \author Franck Perignon.
*/

#include "NumericsMatrix.h"

/** The structure that defines a Friction-Contact (3D or 2D) problem
    \param dimension dimension of the contact space (3D or 2D )
    \param numberOfContacts the number of contacts
    \param M matrix (n X n, with n = 2 or 3*numberOfContacts)
    \param q vector (n)
    \param mu vector of friction coefficients (size: numberOfContacts)
*/
typedef struct
{
  int dimension;
  int numberOfContacts;
  NumericsMatrix* M;
  double* q;
  double* mu;
} FrictionContactProblem;


#ifdef __cplusplus
extern "C"
{
#endif
  int frictionContact_printInFile(FrictionContactProblem*  problem, FILE* file);

  int frictionContact_newFromFile(FrictionContactProblem*  problem, FILE* file);

  void freeFrictionContactProblem(FrictionContactProblem* problem);

#ifdef __cplusplus
}
#endif

#endif

