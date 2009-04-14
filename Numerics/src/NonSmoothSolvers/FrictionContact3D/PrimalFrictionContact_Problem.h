/* Siconos-Numerics version 3.0.0, Copyright INRIA 2005-2008.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */
#ifndef PRIMALFRICTIONCONTACTPROBLEM_H
#define PRIMALFRICTIONCONTACTPROBLEM_H

/*! \page primalFcProblem Primal-Friction-contact problems (2 or 3-dimensional)
 * \section pfcIntro The problem
 * Find \f$(reaction,velocity)\f$ such that:\n
 * \f$
 \left\lbrace
  \begin{array}{l}
  M globalVelocity =  q +  H reaction \\
  velocity = H^T globalVelocity + b\\
  K \ni reaction_n \perp velocity_n \in K^* \\
  \end{array}
  \right.
  \f$\n
  * where
  \f$
  \left\lbrace
  \begin{array}{l}
  K = \{reaction, \|reaction_t\| \leq \mu reaction_n \}
  \end{array}
  \right.
  \f$
  is the Coulomb's Cone \n
    * and with:
    *    - \f$globalVelocity \in R^{n} \f$  the global unknown,
    *    - \f$M \in R^{n \times n } \f$  and \f$q \in R^{n} \f$
    *    - \f$velocity \in R^{m} \f$  and \f$reaction \in R^{m} \f$ the local unknowns,
    *    - \f$b \in R^{m} \f$ is the modified local velocity (\f$ e U_{N,k}\f$)
    *    - \f$M \in R^{n \times n } \f$  and \f$q \in R^{n} \f$
    *    - \f$H \in R^{n \times m } \f$
    \f$ reaction_n\f$ represents the normal part of the reaction while \f$ reaction_t\f$ is its tangential part.

    \f$ \mu \f$ is the friction coefficient (it may be different for each contact).




  \section pfc3DSolversList Available solvers for Friction Contact 3D
  Use the generic function primalFrictionContact3D_driver() to call one the the specific solvers listed below:

  - primalfrictionContact3D_nsgs() : non-smooth Gauss-Seidel solver

  (see the functions/solvers list in PrimalFrictionContact3D_Solvers.h)

  \section pfc3DParam Required and optional parameters
  PrimalFrictionContact3D problems needs some specific parameters, given to the PrimalFrictionContact3D_driver() function thanks to a Solver_Options structure. \n


*/

/*!\file PrimalFrictionContact_Problem.h
  \brief Definition of a structure to handle with friction-contact (2D or 3D) problems.
  \author Vincent Acary.
*/

#include "NumericsMatrix.h"

/** The structure that defines a Friction-Contact (3D or 2D)problem
    \param numberOfContacts, the number of contacts
    \param M matrix (n X n, with n = 2 or 3*numberOfContacts)
    \param H matrix (n X m, with n = 2 or 3*numberOfContacts)
    \param q vector (n)
    \param b vector (m)
    \param mu, vector of friction coefficients (size: numberOfContacts)
    \param isComplete, equal to 0 if some information is missing or wrong for the problem (M or q = NULL, inconsistent sizes), else equal to 1.
*/
typedef struct
{
  int numberOfContacts;
  NumericsMatrix* M;
  NumericsMatrix* H;
  double* q;
  double* b;
  double* mu;
  int isComplete;
} PrimalFrictionContact_Problem;

#endif
