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
#ifndef FRICTIONCONTACTPROBLEM_H
#define FRICTIONCONTACTPROBLEM_H

/*! \page fcProblem Friction-contact problems (2 or 3-dimensional)
  \section fcIntro The problem
  Find \f$(reaction,velocity)\f$ such that:\n

  \f$
  \left\lbrace
  \begin{array}{l}
  M \ reaction + q = velocity \\
  0 \le reaction_n \perp velocity_n \ge 0 \\
  -velocity_t \in \partial\psi_{D_(\mu reaction_n)}(reaction_t)\\
  D_(\mu reaction_n) = \{ reaction_t \mid  \|reaction_t\| \leq \mu reaction_n  \}
  \end{array}
  \right.
  \f$

  \f$ reaction, velocity, q\f$ are vectors of size n and \f$ M \f$ is a nXn matrix, with \f$ n = 2*nc or 3*nc \f$, nc being the number of contacts. \n
  \f$ reaction_n\f$ represents the normal part of the reaction while \f$ reaction_t\f$ is its tangential part.

  \f$ \mu \f$ is the friction coefficient (it may be different for each contact).


  \bf Note: for 2-dimensional friction-contact problems, two formulation are available:\n
  <em>primal problem:</em>\n
  \f$
  \left\lbrace
  \begin{array}{l}
  velocity - M reaction = q \\
  0 \le reaction_n \perp velocity_n \ge 0\\
  -velocity_t \in \partial\psi_{[-\mu reaction_n, \mu reaction_n]}(reaction_t)\\
  \end{array}
  \right.
  \f$

  or \n

  <em>dual problem:</em>\n
  \f$
  \left\lbrace
  \begin{array}{l}
  velocity - M reaction = q \\
  0 \le reaction_n \perp velocity_n \ge 0 \\
  -reaction_t \in \partial\psi_{[-\mu velocity_n, \mu velocity_n]}(velocity_t)\\
  \end{array}
  \right.
  \f$


  \section fc3DSolversList Available solvers for Friction Contact 3D
  Use the generic function frictionContact3D_driver() to call one the the specific solvers listed below:

  - frictionContact3D_nsgs() : non-smooth Gauss-Seidel solver

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

