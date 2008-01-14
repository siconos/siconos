/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2007.
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
#ifndef FRICTIONCONTACT3DSOLVERS_H
#define FRICTIONCONTACT3DSOLVERS_H

#include "NSSTools.h" /* for SparseBlockStructuredMatrix structure */
#include "FrictionContact3D_AlartCurnier.h"
#include "FrictionContact3D_projection.h"
#include "FrictionContact3D_Newton.h"
#include "FrictionContact3D2NCP_Glocker.h"

/*! \page fc3DSolvers Friction-contact problems (3-dimensional)
  \section fc3Dintro The problem
  Find \f$(reaction,velocity)\f$ such that:\n

  \f$
  \left\lbrace
  \begin{array}{l}
  M \ reaction + q = velocity \\
  0 \le reaction_n \perp velocity_n \ge 0 \\
  -velocity_t \in \partial\psi_{[-\mu.reaction_n, \mu.reaction_n]}(reaction_t)\\
  \end{array}
  \right.
  \f$

  \f$ reaction, velocity, q\f$ are vectors of size n and \f$ M \f$ is a nXn matrix, with \f$ n = 3*nc \f$, nc being the number of contacts. \n
  \f$ reaction_n\f$ represents the normal part of the reaction while \f$ reaction_t\f$ is its tangential part.

  \f$ \mu \f$ is the friction coefficient (may be different for each contact).

  \section fc3DSolversList Available solvers
  Use the generic function frictionContact3DSolvers(), to call one the the specific solvers listed below:

  - frictionContact3D_nsgs() : non-smooth Gauss-Seidel solver

  The structure method, argument of frictionContact3DSolvers(), is used to give the name and parameters of the required solver.

  (see the functions/solvers list in FrictionContact3DSolvers.h)


*/

/*!\file FrictionContact3D_Solvers.h
  Subroutines for the resolution of contact problems with friction (3-dimensional case).\n

  \author Mathieu Renouf, Franck Perignon.

*/

/** pointer to function used to call local solver */
typedef void (*SolverPtr)(int, int, double*, int*, double*);

/** pointer to function used to post-processed results after a call to the (local) solver */
typedef void (*PostSolverPtr)(int, double*);

/** pointer to function used to update velocity and compute error */
typedef void (*ComputeErrorPtr)(int, double*, double*, double*);

/** pointer to function used to free memory for objects used in solvers */
typedef void (*FreeSolverPtr)();

#ifdef __cplusplus
extern "C" {
#endif

  /** Non-Smooth Gauss Seidel solver for friction-contact 3D problem
      \param nc, number of contacts (dim of the problem n = 3*nc)
      \param M global matrix (n*n)
      \param q global vector (n)
      \param reaction global vector (n), in-out parameter
      \param velocity global vector (n), in-out parameter
      \param mu, vector of the friction coefficients (size nc)
      \param information about solver result
      \param int vector of parameters (max. iteration number ...)
      \param double vector of parameters (tolerance ...)
  */
  void frictionContact3D_nsgs(int, double*, double*, double*, double*, double*, int *, int *, double*);

  /** Non-Smooth Gauss Seidel solver for friction-contact 3D problem
      \param nc, number of contacts (dim of the problem n = 3*nc)
      \param M global matrix (n*n)
      \param q global vector (n)
      \param reaction global vector (n), in-out parameter
      \param velocity global vector (n), in-out parameter
      \param mu, vector of the friction coefficients (size nc)
      \param information about solver result
      \param int vector of parameters (max. iteration number ...)
      \param double vector of parameters (tolerance ...)
  */
  void frictionContact3D_nsgs_SB(int, SparseBlockStructuredMatrix*, double*, double*, double*, double*, int *, int *, double*);

  /** Check for trivial solution in the friction-contact 3D problem
      \param dim of the problem
      \param q global vector (n)
      \param velocity global vector (n), in-out parameter
      \param reaction global vector (n), in-out parameter
      \param int vector of parameters (max. iteration number ...)
      \param double vector of parameters (tolerance ...)
      \return int =0 if a trivial solution has been found, else = -1
  */
  int checkTrivialCase(int, double*, double*, double*, int*, double*);

#ifdef __cplusplus
}
#endif

#endif
