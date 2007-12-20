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

#include "NSSpack.h" /* for method structure */
#include "NSSTools.h" /* for SparseBlockStructuredMatrix structure */

/*!\file FrictionContact3DSolvers.h
 *
 *  Subroutines for the resolution of contact problems with friction:\n
 *
 *  Try \f$(reaction,velocity)\f$ such that:\n
 *
 *  \f$
 *   \left\lbrace
 *    \begin{array}{l}
 *     M reaction + q = velocity \\
 *     0 \le reaction_n \perp velocity_n \ge 0\\
 *     -velocity_t \in \partial\psi_{[-\mu reaction_n, \mu reaction_n]}(reaction_t)\\
 *     \end{array}
 *    \right.
 *   \f$
 *
 *  here M is an (n x n) matrix, q, reaction and velocity n-vectors.\n
 *
 *  This system of equations and inequalities is solved thanks to FrictionContact solvers.
 *  The present routine is an interface to connect to dedicated solvers, according to the chosen method.
 *
 *  frictionContact3DSolvers is a generic interface allowing the call of one of the FrictionContact3D solvers listed below:
 *
 *  - frictionContact3D_nsgs(...): non-smooth Gauss-Seidel solver
 *
 *  \author Mathieu Renouf, Franck Perignon.
 *
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

  /* General interface to solvers for friction-contact 3D problem
     \param number of contacts, nc (dim of the problem n = 3*nc)
     \param global matrix M (n*n)
     \param global vector q (n)
     \param method
     \param global vector reaction (n), in-out parameter
     \param global vector velocity (n), in-out parameter
     \param vector of the friction coefficients (size nc)
     \return     result (0 if successful otherwise 1).
  */
  int frictionContact3DSolvers(int, double*, double*, method*, double*, double*, double*);

  /* Non-Smooth Gauss Seidel solver for friction-contact 3D problem
     \param number of contacts, nc (dim of the problem n = 3*nc)
     \param global matrix M (n*n)
     \param global vector q (n)
     \param global vector reaction (n), in-out parameter
     \param global vector velocity (n), in-out parameter
     \param vector of the friction coefficients (size nc)
     \param information about solver result
     \param int vector of parameters (max. iteration number ...)
     \param double vector of parameters (tolerance ...)
  */
  void frictionContact3D_nsgs(int, double*, double*, double*, double*, double*, int *, int *, double*);

  /* Non-Smooth Gauss Seidel solver for friction-contact 3D problem
     \param number of contacts, nc (dim of the problem n = 3*nc)
     \param global matrix M (n*n), save as a SparseBlockStructuredMatrix
     \param global vector q (n)
     \param global vector reaction (n), in-out parameter
     \param global vector velocity (n), in-out parameter
     \param vector of the friction coefficients (size nc)
     \param information about solver result
     \param int vector of parameters (max. iteration number ...)
     \param double vector of parameters (tolerance ...)
  */
  void frictionContact3D_nsgs_SB(int, SparseBlockStructuredMatrix*, double*, double*, double*, double*, int *, int *, double*);

  /** Check for trivial solution in the friction-contact 3D problem
     \param dim of the problem
     \param global vector q (n)
     \param global vector velocity (n), in-out parameter
     \param global vector reaction (n), in-out parameter
     \param int vector of parameters (max. iteration number ...)
     \param double vector of parameters (tolerance ...)
     \return int=0 if a trivial solution has been found, else = -1
  */
  int checkTrivialCase(int, double*, double*, double*, int*, double*);

#ifdef __cplusplus
}
#endif

#endif
