/* Siconos-Numerics, Copyright INRIA 2005-2010.
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
#ifndef FRICTIONCONTACT3D_AlartCurnier_H
#define FRICTIONCONTACT3D_AlartCurnier_H

/*!\file FrictionContact3D_AlartCurnier.h
  \brief Typedef and functions declarations related to Alart-Curnier formulation for 3 dimension frictional contact problems.

  Subroutines used when the friction-contact 3D problem is written using Alart-Curnier formulation:

  \f{eqnarray*}
  F(reaction)=\left[\begin{array}{c}
  velocity - M.reaction - q  \\
  1/rn*[velocity_N - (velocity_N - rn*reaction_N)^+]
  1/rt*[velocity_T - proj(velocity_T - rt*reaction_T)]
  \end{array}\right]
  \f}

  where M is an n by n  matrix, q an n-dimensional vector, reaction an n-dimensional  vector and velocity an n-dimensional vector.\n

  We consider a "global" (ie for several contacts) problem, used to initialize the static global variables.
  Then a "local" (ie for one contact => size = 3) problem is built (update function) and solved (solve function).

  Two different storages are available for M: dense and sparse block.

  \author Houari Khenous, Franck Perignon

 */
#include "SparseBlockMatrix.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /** Initialize friction-contact 3D Alart-Curnier formulation
  */
  void frictionContact3D_AC_initialize(FrictionContactProblem*, FrictionContactProblem*);

  /** Update friction-contact 3D problem: formalize local problem for one contact
      \param number (position in global matrix) of the considered contact
      \param global reaction (only the block corresponding to the current contact will be modified,
      the rest is used to formalize the local problem)
  */
  void frictionContact3D_AC_update(int, FrictionContactProblem* , FrictionContactProblem* , double *, SolverOptions* options);

  /** Retrieve global reaction vector using local problem solution
      \param number (position in global matrix) of the considered contact
      \param global reaction
  */
  void frictionContact3D_AC_post(int, double *);

  /** Computes F function used in Newton process for Alart-Curnier formulation
      \param size of the local problem
      \param local reaction
      \param in-out F vector
      \param bool = 1 (true) if the problem is uptodate (ie if F or its jacobian have already been computed for the current local problem)
   */
  void F_AC(int, double *, double *, int);

  /** Computes the jacobian of F function used in Newton process for Alart-Curnier formulation
      \param size of the local problem
      \param local reaction
      \param in-out jacobianF matrix
      \param bool = 1 (true) if the problem is uptodate (ie if F or its jacobian have already been computed for the current local problem)
  */
  void jacobianF_AC(int, double *, double *, int);

  /** Computes FGlobal function with Alart-Curnier formulation, but for all contacts (ie FGlobal = [F(contact)])
      \param global reaction
      \param in-out FGlobal vector
      \param bool = 1 (true) if the problem is uptodate (ie if F or its jacobian have already been computed for the current local problem)
  */
  void computeFGlobal_AC(double*, double*);

  /** free memory for friction contact 3D Alart-Curnier solver */
  void frictionContact3D_AC_free();


  void computeAlartCurnierSTD(double reaction[3], double velocity[3],
                              double mu, double rho[3],
                              double result[3], double A[9], double B[9]);

  void computeAlartCurnierCKPS(double reaction[3], double velocity[3],
                               double mu, double rho[3],
                               double result[3], double A[9], double B[9]);

#ifdef __cplusplus
}
#endif

#endif
