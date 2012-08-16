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
#ifndef FRICTIONCONTACT3D_AlartCurnier_H
#define FRICTIONCONTACT3D_AlartCurnier_H

/*!\file FrictionContact3D_AlartCurnier.h

  \brief Typedef and functions declarations related to Alart-Curnier
  formulation for 3 dimension frictional contact problems.

  Subroutines used when the friction-contact 3D problem is written
  using Alart-Curnier formulation:

  \f{eqnarray*}
  F(reaction)=\left[\begin{array}{c}
  velocity - M.reaction - q  \\
  1/rn*[velocity_N - (velocity_N - rn*reaction_N)^+]
  1/rt*[velocity_T - proj(velocity_T - rt*reaction_T)]
  \end{array}\right]
  \f}

  where M is an n by n matrix, q an n-dimensional vector, reaction an
  n-dimensional vector and velocity an n-dimensional vector.

  We consider a "global" (ie for several contacts) problem, used to
  initialize the static global variables.  Then a "local" (ie for one
  contact => size = 3) problem is built (update function) and solved
  (solve function).

  Two different storages are available for M: dense and sparse block.

  \author Houari Khenous, Franck Perignon

 */
#include "SparseBlockMatrix.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** Initialize friction-contact 3D Alart-Curnier formulation
      \param problem the global problem to solve
      \param localproblem the local problem to solve
  */
  void frictionContact3D_AC_initialize(FrictionContactProblem* problem, FrictionContactProblem* localproblem);

  /** Update friction-contact 3D problem: formalize local problem for one contact
      \param problem the global problem to solve
      \param localproblem the local problem to solve
      \param number (position in global matrix) of the considered contact
      \param reaction global reaction (only the block corresponding to the
      current contact will be modified
      \param options of the solver

      the rest is used to formalize the local problem)
  */
  void frictionContact3D_AC_update(int number, FrictionContactProblem* problem, FrictionContactProblem* localproblem ,
                                   double * reaction, SolverOptions* options);

  /** Retrieve global reaction vector using local problem solution
      \param contactnumber (position in global matrix) of the considered contact
      \param reaction global reaction
  */
  void frictionContact3D_AC_post(int contactnumber, double * reaction);

  /** Computes F function used in Newton process for Alart-Curnier formulation
      \param size of the local problem
      \param localreaction
      \param[in,out] F vector
      \param up2Date boolean = 1 (true) if the problem is uptodate (ie if F or
      its jacobian have already been computed for the current local
      problem)
   */
  void F_AC(int size, double * localreaction, double * F, int up2Date);

  /** Computes the jacobian of F function used in Newton process for
   * Alart-Curnier formulation
      \param size of the local problem
      \param localreaction
      \param[in,out] jacobianF matrix
      \param  up2Date boolean = 1 (true) if the problem is uptodate (ie if F or
      its jacobian have already been computed for the current local
      problem)
  */
  void jacobianF_AC(int size, double * localreaction, double * jacobianF, int up2Date);

  /** Computes FGlobal function with Alart-Curnier formulation, but
   * for all contacts (ie FGlobal = [F(contact)])
      \param reaction : global reaction
      \param[in,out] FGlobal vector
  */
  void computeFGlobal_AC(double* reaction, double* FGlobal);

  /** free memory for friction contact 3D Alart-Curnier solver */
  void frictionContact3D_AC_free();


  void computeAlartCurnierSTD(double reaction[3], double velocity[3],
                              double mu, double rho[3],
                              double result[3], double A[9], double B[9]);

  void computeAlartCurnierCKPS(double reaction[3], double velocity[3],
                               double mu, double rho[3],
                               double result[3], double A[9], double B[9]);

  int LocalNonsmoothNewtonSolver(FrictionContactProblem* localproblem,
                                 double * R, int *iparam, double *dparam);

  int DampedLocalNonsmoothNewtonSolver(FrictionContactProblem* localproblem,
                                       double * R, int *iparam, double *dparam);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
