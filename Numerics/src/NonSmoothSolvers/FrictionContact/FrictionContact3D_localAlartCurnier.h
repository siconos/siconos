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
#ifndef FRICTIONCONTACT3D_localAlartCurnier_H
#define FRICTIONCONTACT3D_localAlartCurnier_H

/*!\file FrictionContact3D_localAlartCurnier.h

  \brief Typedef and functions declarations related to Alart-Curnier
  formulation for 3 dimension frictional contact problems in Local Coordinates.

  Subroutines used when the friction-contact 3D problem is written
  using Alart-Curnier formulation:

  \f{equation*}
  F(reaction)= \begin{pmatrix}
  velocity - M.reaction - q  \\
  1/rn*[velocity_N - (velocity_N - rn*reaction_N)^+]
  1/rt*[velocity_T - proj(velocity_T - rt*reaction_T)]
  \end{pmatrix}
  \f}

  where M is an n by n matrix, q an n-dimensional vector, reaction an
  n-dimensional vector and velocity an n-dimensional vector.

  We consider a "global" (ie for several contacts) problem, used to
  initialize the static global variables.

  Two different storages are available for M: dense and sparse block.

 */
#include "NumericsConfig.h"
#include "FrictionContactProblem.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /* Nonsmooth Newton solver based on the Alart--Curnier function for the
   * local (reduced) frictional contact problem in the dense form
   * \param problem the problem to solve in dense form
   * \param reaction solution and initial guess for reaction
   * \param velocity solution and initial guess for velocity
   * \param info returned info
   * \param options  the solver options
   */
  void frictionContact3D_localAlartCurnier(
    FrictionContactProblem* problem,
    double *reaction,
    double *velocity,
    int *info,
    SolverOptions *options);

  /* Nonsmooth Newton solver based on the Alart--Curnier function for the
   * local (reduced) frictional contact problem in the sparse block form
   * \param problem the problem to solve in sparse block form
   * \param reaction solution and initial guess for reaction
   * \param velocity solution and initial guess for velocity
   * \param info returned info
   * \param options  the solver options
   */
  void frictionContact3D_sparseLocalAlartCurnier(
    FrictionContactProblem* problem,
    double *reaction,
    double *velocity,
    int *info,
    SolverOptions *options);

  /* Init the nonsmooth Newton solver based on the Alart--Curnier function
   * for the local (reduced) frictional contact problem in the sparse block form
   * Mainly, initialize the sparse linear solver
   * \param options  the solver options
   */
  void frictionContact3D_sparseLocalAlartCurnierInit(
    SolverOptions *options);

  void frictionContact3D_AlartCurnierFunction(
    unsigned int problemSize,
    double *reaction3D,
    double *velocity3D,
    double *mu,
    double *rho3D,
    double *output_blocklist3,
    double *output_blocklist3x3_1,
    double *output_blocklist3x3_2);

  void computeAWpB(
    unsigned int problemSize,
    double *blocklist3x3_1,
    double *blockarray3x3,
    double *blocklist3x3_2,
    double *output_blockarray3x3);

  /* Set the default solver options for the LOCALAC Solver
   * \param options  the solver options
   */
  int frictionContact3D_AlartCurnier_setDefaultSolverOptions(
    SolverOptions* options);




#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
