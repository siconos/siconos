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
#ifndef FRICTIONCONTACT3D_GlobalAlartCurnier_H
#define FRICTIONCONTACT3D_GlobalAlartCurnier_H

/*!\file FrictionContact3D_GlobalAlartCurnier.h

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
  initialize the static global variables.

  Two different storages are available for M: dense and sparse block.

 */
#ifdef __cplusplus
extern "C"
{
#endif

  void frictionContact3D_globalAlartCurnierFunction(
    unsigned int problemSize,
    double *reaction3D,
    double *velocity3D,
    double *mu,
    double *rho3D,
    double *result,
    double *result3x3_1,
    double *result3x3_2);

  void frictionContact3D_globalAlartCurnier(
    FrictionContactProblem* problem,
    double *reaction,
    double *velocity,
    int *info,
    SolverOptions *options);

  int frictionContact3D_globalAlartCurnier_setDefaultSolverOptions(
    SolverOptions* options);

#ifdef WITH_MUMPS
  void frictionContact3D_sparseGlobalAlartCurnierInit(
    SolverOptions *SO);

  void frictionContact3D_sparseGlobalAlartCurnier(
    FrictionContactProblem* problem,
    double *reaction,
    double *velocity,
    int *info,
    SolverOptions *options);
#endif

#ifdef __cplusplus
}
#endif

#endif
