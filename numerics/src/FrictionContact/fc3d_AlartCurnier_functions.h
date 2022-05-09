/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#ifndef FRICTIONCONTACT3D_AlartCurnier_functions_H
#define FRICTIONCONTACT3D_AlartCurnier_functions_H

/*!\file fc3d_AlartCurnier_functions.h

  \brief Typedef and functions declarations related to Alart-Curnier
  formulation for 3 dimension frictional contact problems.

  Subroutines used when the friction-contact 3D problem is written
  using Alart-Curnier formulation:

  \rst

  .. math::
     :nowrap:

     F(reaction)=\left[\begin{array}{c}
     velocity - M.reaction - q  \\
     1/rn*[velocity_N - (velocity_N - rn*reaction_N)^+]
     1/rt*[velocity_T - proj(velocity_T - rt*reaction_T)]
     \end{array}\right]

  \endrst

  where M is an n by n matrix, q an n-dimensional vector, reaction an
  n-dimensional vector and velocity an n-dimensional vector.

  We consider a "global" (ie for several contacts) problem, used to
  initialize the static global variables.  Then a "local" (ie for one
  contact => size = 3) problem is built (update function) and solved
  (solve function).

  Two different storages are available for M: dense and sparse block.

 */
#include "NumericsFwd.h"  // for FrictionContactProblem
#include "SiconosConfig.h" // for BUILD_AS_CPP // IWYU pragma: keep

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  void compute_rho_split_spectral_norm_cond(FrictionContactProblem* localproblem, double * rho);

  void compute_rho_split_spectral_norm(FrictionContactProblem* localproblem, double * rho);

  void compute_rho_spectral_norm(FrictionContactProblem* localproblem, double * rho);

  void computeAlartCurnierSTD(double reaction[3], double velocity[3],
                              double mu, double rho[3],
                              double result[3], double A[9], double B[9]);

  void computeAlartCurnierJeanMoreau(double reaction[3], double velocity[3],
                               double mu, double rho[3],
                               double result[3], double A[9], double B[9]);

  /* /\** Computes F function used in Newton process for Alart-Curnier formulation */
  /*     \param size of the local problem */
  /*     \param localreaction */
  /*     \param[in,out] F vector */
  /*     \param up2Date boolean = 1 (true) if the problem is uptodate (ie if F or */
  /*     its jacobian have already been computed for the current local */
  /*     problem) */
  /* *\/ */
  /* void F_AC(int size, double * localreaction, double * F, int up2Date); */

  /* /\** Computes the jacobian of F function used in Newton process for */
  /*  * Alart-Curnier formulation */
  /*     \param size of the local problem */
  /*     \param localreaction */
  /*     \param[in,out] jacobianF matrix */
  /*     \param  up2Date boolean = 1 (true) if the problem is uptodate (ie if F or */
  /*     its jacobian have already been computed for the current local */
  /*     problem) */
  /* *\/ */
  /* void jacobianF_AC(int size, double * localreaction, double * jacobianF, int up2Date); */
  /* /\** Computes FGlobal function with Alart-Curnier formulation, but */
  /*  * for all contacts (ie FGlobal = [F(contact)]) */
  /*     \param reaction : global reaction */
  /*     \param[in,out] FGlobal vector */
  /* *\/ */
  /* void computeFGlobal_AC(double* reaction, double* FGlobal); */

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
