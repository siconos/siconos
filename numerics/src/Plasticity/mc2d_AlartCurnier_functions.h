/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
#ifndef MOHRCOULOMB2D_AlartCurnier_functions_H
#define MOHRCOULOMB2D_AlartCurnier_functions_H

/*!\file fc3d_AlartCurnier_functions.h

  Typedef and functions declarations related to Alart-Curnier
  formulation for Mohr Coulomb 2D contact problems.

  Subroutines used when the Mohr Coulomb problem is written
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
#include "NumericsFwd.h"    // for MohrCoulomb2DProblem
#include "SiconosConfig.h"  // for BUILD_AS_CPP // IWYU pragma: keep

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C" {
#endif

void mc2d_compute_rho_split_spectral_norm_cond(MohrCoulomb2DProblem* localproblem,
                                               double* rho);

void mc2d_compute_rho_split_spectral_norm(MohrCoulomb2DProblem* localproblem, double* rho);

void mc2d_compute_rho_spectral_norm(MohrCoulomb2DProblem* localproblem, double* rho);

void mc2d_computeAlartCurnierSTD(double reaction[3], double velocity[3], double mu, double eta,
                                 double rho[3], double result[3], double A[9], double B[9]);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
