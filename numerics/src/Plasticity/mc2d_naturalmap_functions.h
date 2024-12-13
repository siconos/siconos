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
#ifndef MOHRCOULOMB2D_naturalmap_functions_H
#define MOHRCOULOMB2D_naturalmap_functions_H

/*!\file fc3d_AlartCurnier_functions.h

  Typedef and functions declarations related to naturalmap map
  formulation for Mohr Coulomb 2D contact problems.


 */
#include "NumericsFwd.h"    // for MohrCoulomb2DProblem
#include "SiconosConfig.h"  // for BUILD_AS_CPP // IWYU pragma: keep

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C" {
#endif

void mc2d_computeNaturalMap(double reaction[3], double velocity[3], double eta, double theta,
                            double rho[3], double result[3], double A[9], double B[9]);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
