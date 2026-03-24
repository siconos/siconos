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
#ifndef ROLLING_FRICTION_naturalmap_functions_H
#define ROLLING_FRICTION_naturalmap_functions_H

/*!\file rolling_naturalmap_functions.h

  Typedef and functions declarations related to naturalmap map
  formulation for rolling friction contact problems.


 */
#include "NumericsFwd.h"    // for Plasticity2DProblem
#include "SiconosConfig.h"  // for BUILD_AS_CPP // IWYU pragma: keep

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C" {
#endif

void rolling_friction_3D_computeNaturalMap(double reaction[5], double velocity[5], double mu,
                                           double mur, double rho[1], double result[5],
                                           double A[25], double B[25]);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
