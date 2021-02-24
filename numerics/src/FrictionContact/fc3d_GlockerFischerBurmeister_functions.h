/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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
#ifndef FRICTIONCONTACT3D_GlockerFischerBurmeister_functions_H
#define FRICTIONCONTACT3D_GlockerFischerBurmeister_functions_H

/*!\file fc3d_GlockerFischerBurmeister_functions.h

  \brief Typedef and functions declarations related GlockerFischerBurmeister
 */
#include "SiconosConfig.h" // for BUILD_AS_CPP // IWYU pragma: keep

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif
typedef void (*UpdateSolverPtr)(int, double*);


  void F_GlockerFischerBurmeister(int sizeF, double* reaction, double* FVector, int up2Date);


  void jacobianF_GlockerFischerBurmeister(int sizeF, double* reaction, double* jacobianFMatrix, int up2Date);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
