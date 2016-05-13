/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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

#ifndef Fixe_H
#define Fixe_H

#include "SiconosConfig.h"
#include "NonSmoothNewton.h"

/*!\file NCP_FixedP.h

  \author Houari Khenous, 31/05/2008
*/

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  int Fixe(int n, double* z, int* iparam, double* dparam);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
