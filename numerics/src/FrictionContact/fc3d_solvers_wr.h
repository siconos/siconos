/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#ifndef FRICTIONCONTACT3DSOLVERS_WR_H
#define FRICTIONCONTACT3DSOLVERS_WR_H

/*!\file gfc3d_Solvers.h
  Subroutines for the resolution of contact problems with friction (3-dimensional case).

*/
#include "NumericsFwd.h"  // for FrictionContactProblem, GlobalFrictionConta...
#include "SiconosConfig.h" // for BUILD_AS_CPP // IWYU pragma: keep

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  GlobalFrictionContactProblem * fc3d_reformulation_global_problem(FrictionContactProblem* problem);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
