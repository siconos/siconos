/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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
#ifndef FRICTIONCONTACT3D_local_problem_tools_H
#define FRICTIONCONTACT3D_local_problem_tools_H

/*!\file 

 */
#include "NumericsFwd.h"  // for FrictionContactProblem
#include "SiconosConfig.h" // for BUILD_AS_CPP // IWYU pragma: keep

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif
  FrictionContactProblem* fc3d_local_problem_allocate(FrictionContactProblem* problem);
  void fc3d_local_problem_free(FrictionContactProblem* localproblem,
                               FrictionContactProblem* problem);
  void fc3d_local_problem_compute_q(FrictionContactProblem * problem, FrictionContactProblem * localproblem, double *reaction, int contact);
  void fc3d_local_problem_fill_M(FrictionContactProblem * problem, FrictionContactProblem * localproblem, int contact);
  

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
