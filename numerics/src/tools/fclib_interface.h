/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.

 * Copyright 2016 INRIA.

 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at

 * http://www.apache.org/licenses/LICENSE-2.0

 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#ifndef FCLIB_INTERFACE_H
#define FCLIB_INTERFACE_H

#include "SiconosConfig.h"

#if defined(WITH_FCLIB)
#include <fclib.h>
#endif

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

#if defined(WITH_FCLIB)
  FrictionContactProblem* from_fclib_local(const struct fclib_local *fclib_problem);


  FrictionContactProblem* frictionContact_fclib_read(const char *path);

  int frictionContact_fclib_write(FrictionContactProblem* problem,
                                  char * title, char * description,
                                  char * mathInfo,
                                  const char *path, int ndof);

  GlobalFrictionContactProblem* from_fclib_global(const struct fclib_global *fclib_problem);


  GlobalFrictionContactProblem* globalFrictionContact_fclib_read(const char *path);


  int globalFrictionContact_fclib_write(GlobalFrictionContactProblem* problem,
                                        char * title, char * description,
                                        char * mathInfo,
                                        const char *path);

#endif

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif


