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

/*!\file fclib_interface.h
 * \brief interface to fclib
 */

#ifndef FCLIB_INTERFACE_H
#define FCLIB_INTERFACE_H

#include "SiconosConfig.h"  // for BUILD_AS_CPP, WITH_FCLIB

#if defined(WITH_FCLIB)
#include "NumericsFwd.h"    // for FrictionContactProblem, GlobalFrictionCon...

typedef struct fclib_local fclib_local;
typedef struct fclib_global fclib_global;
typedef struct fclib_global_rolling fclib_global_rolling;
typedef struct fclib_solution fclib_solution;

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  FrictionContactProblem* from_fclib_local(const fclib_local *fclib_problem);


  FrictionContactProblem* frictionContact_fclib_read(const char *path);

  int frictionContact_fclib_write(FrictionContactProblem* problem,
                                  char * title, char * description,
                                  char * mathInfo,
                                  const char *path, int ndof);

  int frictionContact_fclib_write_guess( double * reaction, double * velocity,
                                         const char *path);

  GlobalFrictionContactProblem* from_fclib_global(const fclib_global *fclib_problem);


  GlobalFrictionContactProblem* globalFrictionContact_fclib_read(const char *path);


  int globalFrictionContact_fclib_write(GlobalFrictionContactProblem* problem,
                                        char * title, char * description,
                                        char * mathInfo,
                                        const char *path);

  GlobalRollingFrictionContactProblem* from_fclib_global_rolling(const fclib_global_rolling *fclib_problem);

  GlobalRollingFrictionContactProblem* globalRollingFrictionContact_fclib_read(const char *path);

  int globalRollingFrictionContact_fclib_write(GlobalRollingFrictionContactProblem* problem,
                                               char * title,
                                               char * description,
                                               char * mathInfo,
                                               const char *path);
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif
#endif


#endif
