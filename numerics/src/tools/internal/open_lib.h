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

/** \file open_lib.h
 * \brief function to open library and find functions in then (useful for
 * plugin */

#ifndef OPEN_LIB_H
#define OPEN_LIB_H

#include "SiconosConfig.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** open a library and return an handle (casted as void*)
   * \param lib_name name of the library
   * \param flags additional flags (for dlopen)
   */
  void* open_library(const char* lib_name, const int flags);

  /** get the address of a function in an already opened lib
   * \param plugin handle to the library
   * \param func name of function
   */
  void* get_function_address(void* plugin, const char* func);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
