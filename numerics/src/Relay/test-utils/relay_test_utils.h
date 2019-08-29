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
#ifndef RELAY_TEST_UTILS_H
#define RELAY_TEST_UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include "relay_cst.h"
#include "SolverOptions.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif
  void relay_fillParamWithRespectToSolver(SolverOptions *options, char * solvername, RelayProblem* problem);
  int relay_test_function(FILE * f, int solverId);

  /** Defines the list of data files
      returns an array of char
  */
  char ** data_collection(void);
  
  /** Defines the set of tests (including parameters)
      returns an 'array' of tests, some kind of dict.
  */
  char *** test_collection(int, char **);

  

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif


