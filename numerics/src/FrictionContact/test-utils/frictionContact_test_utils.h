/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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
#ifndef FRICTIONCONTACT_TEST_UTILS_H
#define FRICTIONCONTACT_TEST_UTILS_H

#include "Friction_cst.h"
#include "test_utils.h" // for TestCase
#include "SiconosConfig.h"  // for BUILD_AS_CPP // IWYU pragma: keep

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

#ifdef HAVE_GAMS_C_API
  /** Extra setup for options related to GAMS solvers family.
   */
  void frictionContact_test_gams_opts(SolverOptions * options);
#endif

  /** Solves fc3d using parameters and reference from a pre-defined TestCase
      return 1 if the test has succeeded.
  */
  int frictionContact_test_function(TestCase*);

  /** Solves gfc3d using parameters and reference from a pre-defined TestCase
      return 1 if the test has succeeded.
  */
  int globalFrictionContact_test_function(TestCase*);

  /** Solves rfc3d using parameters and reference from a pre-defined TestCase
      return 1 if the test has succeeded.
  */
  int rollingFrictionContact_test_function(TestCase*);

  /** Creates a test collection (a 'list' of tests, each test being a TestCase, i.e. options + input data).

      this function must be implemented for each tests collection (see e.g. test_nsgs_1.c, test_fp_1.c and so on)

      \param n_data number of ref files
      \param data_collection 'list' of ref files
      \param[out] number of tests
      \return an array of tests, to be executed with run_test_collection
  */
  TestCase * build_test_collection(int n_data, const char ** data_collection, int*);


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
