// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef FRICTIONCONTACT_TEST_UTILS_H
#define FRICTIONCONTACT_TEST_UTILS_H

#include "test_utils.h"  // for TestCase

#if defined(__cplusplus)
extern "C" {
#endif

#ifdef HAVE_GAMS_C_API
/** Extra setup for options related to GAMS solvers family.
 */
void frictionContact_test_gams_opts(SolverOptions* options);
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

/** Creates a test collection (a 'list' of tests, each test being a TestCase, i.e. options +
   input data).

    this function must be implemented for each tests collection (see e.g. test_nsgs_1.c,
   test_fp_1.c and so on)

    \param n_data number of ref files
    \param data_collection 'list' of ref files
    \param[out] number of tests
    \return an array of tests, to be executed with run_test_collection
*/

int globalRollingFrictionContact_test_function(TestCase*);

/** Solves grfc3d using parameters and reference from a pre-defined TestCase
    return 1 if the test has succeeded.
*/

TestCase* build_test_collection(int n_data, const char** data_collection, int*);

#if defined(__cplusplus)
}
#endif

#endif
