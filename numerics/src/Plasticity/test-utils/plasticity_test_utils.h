/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef PLASTICITY_TEST_UTILS_H
#define PLASTICITY_TEST_UTILS_H

#include "test_utils.h"  // for TestCase

#if defined(__cplusplus)
extern "C" {
#endif

/** Solves plasticity_2d using parameters and reference from a pre-defined TestCase
    return 1 if the test has succeeded.
*/
int plasticity_test_function(TestCase*);

/** Creates a test collection (a 'list' of tests, each test being a TestCase, i.e. options +
   input data).

    this function must be implemented for each tests collection (see e.g. test_nsgs_1.c)

    \param n_data number of ref files
    \param data_collection 'list' of ref files
    \param[out] number of tests
    \return an array of tests, to be executed with run_test_collection
*/
TestCase* build_test_collection(int n_data, const char** data_collection, int*);

#if defined(__cplusplus)
}
#endif

#endif
