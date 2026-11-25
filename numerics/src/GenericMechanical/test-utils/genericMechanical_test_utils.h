// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef GENERICMECHANICAL_TEST_UTILS_H
#define GENERICMECHANICAL_TEST_UTILS_H

#include "test_utils.h"

#if defined(__cplusplus)
extern "C" {
#endif

/** Creates a test collection (a 'list' of tests, each test being a TestCase, i.e. options +
   input data).

    this function must be implemented for each tests collection (see e.g. test_nsgs_1.c,
   test_fp_1.c and so on)

    \param n_data number of ref files
    \param data_collection 'list' of ref files
    \param[out] number of tests
    \return an array of tests, to be executed with run_test_collection
*/
TestCase* build_test_collection(int n_data, const char** data_collection, int*);

/** Solves gmp using parameters and reference from a pre-defined TestCase
    return 1 if the test has succeeded.
*/
int gmp_test_function(TestCase*);

#if defined(__cplusplus)
}
#endif

#endif
