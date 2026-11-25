// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef COHESIVEFRICTIONCONTACT_TEST_UTILS_H
#define COHESIVEFRICTIONCONTACT_TEST_UTILS_H

#include "test_utils.h"  // for TestCase

#if defined(__cplusplus)
extern "C" {
#endif

/** Solves a cohesive friction-contact problem read from current->filename with
    current->options.
    \return 0 if the solver succeeded and the solution is finite.
*/
int cohesiveFrictionContact_test_function(TestCase* current);

/** Creates a test collection (a 'list' of tests, each test being a TestCase, i.e. options +
    input data). Must be implemented by each collection (see e.g. test_cohesive_nsgs_1.c).

    \param n_data number of data files
    \param data_collection 'list' of data files
    \param[out] number_of_tests number of tests
    \return an array of tests, to be executed with run_test_collection
*/
TestCase* build_test_collection(int n_data, const char** data_collection,
                                int* number_of_tests);

#if defined(__cplusplus)
}
#endif

#endif
