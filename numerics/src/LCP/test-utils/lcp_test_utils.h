/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef LCP_TEST_UTILS_H
#define LCP_TEST_UTILS_H

#include "test_utils.h"  // for TestCase

#if defined(__cplusplus)
extern "C" {
#endif
/** Solve lcp using parameters and reference from a pre-defined TestCase
    return 1 if the test has succeeded.
 */
int lcp_test_function(TestCase*);

#if defined(__cplusplus)
}
#endif

#endif
