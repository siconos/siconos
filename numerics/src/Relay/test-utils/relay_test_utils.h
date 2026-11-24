/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef RELAY_TEST_UTILS_H
#define RELAY_TEST_UTILS_H

#include "test_utils.h"  // for TestCase

#if defined(__cplusplus)
extern "C" {
#endif
void relay_fillParamWithRespectToSolver(SolverOptions* options, char* solvername,
                                        RelayProblem* problem);
int relay_test_function(TestCase*);
#if defined(__cplusplus)
}
#endif

#endif
