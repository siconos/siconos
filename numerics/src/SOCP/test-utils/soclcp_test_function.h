// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef SOCLCP_TEST_FUNCTION_H
#define SOCLCP_TEST_FUNCTION_H

#include <stdio.h>  // for FILE

#include "NumericsFwd.h"  // for SolverOptions

#if defined(__cplusplus)
extern "C" {
#endif

int soclcp_test_function(FILE* f, SolverOptions* options);

#if defined(__cplusplus)
}
#endif

#endif
