// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef MLCP_DIRECT_ENUM_H
#define MLCP_DIRECT_ENUM_H

/*
 * who use the direct solver:
 * First, you have to call the mlcp_direct_enum_init().
 * Use a lot the mlcp_direct_enum function
 * mlcp_direct_enum_reset at the end.
 *
 *
 */

#include "NumericsFwd.h"  // for MixedLinearComplementarityProblem, SolverOp...
int mlcp_direct_enum_getNbIWork(MixedLinearComplementarityProblem* problem,
                                SolverOptions* options);
int mlcp_direct_enum_getNbDWork(MixedLinearComplementarityProblem* problem,
                                SolverOptions* options);

void mlcp_direct_enum_init(MixedLinearComplementarityProblem* problem, SolverOptions* options);
void mlcp_direct_enum_reset(void);

#endif  // MLCP_DIRECT_ENUM_H
