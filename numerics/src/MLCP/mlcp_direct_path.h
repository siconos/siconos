// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef MLCP_DIRECT_PATH_H
#define MLCP_DIRECT_PATH_H

/*
 * who use the direct solver:
 * First, you have to call the mlcp_direct_path_init().
 * Use a lot the mlcp_direct_path function
 * mlcp_direct_path_reset at the end.
 *
 *
 */

#include "NumericsFwd.h"  // for MixedLinearComplementarityProblem, SolverOp...
int mlcp_direct_path_getNbIWork(MixedLinearComplementarityProblem* problem,
                                SolverOptions* options);
int mlcp_direct_path_getNbDWork(MixedLinearComplementarityProblem* problem,
                                SolverOptions* options);

void mlcp_direct_path_init(MixedLinearComplementarityProblem* problem, SolverOptions* options);
void mlcp_direct_path_reset(void);

#endif  // MLCP_DIRECT_PATH_H
