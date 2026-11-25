// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef MLCP_SIMPLEX_H
#define MLCP_SIMPLEX_H

/*
 *
 */

#include "NumericsFwd.h"  // for MixedLinearComplementarityProblem, SolverOp...
void mlcp_simplex_init(MixedLinearComplementarityProblem* problem, SolverOptions* options);
void mlcp_simplex_reset(void);

#endif  // MLCP_SIMPLEX_H
