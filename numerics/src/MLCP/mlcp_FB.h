// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef MLCP_FB_H
#define MLCP_FB_H

#include "NumericsFwd.h"  // for MixedLinearComplementarityProblem, SolverOp...

void mlcp_FB_init(MixedLinearComplementarityProblem* problem, SolverOptions* options);
void mlcp_FB_reset(void);

int mlcp_FB_getNbIWork(MixedLinearComplementarityProblem* problem, SolverOptions* options);
int mlcp_FB_getNbDWork(MixedLinearComplementarityProblem* problem, SolverOptions* options);

#endif  // MLCP_FB_H
