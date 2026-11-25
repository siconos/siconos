// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "LinearComplementarityProblem.h"       // for LinearComplement...
#include "MixedLinearComplementarityProblem.h"  // for MixedLinearComplement...

#ifndef MLCP_TO_LCP_H
#define MLCP_TO_LCP_H

#if defined(__cplusplus)
extern "C" {
#endif

LinearComplementarityProblem* mlcp_to_lcp(MixedLinearComplementarityProblem* problem);

#if defined(__cplusplus)
}
#endif

#endif
