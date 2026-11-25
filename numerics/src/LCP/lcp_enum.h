// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file lcp_enum.h
 * \brief functions for solving LCP using enumeration */

#ifndef LCP_ENUM_H
#define LCP_ENUM_H

#include "NumericsFwd.h"  // for LinearComplementarityProblem, SolverOptions

#if defined(__cplusplus)
extern "C" {
#endif
int lcp_enum_getNbIWork(LinearComplementarityProblem* problem, SolverOptions* options);
int lcp_enum_getNbDWork(LinearComplementarityProblem* problem, SolverOptions* options);

#if defined(__cplusplus)
}
#endif

#endif  // LCP_ENUM_H
