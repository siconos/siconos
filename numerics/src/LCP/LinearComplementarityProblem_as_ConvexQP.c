// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "LinearComplementarityProblem_as_ConvexQP.h"

#include <math.h>  // for fmax

#include "ConvexQP.h"                      // for ConvexQP
#include "LinearComplementarityProblem.h"  // for LinearComplementarityProblem
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "numerics_errors.h"
#include "siconos_debug.h"  // for DEBUG_PRINT

void Projection_ConvexQP_LCP(void* cqpIn, double* x, double* PX) {
  DEBUG_PRINT("Projection_ConvexQP_LCP(void *cqpIn, double *x, double *PX)\n")

  ConvexQP* cqp = (ConvexQP*)cqpIn;
  LinearComplementarityProblem_as_ConvexQP* pb =
      (LinearComplementarityProblem_as_ConvexQP*)cqp->env;
  LinearComplementarityProblem* lcp = pb->lcp;

  int n = lcp->size;
  for (int i = 0; i < n; ++i) {
    PX[i] = fmax(0, x[i]);
  }
}
