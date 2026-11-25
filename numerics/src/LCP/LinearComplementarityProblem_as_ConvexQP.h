// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef LCP_AS_CONVEXQP_H
#define LCP_AS_CONVEXQP_H

/*!\file LinearComplementarityProblem_as_ConvexQP.h
  \brief Definition of a structure to handle with LCP problems.
*/

#include "NumericsFwd.h"  // for ConvexQP, LinearComplementarityProblem, Sol...

/** \struct LinearComplementarityProblem_as_ConvexQP LinearComplementarityProblem_as_ConvexQP.h
 *
 */
struct LinearComplementarityProblem_as_ConvexQP {
  /* the ConvexQP associated with the FC3D problem */
  ConvexQP* cqp;
  /* the FC3D associated with the ConvexQP  */
  LinearComplementarityProblem* lcp;
  /* the SolverOptions that might be used to pass some numerical parameters */
  SolverOptions* options;
};

#if defined(__cplusplus)
extern "C" {
#endif

void Projection_ConvexQP_LCP(void* cqpIn, double* x, double* PX);

#if defined(__cplusplus)
}
#endif

#endif
