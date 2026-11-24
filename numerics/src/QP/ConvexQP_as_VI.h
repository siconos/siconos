/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef CONVEXQP_AS_VI_H
#define CONVEXQP_AS_VI_H

/*!\file ConvexQP_as_VI.h
  \brief Definition of a structure to handle with ConvexQP problems.
*/

#include "NumericsFwd.h"

/** \struct ConvexQP_as_VI ConvexQP_as_VI.h
 *
 */
struct ConvexQP_as_VI {
  /* the VI associated with the CQP problem */
  VariationalInequality* vi;
  /* the CQP associated with the VI  */
  ConvexQP* cqp;
  /* the SolverOptions that might be used to pass some numerical parameters */
  SolverOptions* options;
};

#if defined(__cplusplus)
extern "C" {
#endif

void Function_VI_CQP(void* self, int n, double* x, double* F);

void Projection_VI_CQP(void* viIn, double* x, double* PX);

#if defined(__cplusplus)
}
#endif

#endif
