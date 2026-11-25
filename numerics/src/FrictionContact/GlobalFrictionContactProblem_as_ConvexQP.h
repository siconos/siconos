/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef GLOBALFRICTIONCONTACTPROBLEM_AS_CONVEXQP_H
#define GLOBALFRICTIONCONTACTPROBLEM_AS_CONVEXQP_H

/*!\file GlobalFrictionContactProblem_as_ConvexQP.h
  \brief Definition of a structure to handle with friction-contact (2D or 3D) problems.
*/

#include "NumericsFwd.h"  // for ConvexQP, GlobalFrictionContactProblem, Sol...

/** \struct GlobalFrictionContactProblem_as_ConvexQP GlobalFrictionContactProblem_as_ConvexQP.h
 *
 */
struct GlobalFrictionContactProblem_as_ConvexQP {
  /* the ConvexQP associated with the FC3D problem */
  ConvexQP* cqp;
  /* the GFC3D associated with the ConvexQP  */
  GlobalFrictionContactProblem* gfc3d;
  /* the SolverOptions that might be used to pass some numerical parameters */
  SolverOptions* options;
};

#if defined(__cplusplus)
extern "C" {
#endif

void Projection_ConvexQP_GFC3D_DualCone(void* cqpIn, double* x, double* PX);

#if defined(__cplusplus)
}
#endif

#endif
