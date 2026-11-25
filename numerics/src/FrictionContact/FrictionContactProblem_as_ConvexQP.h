/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef FRICTIONCONTACTPROBLEM_AS_CONVEXQP_H
#define FRICTIONCONTACTPROBLEM_AS_CONVEXQP_H

/*!\file FrictionContactProblem_as_ConvexQP.h
  \brief Definition of a structure to handle with friction-contact (2D or 3D) problems.
*/

#include "NumericsFwd.h"  // for ConvexQP, FrictionContactProblem, SolverOpt...
/** \struct FrictionContactProblem_as_ConvexQP FrictionContactProblem_as_ConvexQP.h
 *
 */
struct FrictionContactProblem_as_ConvexQP {
  /* the ConvexQP associated with the FC3D problem */
  ConvexQP* cqp;
  /* the FC3D associated with the ConvexQP  */
  FrictionContactProblem* fc3d;
  /* the SolverOptions that might be used to pass some numerical parameters */
  SolverOptions* options;
};

#if defined(__cplusplus)
extern "C" {
#endif

void Projection_ConvexQP_FC3D_Cylinder(void* cqpIn, double* x, double* PX);

void Projection_ConvexQP_FC3D_Disk(void* cqpIn, double* x, double* PX);

#if defined(__cplusplus)
}
#endif

#endif
