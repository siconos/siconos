/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef FRICTIONCONTACTPROBLEM_AS_VI_H
#define FRICTIONCONTACTPROBLEM_AS_VI_H

/*!\file FrictionContactProblem_as_VI.h
  \brief Definition of a structure to handle with friction-contact (2D or 3D) problems.
*/

#include "NumericsFwd.h"  // for FrictionContactProblem, SolverOptions, Vari...

/** \struct FrictionContactProblem_as_VI FrictionContactProblem_as_VI.h
 *
 */
struct FrictionContactProblem_as_VI {
  /* the VI associated with the FC3D problem */
  VariationalInequality* vi;
  /* the FC3D associated with the VI  */
  FrictionContactProblem* fc3d;
  /* the SolverOptions that might be used to pass some numerical parameters */
  SolverOptions* options;
};

#if defined(__cplusplus)
extern "C" {
#endif

void Function_VI_FC3D(void* self, int n, double* x, double* F);

void Projection_VI_FC3D(void* viIn, double* x, double* PX);

void Function_VI_FC3D_Cylinder(void* self, int n, double* x, double* F);

void Projection_VI_FC3D_Cylinder(void* viIn, double* x, double* PX);

void Function_VI_FC3D_Disk(void* self, int n, double* x, double* F);

void Projection_VI_FC3D_Disk(void* viIn, double* x, double* PX);

#if defined(__cplusplus)
}
#endif

#endif
