/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef GLOBALFRICTIONCONTACTPROBLEM_AS_VI_H
#define GLOBALFRICTIONCONTACTPROBLEM_AS_VI_H

/*!\file GlobalFrictionContactProblem_as_VI.h
  \brief Definition of a structure to handle with global friction-contact (2D or 3D) problems.
*/

#include "NumericsFwd.h"  // for GlobalFrictionContactProblem, VariationalIn...

/** \struct GlobalFrictionContactProblem_as_VI GlobalFrictionContactProblem_as_VI.h
 *
 */
struct GlobalFrictionContactProblem_as_VI {
  /* the VI associated with the FC3D probelem */
  VariationalInequality* vi;
  /* the FC3D associated with the VI  */
  GlobalFrictionContactProblem* gfc3d;
};

#if defined(__cplusplus)
extern "C" {
#endif

void Function_VI_GFC3D(void* self, int n, double* x, double* F);

void Projection_VI_GFC3D(void* viIn, double* x, double* PX);

void Function_VI_GFC2D(void* self, int n, double* x, double* F);

void Projection_VI_GFC2D(void* viIn, double* x, double* PX);

#if defined(__cplusplus)
}
#endif

#endif
