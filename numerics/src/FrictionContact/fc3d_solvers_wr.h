/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef FRICTIONCONTACT3DSOLVERS_WR_H
#define FRICTIONCONTACT3DSOLVERS_WR_H

/*!\file gfc3d_Solvers.h
  Subroutines for the resolution of contact problems with friction (3-dimensional case).

*/
#include "NumericsFwd.h"  // for FrictionContactProblem, GlobalFrictionConta...

#if defined(__cplusplus)
extern "C" {
#endif

/** Expand a local FC3D problem into a GlobalFrictionContactProblem. */
GlobalFrictionContactProblem* fc3d_expand_to_global_problem(FrictionContactProblem* problem);

/* Backward compatibility function name (deprecated) */
GlobalFrictionContactProblem* fc3d_reformulation_global_problem(
    FrictionContactProblem* problem);

#if defined(__cplusplus)
}
#endif

#endif
