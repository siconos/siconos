/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#ifndef ROLLINGFRICTIONCONTACT2D_local_problem_tools_H
#define ROLLINGFRICTIONCONTACT2D_local_problem_tools_H

/*!\file

 */
#include "NumericsFwd.h"  // for RollingFrictionContactProblem

#if defined(__cplusplus)
extern "C" {
#endif
RollingFrictionContactProblem* rolling_friction_2d_local_problem_allocate(
    RollingFrictionContactProblem* problem);
void rolling_friction_2d_local_problem_free(RollingFrictionContactProblem* localproblem,
                                            RollingFrictionContactProblem* problem);
void rolling_friction_2d_local_problem_compute_q(RollingFrictionContactProblem* problem,
                                                 RollingFrictionContactProblem* localproblem,
                                                 double* reaction, int contact);
void rolling_friction_2d_local_problem_fill_M(RollingFrictionContactProblem* problem,
                                              RollingFrictionContactProblem* localproblem,
                                              int contact);

#if defined(__cplusplus)
}
#endif

#endif
