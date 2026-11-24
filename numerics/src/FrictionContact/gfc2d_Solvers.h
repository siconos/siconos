// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef GLOBALFRICTIONCONTACT2DSOLVERS_H
#define GLOBALFRICTIONCONTACT2DSOLVERS_H

/*!\file gfc2d_Solvers.h
  Subroutines for the resolution of contact problems with friction (2-dimensional case).
*/

#include "GlobalFrictionContactProblem.h"  // IWYU pragma: keep

#if defined(__cplusplus)
extern "C" {
#endif

/**
    Non-Smooth Gauss Seidel solver with reformulation for friction-contact 2D problem

    \param problem the friction-contact 2D problem to solve
    \param reaction global vector (n), in-out parameter
    \param velocity global vector (n), in-out parameters
    \param globalVelocity global vector (m), in-out parameters
    \param info return 0 if the solution is found
    \param options the solver options :
    iparam[0] : Maximum iteration number
    dparam[0] : tolerance
    dparam[1] : (out) error
  */
void gfc2d_nsgs_wr(GlobalFrictionContactProblem* problem, double* reaction, double* velocity,
                   double* globalVelocity, int* info, SolverOptions* options);

#if defined(__cplusplus)
}
#endif

#endif
