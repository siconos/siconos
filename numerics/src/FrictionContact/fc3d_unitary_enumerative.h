// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef FRICTIONCONTACT3DUNITARY_ENUMERATIVE_H
#define FRICTIONCONTACT3DUNITARY_ENUMERATIVE_H

/*!\file fc3d_unitary_enumerative.h
  \brief Typedef and functions declarations related to the quartic solver for 3 dimension
  frictional contact problems.

  Each solver must have 4 functions in its interface:
  - initialize: link local static variables to the global ones (M,q,...)
  - update: link/fill the local variables corresponding to sub-blocks of the full problem, for
  a specific contact
  - solve: solve the local problem
  - free

*/
#include "NumericsFwd.h"  // for FrictionContactProblem, SolverOptions

#if defined(__cplusplus)
extern "C" {
#endif
void fc3d_unitary_enumerative_free(FrictionContactProblem* problem, SolverOptions* dummy2);
void fc3d_unitary_enumerative_initialize(FrictionContactProblem* problem);
/*API for the nsgs*/
int fc3d_unitary_enumerative_solve(FrictionContactProblem* problem, double* reaction,
                                   SolverOptions* options);
int fc3d_unitary_enumerative_solve_sliding(FrictionContactProblem* problem, double* reaction,
                                           SolverOptions* options);
int fc3d_unitary_enumerative_test_non_sliding(FrictionContactProblem* problem,
                                              double* reaction, double* velocity,
                                              SolverOptions* options);
int fc3d_unitary_enumerative(FrictionContactProblem* problem, double* reaction,
                             double* velocity, int* info, SolverOptions* options);
int fc3d_unitary_enumerative_solve_poly_nu_sliding(FrictionContactProblem* problem,
                                                   double* reaction, SolverOptions* options);
#if defined(__cplusplus)
}
#endif

#endif
