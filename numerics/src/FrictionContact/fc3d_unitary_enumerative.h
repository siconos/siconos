/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#ifndef FRICTIONCONTACT3DUNITARY_ENUMERATIVE_H
#define FRICTIONCONTACT3DUNITARY_ENUMERATIVE_H

/*!\file fc3d_unitary_enumerative.h
  \brief Typedef and functions declarations related to the quartic solver for 3 dimension frictional contact problems.

  Each solver must have 4 functions in its interface:
  - initialize: link local static variables to the global ones (M,q,...)
  - update: link/fill the local variables corresponding to sub-blocks of the full problem, for a specific contact
  - solve: solve the local problem
  - free

*/
#include "SparseBlockMatrix.h"
#include "SolverOptions.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif
  void fc3d_unitary_enumerative_free(FrictionContactProblem* problem);
  void fc3d_unitary_enumerative_initialize(FrictionContactProblem* problem);
  /*API for the nsgs*/
  int fc3d_unitary_enumerative_solve(FrictionContactProblem* problem, double * reaction, SolverOptions* options);
  int fc3d_unitary_enumerative_solve_sliding(FrictionContactProblem* problem, double * reaction, SolverOptions* options);
  int fc3d_unitary_enumerative_test_non_sliding(FrictionContactProblem* problem, double * reaction, double * velocity, SolverOptions* options);
  int fc3d_unitary_enumerative(FrictionContactProblem* problem, double * reaction, double * velocity, int *info, SolverOptions* options);
  int fc3d_unitary_enumerative_setDefaultSolverOptions(SolverOptions* options);
  int fc3d_unitary_enumerative_solve_poly_nu_sliding(FrictionContactProblem* problem, double * reaction, SolverOptions* options);
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif
