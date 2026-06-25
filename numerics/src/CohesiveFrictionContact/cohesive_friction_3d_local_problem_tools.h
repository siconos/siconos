/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
#ifndef COHESIVEFRICTIONCONTACT3D_local_problem_tools_H
#define COHESIVEFRICTIONCONTACT3D_local_problem_tools_H

/*!\file

 */
#include "NumericsFwd.h"  // for FrictionContactProblem
#include "fc3d_local_problem_tools.h"
#if defined(__cplusplus)
extern "C" {
#endif

/** pointer to function used to update local problem */
/* typedef void (*ContactUpdatePtr)(int, CohesiveFrictionContactProblem*,
 * FrictionContactProblem*, double*, */
/*                           SolverOptions*); */
typedef void (*CohesiveUpdatePtr)(int, CohesiveFrictionContactProblem*,
                                  CohesiveFrictionContactProblem*, double*, SolverOptions*);

/** pointer to function used to call local solver */
typedef int (*CohesiveSolverPtr)(CohesiveFrictionContactProblem*, double*, SolverOptions*);

/** pointer to function used to free memory for objects used in nsgs solvers */
typedef void (*CohesiveFreeLocalSolverPtr)(CohesiveFrictionContactProblem*, SolverOptions*);

struct CohesiveLocalProblemFunctionToolkit {
  SolverPtr local_solver_contact;
  FreeLocalSolverPtr free_local_solver_contact;

  CohesiveSolverPtr local_solver_cohesion;
  CohesiveFreeLocalSolverPtr free_local_solver_cohesion;

  CohesiveUpdatePtr update_local_problem;
  /* ContactUpdatePtr update_local_problem_contact; */

  PostSolverPtr post_processed_local_result;

  CopyLocalReactionPtr copy_local_reaction;
  PerformRelaxationPtr perform_relaxation;
  LightErrorSquaredPtr light_error_squared;
  SquaredNormPtr squared_norm;
};

struct CohesiveLocalProblemFunctionToolkit* cohesiveLocalProblemFunctionToolkit_new(void);

void cohesiveLocalProblemFunctionToolkit_display(struct CohesiveLocalProblemFunctionToolkit*);

CohesiveFrictionContactProblem* cohesive_friction_3d_local_problem_allocate(
    NM_types storageType);

void cohesive_friction_3d_local_problem_free(CohesiveFrictionContactProblem* local_problem,
                                             CohesiveFrictionContactProblem* problem);

void cohesive_friction_3d_local_problem_compute_q(CohesiveFrictionContactProblem* problem,
                                                  CohesiveFrictionContactProblem* localproblem,
                                                  double* reaction, int contact);
void cohesive_friction_3d_local_problem_fill_M(CohesiveFrictionContactProblem* problem,
                                               CohesiveFrictionContactProblem* localproblem,
                                               int contact);

#if defined(__cplusplus)
}
#endif

#endif
