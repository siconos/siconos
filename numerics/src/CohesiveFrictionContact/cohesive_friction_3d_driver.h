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
#ifndef COHESIVE_FRICTION_3D_DRIVER_H
#define COHESIVE_FRICTION_3D_DRIVER_H

/*!\file cohesive_friction_3d_driver.h
  Driver routine for solving 3D cohesive friction-contact problems using
  the solver registration system.

  This driver uses dynamic dispatch via the solver registry, eliminating
  the need for giant switch statements.
*/

#include "CohesiveFrictionContactProblem.h"
#include "SolverOptions.h"
#include "solver_registry.h"

#if defined(__cplusplus)
extern "C" {
#endif

/**
 * Check for trivial case (no contact)
 *
 * \param[in] problem the cohesive friction-contact problem
 * \param[out] velocity the velocity vector
 * \param[out] reaction the reaction vector
 * \param[in] options solver options
 * \return NUMERICS_OK if trivial case, error code otherwise
 */
int cohesive_friction_3d_checkTrivialCase(CohesiveFrictionContactProblem *problem,
                                          double *velocity, double *reaction,
                                          SolverOptions *options);

/**
 * Driver for solving a 3D Cohesive Friction-Contact problem
 *
 * This is the main entry point for solving cohesive friction-contact problems.
 * It uses the solver registration system for dynamic dispatch.
 *
 * \param[in,out] problem the cohesive friction-contact problem to solve
 * \param[out] reaction the reaction vector (solution)
 * \param[out] velocity the velocity vector
 * \param[in,out] options solver options (contains solver ID, tolerances, etc.)
 * \return 0 if successful, error code otherwise
 *
 * \note The reaction and velocity arrays must be pre-allocated with size
 *       dimension * numberOfContacts.
 */
int cohesive_friction_3d_driver(CohesiveFrictionContactProblem *problem,
                                double *reaction,
                                double *velocity,
                                SolverOptions *options);

/**
 * Create solver options for cohesive friction 3D problems
 *
 * This function creates and initializes solver options using the
 * registration system. It validates that the solver exists and is
 * appropriate for cohesive friction 3D problems.
 *
 * \param[in] solver_id the solver identifier to use
 * \return pointer to SolverOptions if successful, NULL otherwise
 *
 * \note The returned options must be freed with solver_options_delete()
 */
SolverOptions* cohesive_friction_3d_solver_options_create(solver_id_t solver_id);

/**
 * List all available cohesive friction 3D solvers
 *
 * Prints a formatted table of all registered solvers suitable for
 * cohesive friction 3D problems.
 */
void cohesive_friction_3d_list_available_solvers(void);

/**
 * Print detailed information about a specific solver
 *
 * \param[in] solver_id the solver identifier
 */
void cohesive_friction_3d_print_solver_info(solver_id_t solver_id);


#if defined(__cplusplus)
}
#endif

#endif  // COHESIVE_FRICTION_3D_DRIVER_H
