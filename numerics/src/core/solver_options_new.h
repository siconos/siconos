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

/*!\file solver_options_new.h
 * \brief Enhanced solver options creation using registration system
 *
 * This header provides an improved solver_options_create function that
 * uses the solver registration system to:
 * - Get default parameters (max_iter, tolerance) from registry
 * - Call solver-specific initialization if provided
 * - Support creating options by name as well as by ID
 */

#ifndef SOLVER_OPTIONS_NEW_H
#define SOLVER_OPTIONS_NEW_H

#include "SolverOptions.h"
#include "utils/solver_registry.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Create solver options using registration system.
 * 
 * This is the NEW implementation that uses the solver registry to:
 * 1. Look up solver metadata (name, defaults)
 * 2. Set default parameters from registration
 * 3. Call solver-specific initialization if available
 * 
 * @param solver_id The solver ID (e.g., FC3D_NSGS)
 * @return Pointer to created options, or NULL if solver not found
 */
SolverOptions* solver_options_create_registered(solver_id_t solver_id);

/**
 * Create solver options by name.
 * 
 * Convenience function to create options using solver name instead of ID.
 * Useful for creating options from user input (command line, config files).
 * 
 * @param solver_name The solver name (e.g., "FC3D_NSGS")
 * @return Pointer to created options, or NULL if name not found
 */
SolverOptions* solver_options_create_by_name(const char* solver_name);

/**
 * Create options and apply solver initialization.
 * 
 * This is a convenience wrapper that creates options AND calls the solver's
 * init function (if provided). Use this when you want fully initialized options.
 * 
 * @param solver_id The solver ID
 * @param problem Optional problem pointer (for problem-specific init)
 * @return Pointer to created and initialized options, or NULL on failure
 */
SolverOptions* solver_options_create_and_init(solver_id_t solver_id, void* problem);

/**
 * Print solver options with defaults from registry.
 * 
 * Shows the current options AND the registered defaults for comparison.
 * Useful for debugging and user documentation.
 * 
 * @param options The options to print
 */
void solver_options_print_with_defaults(const SolverOptions* options);

/**
 * Reset options to registered defaults.
 * 
 * Restores all parameters to their registered default values.
 * Useful for "reset to defaults" functionality.
 * 
 * @param options The options to reset (modified in place)
 */
void solver_options_reset_to_defaults(SolverOptions* options);

/**
 * List all available solvers with their default parameters.
 * 
 * Prints a formatted table of all registered solvers and their defaults.
 * Useful for help messages and documentation.
 */
void solver_options_list_all_defaults(void);

#ifdef __cplusplus
}
#endif

#endif /* SOLVER_OPTIONS_NEW_H */
