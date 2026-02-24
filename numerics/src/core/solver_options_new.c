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

/*!\file solver_options_new.c
 * \brief Implementation of registration-based solver options creation
 */

#include "solver_options_new.h"
#include "numerics_verbose.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Default sizes for iparam and dparam arrays */
#define DEFAULT_IPARAM_SIZE 20
#define DEFAULT_DPARAM_SIZE 10

SolverOptions* solver_options_create_registered(solver_id_t solver_id) {
    /* Look up solver in registry */
    const SolverEntry* solver = solver_registry_lookup(solver_id);
    
    if (!solver) {
        fprintf(stderr, "solver_options_create_registered: solver ID %d not found\n", 
                solver_id);
        return NULL;
    }
    
    /* Allocate options structure */
    SolverOptions* options = (SolverOptions*)calloc(1, sizeof(SolverOptions));
    if (!options) {
        fprintf(stderr, "solver_options_create_registered: memory allocation failed\n");
        return NULL;
    }
    
    /* Set solver ID */
    options->solverId = solver_id;
    
    /* Allocate parameter arrays */
    options->iSize = DEFAULT_IPARAM_SIZE;
    options->dSize = DEFAULT_DPARAM_SIZE;
    options->iparam = (int*)calloc(options->iSize, sizeof(int));
    options->dparam = (double*)calloc(options->dSize, sizeof(double));
    
    if (!options->iparam || !options->dparam) {
        fprintf(stderr, "solver_options_create_registered: parameter allocation failed\n");
        solver_options_delete(options);
        return NULL;
    }
    
    /* Set default parameters from registry */
    options->iparam[SICONOS_IPARAM_MAX_ITER] = solver->default_max_iter;
    options->dparam[SICONOS_DPARAM_TOL] = solver->default_tol;
    
    /* Initialize other common defaults */
    options->iparam[SICONOS_IPARAM_ITER_DONE] = 0;
    options->dparam[SICONOS_DPARAM_RESIDU] = 0.0;
    
    numerics_printf_verbose(2, "Created options for solver '%s': max_iter=%d, tol=%.2e",
                            solver->name, solver->default_max_iter, solver->default_tol);
    
    return options;
}

SolverOptions* solver_options_create_by_name(const char* solver_name) {
    if (!solver_name) {
        fprintf(stderr, "solver_options_create_by_name: NULL solver name\n");
        return NULL;
    }
    
    /* Look up solver by name */
    const SolverEntry* solver = solver_registry_lookup_by_name(solver_name);
    
    if (!solver) {
        fprintf(stderr, "solver_options_create_by_name: solver '%s' not found\n", 
                solver_name);
        fprintf(stderr, "Use solver_options_list_all_defaults() to see available solvers\n");
        return NULL;
    }
    
    return solver_options_create_registered(solver->id);
}

SolverOptions* solver_options_create_and_init(solver_id_t solver_id, void* problem) {
    /* Create options with defaults */
    SolverOptions* options = solver_options_create_registered(solver_id);
    
    if (!options) {
        return NULL;
    }
    
    /* Look up solver for init function */
    const SolverEntry* solver = solver_registry_lookup(solver_id);
    
    if (solver && solver->init) {
        int init_status = solver->init(problem, options);
        
        if (init_status != NUMERICS_OK) {
            fprintf(stderr, "solver_options_create_and_init: init failed for solver '%s'\n",
                    solver->name);
            solver_options_delete(options);
            return NULL;
        }
        
        numerics_printf_verbose(2, "Initialized solver '%s' options", solver->name);
    }
    
    return options;
}

void solver_options_print_with_defaults(const SolverOptions* options) {
    if (!options) {
        printf("solver_options_print_with_defaults: NULL options\n");
        return;
    }
    
    /* Look up solver for default values */
    const SolverEntry* solver = solver_registry_lookup(options->solverId);
    
    printf("\nSolver Options (ID=%d):\n", options->solverId);
    printf("=======================\n");
    
    if (solver) {
        printf("Name:        %s\n", solver->name);
        printf("Description: %s\n", solver->description);
        printf("\nDefault Parameters:\n");
        printf("  max_iter:  %d (registry default: %d)\n",
               options->iparam[SICONOS_IPARAM_MAX_ITER],
               solver->default_max_iter);
        printf("  tolerance: %.2e (registry default: %.2e)\n",
               options->dparam[SICONOS_DPARAM_TOL],
               solver->default_tol);
    } else {
        printf("WARNING: Solver ID %d not found in registry\n", options->solverId);
        printf("Current values:\n");
        printf("  max_iter:  %d\n", options->iparam[SICONOS_IPARAM_MAX_ITER]);
        printf("  tolerance: %.2e\n", options->dparam[SICONOS_DPARAM_TOL]);
    }
    
    printf("\nCurrent State:\n");
    printf("  iterations done: %d\n", options->iparam[SICONOS_IPARAM_ITER_DONE]);
    printf("  residual:        %.2e\n", options->dparam[SICONOS_DPARAM_RESIDU]);
    printf("\n");
}

void solver_options_reset_to_defaults(SolverOptions* options) {
    if (!options) return;
    
    const SolverEntry* solver = solver_registry_lookup(options->solverId);
    
    if (!solver) {
        fprintf(stderr, "solver_options_reset_to_defaults: solver ID %d not found\n",
                options->solverId);
        return;
    }
    
    /* Reset to registered defaults */
    options->iparam[SICONOS_IPARAM_MAX_ITER] = solver->default_max_iter;
    options->dparam[SICONOS_DPARAM_TOL] = solver->default_tol;
    options->iparam[SICONOS_IPARAM_ITER_DONE] = 0;
    options->dparam[SICONOS_DPARAM_RESIDU] = 0.0;
    
    numerics_printf_verbose(1, "Reset options for '%s' to defaults: max_iter=%d, tol=%.2e",
                            solver->name, solver->default_max_iter, solver->default_tol);
}

void solver_options_list_all_defaults(void) {
    printf("\n");
    printf("+------------------------------------------------------------------------------\n");
    printf("| Available Solvers with Default Parameters                                     \n");
    printf("+------------------------------------------------------------------------------\n");
    printf("| %-20s %-8s %-12s %-12s\n", "Name", "ID", "Max Iter", "Tolerance");
    printf("+--------------------+--------+------------+------------\n");
    
    size_t count;
    const SolverEntry** solvers = solver_registry_get_all(&count);
    
    for (size_t i = 0; i < count; i++) {
        const SolverEntry* s = solvers[i];
        printf("| %-20s %-8d %-12d %.2e\n",
               s->name, s->id, s->default_max_iter, s->default_tol);
    }
    
    printf("+--------------------+--------+------------+------------\n");
    printf("| Total: %zu solver(s) registered\n", count);
    printf("+------------------------------------------------------------------------------\n");
    printf("\n");
}
