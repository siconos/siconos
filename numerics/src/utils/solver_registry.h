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

/*!\file solver_registry.h
 * \brief Solver registration system for dynamic solver dispatch
 *
 * This is a simplified table-based registration system that works
 * reliably across platforms.
 */

#ifndef SOLVER_REGISTRY_H
#define SOLVER_REGISTRY_H

#include "SolverOptions.h"
#include "numerics_errors.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef int solver_id_t;

typedef int (*SolverInitFunc)(void* problem, SolverOptions* options);
typedef int (*SolverSolveFunc)(void* problem, double* reaction,
                               double* velocity, SolverOptions* options);
typedef void (*SolverFreeFunc)(void* problem, SolverOptions* options);
typedef double (*SolverErrorFunc)(void* problem, double* reaction,
                                  double* velocity, double* work);
typedef void (*SolverSetDefaultFunc)(SolverOptions* options);

typedef struct {
    solver_id_t id;
    const char* name;
    const char* description;
    SolverInitFunc init;
    SolverSolveFunc solve;
    SolverFreeFunc free;
    SolverErrorFunc error;
    SolverSetDefaultFunc set_default;  /* Called during options creation (no problem needed) */
    int default_max_iter;
    double default_tol;
    int is_local_solver;
} SolverEntry;

/* Maximum number of solvers that can be registered */
#define SOLVER_REGISTRY_MAX 128

/* Registration function - called by solvers to register themselves */
int solver_registry_register(const SolverEntry* entry);

/* Main registration macro - all solvers must provide a set_default function */
#define REGISTER_SOLVER(id, name_str, desc_str, init_fn, solve_fn, free_fn, err_fn, \
                         set_default_fn, max_iter_val, tol_val, local_val) \
    static int _solver_reg_##id(void) __attribute__((constructor)); \
    static int _solver_reg_##id(void) { \
        static const SolverEntry entry = { \
            id, name_str, desc_str, init_fn, solve_fn, free_fn, err_fn, set_default_fn, \
            max_iter_val, tol_val, local_val \
        }; \
        return solver_registry_register(&entry); \
    }

/* Registry lookup functions */
const SolverEntry* solver_registry_lookup(solver_id_t id);
const SolverEntry* solver_registry_lookup_by_name(const char* name);
const SolverEntry** solver_registry_get_all(size_t* count);
size_t solver_registry_count(void);
int solver_registry_exists(solver_id_t id);
void solver_registry_print(void);

#ifdef __cplusplus
}
#endif

#endif /* SOLVER_REGISTRY_H */
