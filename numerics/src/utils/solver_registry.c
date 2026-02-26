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

/*!\file solver_registry.c
 * \brief Implementation of the solver registration system
 */

#include "solver_registry.h"
#include <stdio.h>
#include <string.h>
#include "utils/numerics_errors.h"

/* Static registry table */
static const SolverEntry* registry[SOLVER_REGISTRY_MAX];
static size_t registry_count = 0;

int solver_registry_register(const SolverEntry* entry) {
    CHECK_NULL(entry);
    if (registry_count >= SOLVER_REGISTRY_MAX) return -1;
    
    /* Check for duplicate ID */
    for (size_t i = 0; i < registry_count; i++) {
        if (registry[i]->id == entry->id) {
            fprintf(stderr, "solver_registry: duplicate solver ID %d\n", entry->id);
            return -1;
        }
    }
    
    registry[registry_count++] = entry;
    return 0;
}

const SolverEntry* solver_registry_lookup(solver_id_t id) {
    for (size_t i = 0; i < registry_count; i++) {
        if (registry[i]->id == id) {
            return registry[i];
        }
    }
    return NULL;
}

const SolverEntry* solver_registry_lookup_by_name(const char* name) {
    if (!name) return NULL;
    
    for (size_t i = 0; i < registry_count; i++) {
        if (registry[i]->name && strcmp(registry[i]->name, name) == 0) {
            return registry[i];
        }
    }
    return NULL;
}

size_t solver_registry_count(void) {
    return registry_count;
}

const SolverEntry** solver_registry_get_all(size_t* count) {
    if (count) *count = registry_count;
    return registry;
}

int solver_registry_exists(solver_id_t id) {
    return solver_registry_lookup(id) != NULL;
}

void solver_registry_print(void) {
    printf("\nRegistered solvers (%zu):\n", registry_count);
    printf("%-10s %-30s %-50s %s\n", "ID", "Name", "Description", "Local");
    printf("%-10s %-30s %-50s %s\n", "--", "----", "-----------", "-----");
    
    for (size_t i = 0; i < registry_count; i++) {
        const SolverEntry* e = registry[i];
        printf("%-10d %-30s %-50s %s\n",
               e->id,
               e->name ? e->name : "N/A",
               e->description ? e->description : "",
               e->is_local_solver ? "yes" : "no");
    }
}
