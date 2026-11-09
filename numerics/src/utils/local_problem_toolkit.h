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

/*!\file local_problem_toolkit.h
 * \brief Shared helpers for NSGS local-problem toolkit allocation.
 *
 * Several problem families (FrictionContact, CohesiveFrictionContact,
 * Plasticity) use a small toolkit structure that stores function pointers for
 * the local solver, update, post-processing, relaxation and error helpers.
 * The allocation pattern is identical in every module: malloc the structure
 * and zero-initialize all fields.  This header provides a single inline helper
 * so the boilerplate is not duplicated.
 */

#ifndef LOCAL_PROBLEM_TOOLKIT_H
#define LOCAL_PROBLEM_TOOLKIT_H

#include <stdlib.h>  // for malloc
#include <string.h>  // for memset

/** Allocate and zero-initialize a local-problem toolkit structure.
 *
 * \param size size in bytes of the toolkit structure (use sizeof)
 * \return a pointer to the zero-initialized structure, or NULL on allocation failure
 */
static inline void* local_problem_toolkit_new(size_t size) {
  void* toolkit = malloc(size);
  if (toolkit) memset(toolkit, 0, size);
  return toolkit;
}

#endif /* LOCAL_PROBLEM_TOOLKIT_H */
