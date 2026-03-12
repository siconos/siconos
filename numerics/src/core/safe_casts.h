/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2026 INRIA.
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

/* @file safe_casts.h
 * @brief Tools to handle conversions safely
 *
 * Particularly important for casts between int, unsigned,
 * blas int and so on
 */

#ifndef SAFE_CASTS_H
#define SAFE_CASTS_H

#include <assert.h>
#include <cs.h>  // For CS_INT
#include <limits.h>
#include <openblas_config.h>
#include <stddef.h>

#include "SiconosBlas.h"  // For BLASINT_MAX

static_assert(sizeof(blasint) == SIZEOF_BLASINT, "Error in blasint size detection");

// The routines below must be used for all conversions between integers.
// Implicit casts must be avoided. See the 10000 warnings in devmode.

/** convert to int format used in blas */
static inline blasint to_blasint(size_t value) {
  assert(value <= (size_t)BLASINT_MAX);
  return (blasint)value;
}

/** convert to int format used in cs.h */
static inline CS_INT to_csint(size_t value) {
  assert(value <= (size_t)CS_INT_MAX);
  return (CS_INT)value;
}
static inline int to_int(size_t value) {
  assert(value <= (size_t)INT_MAX);
  return (int)value;
}
static inline size_t to_size_t(int value) {
  assert(value >= 0);
  return (size_t)value;
}

static inline size_t csint_to_size_t(CS_INT value) {
  assert(value >= 0);
  assert((uint64_t)value <= SIZE_MAX);
  return (size_t)value;
}

static inline int csint_to_int(CS_INT value) {
  assert(value >= INT_MIN);
  assert(value <= INT_MAX);
  return (int)value;
}
#endif
