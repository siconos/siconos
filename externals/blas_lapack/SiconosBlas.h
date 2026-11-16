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

#ifndef SiconosBlas_H
#define SiconosBlas_H

#if defined(__cplusplus)
extern "C" {
#endif

// tells include-what-you-use to keep this file
// and not to suggest cblas.h or alike.
// IWYU pragma: begin_exports
#if defined(HAS_MKL_CBLAS)
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif
// IWYU pragma: end_exports

#ifdef __cplusplus
}
#undef restrict
#define restrict __restrict
#endif

// blasint: integer type of the cblas interface.
// OpenBLAS defines it (openblas_config.h, included by its cblas.h).
// Define it for the other implementations (Accelerate, MKL, reference cblas ...).
#if !defined(OPENBLAS_CONFIG_H)
#if defined(HAS_MKL_CBLAS)
typedef MKL_INT blasint;
#elif defined(CBLAS_INT)
typedef CBLAS_INT blasint;
#elif SIZEOF_BLASINT == 8
typedef long long blasint;
#else
typedef int blasint;
#endif
#endif

// BLASINT_SIZE set by cmake
#if SIZEOF_BLASINT == 8
#define BLASINT_MAX LLONG_MAX
#else
#define BLASINT_MAX INT_MAX
#endif

#endif  // SiconosBlas_H
