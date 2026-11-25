/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
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
#elif defined(HAS_ACCELERATE)
// Accelerate: the interface of macOS >= 13.3 (the former cblas.h is deprecated).
#ifndef ACCELERATE_NEW_LAPACK
#define ACCELERATE_NEW_LAPACK
#endif
#include <cblas_new.h>
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
