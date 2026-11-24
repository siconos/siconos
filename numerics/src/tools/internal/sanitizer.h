// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef SN_SANITIZER_H
#define SN_SANITIZER_H

/*!\file sanitizer.h
 * \brief sanitizer helpers
 */

#include "SiconosBlas.h"

#if defined(__has_feature)
#if __has_feature(memory_sanitizer)
#define MSAN_INIT_VAR(v, n) \
  for (size_t ii = 0; ii < (size_t)n; ++ii) v[ii] = 0.;

#else

#define MSAN_INIT_VAR(v, n)

#endif

#else

#define MSAN_INIT_VAR(v, n)

#endif

/** wrapper around cblas_dcopy which initialize the destination in order to make msan happy
 * \param n length of the data
 * \param src source
 * \param inc_src increment on the source
 * \param dest destination
 * \param inc_src increment on the destination
 */
static inline void cblas_dcopy_msan(int n, const double* src, int inc_src, double* dest,
                                    int inc_dest) {
  MSAN_INIT_VAR(dest, n);
  cblas_dcopy(n, src, inc_src, dest, inc_dest);
}

#endif
