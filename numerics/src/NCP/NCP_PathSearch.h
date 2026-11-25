// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef NCP_PATHSEARCH_H
#define NCP_PATHSEARCH_H

/*! \file NCP_PathSearch.h
 * \brief Functions for the pathsearch for NCP
 *
 */

#if defined(__cplusplus)
#undef restrict
#define restrict __restrict
#endif

#if defined(__cplusplus)
extern "C" {
#endif

/** compute x from z for an NCP
 * \param n size of the problem
 * \param z current iterate
 * \param F value of the NCP function
 * \param[out] x current newton iterate
 * */
static inline void ncp_pathsearch_compute_x_from_z(unsigned n, double* restrict z,
                                                   double* restrict F, double* restrict x) {
  /* init value of x */
  /* see Linear Algebra Enhancements to the PATH Solver by Li, Ferris and Munson */
  for (unsigned i = 0; i < n; ++i) {
    /* XXX F[i] > 0.0 or F[i] > DBL_EPSILON ? */
    if ((z[i] <= 0.0) && (F[i] > 0.0))
      x[i] = -F[i];
    else
      x[i] = z[i];
  }
}

#if defined(__cplusplus)
}
#endif

#endif
