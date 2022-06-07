/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#ifndef NCP_PATHSEARCH_H
#define NCP_PATHSEARCH_H

/*! \file NCP_PathSearch.h
 * \brief Functions for the pathsearch for NCP
 *
 */

#include "SiconosConfig.h" // for BUILD_AS_CPP // IWYU pragma: keep

#if defined(__cplusplus)
#undef restrict
#define restrict __restrict
#endif

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

  /** compute x from z for an NCP
   * \param n size of the problem
   * \param z current iterate
   * \param F value of the NCP function
   * \param[out] x current newton iterate
   * */
  static inline void ncp_pathsearch_compute_x_from_z(unsigned n, double* restrict z, double* restrict F,double* restrict x)
  {
    /* init value of x */
    /* see Linear Algebra Enhancements to the PATH Solver by Li, Ferris and Munson */
    for (unsigned i = 0; i < n; ++i)
    {
      /* XXX F[i] > 0.0 or F[i] > DBL_EPSILON ? */
      if ((z[i] <= 0.0) && (F[i] > 0.0))
        x[i] = -F[i];
      else
        x[i] = z[i];
    }
  }

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif


#endif
