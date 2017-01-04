/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
 * \author Olivier Huber
 */

#include "PathSearch.h"

#include "SiconosBlas.h"
#include "NumericsMatrix.h"
#include "NonlinearComplementarityProblem.h"
#include "LinearComplementarityProblem.h"

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


  /** update the lcp subproblem: M, q and r
   * \param problem the NCP problem to solve
   * \param lcp_subproblem the lcp problem to fill
   * \param n size of the NCP problem
   * \param x_plus positive part of x
   * \param x current newton iterate
   * \param r value of the normal map
   */
  static inline void ncp_pathsearch_update_lcp_data(NonlinearComplementarityProblem* problem, LinearComplementarityProblem* lcp_subproblem, unsigned n, double* restrict x_plus, double* restrict x, double* restrict r)
  {
    /* compute M = nabla F(x_plus) */
    problem->compute_nabla_F(problem->env, n, x_plus, problem->nabla_F);

    /* r = F_+(x) = F(x_+) + x - x_+ */
    /* the real q = q - r = x_+ - x - M x_plus */
    /* q = -M x_plus */
    NM_gemv(-1.0, problem->nabla_F, x_plus, 0.0, lcp_subproblem->q);

    /* first compute r = x - x_+ */
    cblas_dcopy(n, x, 1, r, 1); /* r = x */
    cblas_daxpy(n, -1.0, x_plus, 1, r, 1); /* r -= x_plus */
    /* we factorized computations */
    cblas_daxpy(n, -1.0, r, 1, lcp_subproblem->q, 1); /* q -= x - x_plus */
    /* Finish r */
    problem->compute_F(problem->env, n, x_plus, x); /* compute F(x_plus) */
    cblas_daxpy(n, 1.0, x, 1, r, 1); /* r += F(x_plus) */

  }

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif


#endif
