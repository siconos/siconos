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

/*!\file lcp_avi_caoferris.c
 \brief Solve an LCP by reformulating it as an AVI and the solver by Cao and
Ferris solves the subsequent AVI.
 \author Olivier Huber
*/

#include "LinearComplementarityProblem.h"
#include "avi_caoferris.h"
#include "AVI_Solvers.h"
#include "LCP_Solvers.h"
#include <assert.h>

void lcp_avi_caoferris(LinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options)
{
  unsigned int n = problem->size;
  assert(n > 0);

  double* d_vec = (double *)malloc(n*sizeof(double));
  for (unsigned i = 0; i < n; ++i) d_vec[i] = -1.0;

  /* Set of active constraint is trivial */
  unsigned* A = (unsigned*)malloc(n*sizeof(unsigned));
  for (unsigned i = 0; i < n; ++i) A[i] = i + 1;

  /* Call directly the 3rd stage */
  *info = avi_caoferris_stage3(problem, w, z, d_vec, n, A, options);

  /* free allocated stuff */
  free(A);
  free(d_vec);
}
