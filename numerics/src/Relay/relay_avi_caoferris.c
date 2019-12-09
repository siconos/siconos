/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
*/

#include "AVI_Solvers.h"
#include "Relay_Solvers.h"
#include "avi_caoferris.h"
#include "relay_cst.h"
#include "LinearComplementarityProblem.h"
#include "NumericsMatrix.h"
//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"


void relay_avi_caoferris(RelayProblem* problem, double *z, double *w, int *info, SolverOptions* options)
{
  unsigned int n = problem->size;
  assert(n > 0);
  unsigned int s = 2*n;

  /* Copy the data from Relay problem */
  LinearComplementarityProblem lcplike_pb;
  lcplike_pb.size = s;
  NumericsMatrix num_mat;
  NM_fill(&num_mat, NM_DENSE, s, s, calloc(s*s, sizeof(double)));

  lcplike_pb.M = &num_mat;

  lcplike_pb.q = (double *)malloc(s*sizeof(double));
  double* b_bar = (double *)malloc(s*sizeof(double));

  /* We can always choose the extreme point such that the matrix of active
   * constrains is the identity and the other matrix is minus the identity.
   * B_A = Id, B_I = -Id, b_A = lb, b_I = -ub
   * TODO: implement the case when the user gives an hint about the solution
   *
   * Siconos/Kernel is giving us column-major matrix and the solver is
   * expecting matrix in this format
   * */

  double tmp;
  for (unsigned int i = 0; i < n; ++i)
  {
    tmp = 0.0;
    lcplike_pb.M->matrix0[i*(s+1)] = 1.0;
    lcplike_pb.M->matrix0[(i + n)*(s+1)] = -1.0;
    lcplike_pb.q[i] = problem->q[i];
    lcplike_pb.q[i+n] = - problem->lb[i] + problem->ub[i];
    for (unsigned j = 0; j < n; ++j)
    {
      lcplike_pb.M->matrix0[i + (j+n)*s] = problem->M->matrix0[i + j*n];
      tmp += problem->M->matrix0[i + j*n]*problem->lb[j];
    }
    /* \bar{a} =  -\tilde{a} + Ay_e */
    lcplike_pb.q[i] += tmp;
  }
  double* d_vec = (double *)malloc(s*sizeof(double));
  for (unsigned i = 0; i<n; ++i)
  {
    d_vec[i] = -1.0;
    d_vec[i+n] = 0;
  }

  /* Set of active constraint is trivial */
  unsigned* A = (unsigned*)malloc(n*sizeof(unsigned));
  for (unsigned i = 0; i < n; ++i) A[i] = i + 1;

  double* u_vec = (double *)calloc(s, sizeof(double));
  double* s_vec = (double *)calloc(s, sizeof(double));
  /* Call directly the 3rd stage 
   * Here w is used as u and z as s in the AVI */
  *info = avi_caoferris_stage3(&lcplike_pb, u_vec, s_vec, d_vec, n, A, options);

  /* Update z  */
  /* XXX why no w ?  */
  DEBUG_PRINT_VEC_INT(A, n);
  for (unsigned i = 0; i<n; ++i) z[i] = s_vec[A[i]-1] + problem->lb[i];
  /* free allocated stuff */
  free(u_vec);
  free(s_vec);
  free(A);
  free(d_vec);
  NM_clear(lcplike_pb.M);
  free(lcplike_pb.q);
  free(b_bar);
}

int relay_avi_caoferris_setDefaultSolverOptions(SolverOptions* options)
{
  solver_options_set(options, SICONOS_RELAY_AVI_CAOFERRIS);
  return 0;
}

