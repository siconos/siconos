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

#include "AVI_Solvers.h"
#include "Relay_Solvers.h"
#include "avi_caoferris.h"
#include "relay_cst.h"
#include "AffineVariationalInequalities.h"
#include "SiconosSets.h"
#include "NumericsMatrix.h"
//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"
#include <stdlib.h>

void relay_avi_caoferris_test(RelayProblem* problem, double *z, double *w, int *info, SolverOptions* options)
{
  unsigned int n = problem->size;
  assert(n > 0);
  assert(problem->M);
  unsigned int s = 2*n;

  /* Copy the data from Relay problem to an AVI struct */
  AffineVariationalInequalities avi_pb;
  avi_pb.size = n;
  avi_pb.M = problem->M;
  avi_pb.q = problem->q;
  polyhedron poly;
  avi_pb.poly.split = &poly;

  poly.id = SICONOS_SET_POLYHEDRON;
  poly.size_ineq = s;
  poly.size_eq = 0;
  poly.H = NM_create_from_data(NM_DENSE, s, n, calloc(s*n, sizeof(double)));
  double* H = poly.H->matrix0;
  poly.K = (double *)malloc(s*sizeof(double));
  poly.Heq = NULL;
  poly.Keq = NULL;

  DEBUG_PRINT_VEC(problem->lb, n);
  DEBUG_PRINT_VEC(problem->ub, n);

  int starting_constraint = rand() % s;
  for (unsigned i = 0, j = starting_constraint; i < s; ++i, j = (j+1) % s)
  {
    if (j >= n)
    {
      H[i + s*(j-n)] = 1.0;
      poly.K[i] = problem->lb[j-n];
    }
    else
    {
      H[i + s*j] = -1.0;
      poly.K[i] = -problem->ub[j];
    }
  }
  DEBUG_PRINT("H matrix\n");
  DEBUG_EXPR_WE(for (unsigned i = 0; i < s; ++i)
      { for(unsigned j = 0 ; j < n; ++j)
      { DEBUG_PRINTF("% 2.2e ", H[i + j*s]) }
      DEBUG_PRINT("\n")});

  DEBUG_PRINT("K vector\n");
  DEBUG_EXPR_WE(for (unsigned i = 0; i < s; ++i)
      { DEBUG_PRINTF("% 2.2e ", poly.K[i])
      DEBUG_PRINT("\n")});

 /* Call directly the 3rd stage 
   * Here w is used as u and z as s in the AVI */
  *info = avi_caoferris(&avi_pb, z, w, options);

  free_polyhedron(&poly);
}

int relay_avi_caoferris_test_setDefaultSolverOptions(SolverOptions* options)
{
  solver_options_set(options, SICONOS_RELAY_AVI_CAOFERRIS_TEST);
  return 0;
}

