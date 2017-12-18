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

/*!\file vertex_extraction.c
 * \brief interface to various LP solvers
 */

#include "vertex_extraction.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#ifdef __cplusplus
#undef restrict
#define restrict __restrict
#endif

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"

#ifdef WITH_LPSOLVE
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#include "lp_lib.h"

#include "NumericsMatrix.h"


void siconos_find_vertex(const polyhedron* P, unsigned size, lapack_int* basis)
{
  unsigned nrows = P->size_ineq;
  assert(P->H->matrix0);
  double* restrict H = P->H->matrix0;
  double* restrict K = P->K;
  lprec *lp;
  lp = make_lp(nrows, nrows+size);
  set_verbose(lp, CRITICAL);
  set_minim(lp);

  int* rowno = (int*) malloc(nrows * sizeof(int));
  for (unsigned i = 0; i < nrows; ++i)
  {
    rowno[i] = i+1;
    set_mat(lp, i+1, size+i+1, -1.0);
    //assert(K[nrows-1-i] >= 0.0);
    set_mat(lp, 0, size+1+i, fabs(K[nrows-1-i]) + 0.0001*(rand()%100));
    set_constr_type(lp, i+1, EQ);
  }
//    set_constr_type(lp, nrows+1, GE);

  for (unsigned i = 1, j = 0; i <= size; ++i, j+= nrows)
  {
    set_columnex(lp, i, nrows, &H[j], rowno);
    set_unbounded(lp, i);
  }
  set_rh_vec(lp, &K[-1]);
//  set_mat(lp, nrows+1, size+1, 1.0);

#ifdef DEBUG_STDOUT
  print_lp(lp);
#endif
  /* Solve the LP */
  int info = solve(lp);
  if (info != 0)
  {
    printf("find_vertex_lpsolve: failure in the LP solver: info = %d\n", info);
    exit(EXIT_FAILURE);
  }

  int* lp_basis = NULL;
  if (sizeof(lapack_int) != sizeof(int))
  {
    lp_basis = (int*)malloc((nrows + 1)*sizeof(int));
  }
  else
  {
    lp_basis = (int*)basis;
  }

  get_basis(lp, lp_basis, FALSE);

  if (sizeof(lapack_int) != sizeof(int))
  {
    for (size_t i = 0; i < (nrows + 1); ++i)
    {
      basis[i] = (lapack_int)lp_basis[i];
    }
  }

#ifdef DEBUG_STDOUT
  print_solution(lp, 3);
#endif
  delete_lp(lp);
  free(rowno);
}

#else

void siconos_find_vertex(const polyhedron* P, unsigned size, int* basis)
{
  printf("You need to compile the lp_solve support");
  exit(EXIT_FAILURE);
}
#pragma  GCC diagnostic pop

#endif
