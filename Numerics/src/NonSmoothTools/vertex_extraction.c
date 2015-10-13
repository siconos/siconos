/* Siconos-Numerics, Copyright INRIA 2005-2015
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
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

#include "SiconosSets.h"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"

#ifdef WITH_LPSOLVE
#include "lp_lib.h"

void siconos_find_vertex(const polyhedron* P, unsigned size, int* basis)
{
  unsigned nrows = P->size_ineq;
  double* restrict H = P->H;
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
    printf("find_vertex_lpsolve: failure in the LP solver\n");
    exit(EXIT_FAILURE);
  }

  get_basis(lp, basis, FALSE);
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

#endif
