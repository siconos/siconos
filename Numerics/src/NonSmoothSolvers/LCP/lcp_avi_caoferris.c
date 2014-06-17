/* Siconos-Numerics, Copyright INRIA 2005-2012.
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

/*!\file lcp_avi_caoferris.c
 \brief Solve an LCP by reformulating it as an AVI and the solver by Cao and
Ferris solves the subsequent AVI.
 \author Olivier Huber
*/

#include "LinearComplementarityProblem.h"
#include "avi_caoferris.h"
#include "AVI_Solvers.h"
#include "LCP_Solvers.h"

void lcp_avi_caoferris(LinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options)
{
  unsigned int n = problem->size;

  /* Copy the data from LCP problem */
  AffineVariationalInequalities avi_pb;
  avi_pb.size = n;
  avi_pb.M = problem->M;
  avi_pb.q = problem->q;
  avi_pb.d = (double *)malloc(n*sizeof(double));
  for (unsigned int i = 0; i<n; ++i) avi_pb.d[i] = -1;
  avi_pb.poly = NULL;

  /* Set of active constraint is trivial */
  unsigned int * A = (unsigned int *)malloc(n*sizeof(unsigned int));
  for (unsigned int i = 0; i<n; ++i) A[i] = i+1;

  /* Call directly the 3rd stage */
  *info = avi_caoferris_stage3(&avi_pb, w, z, n, A, options);

  /* free allocated stuff */
  free(A);
  free(avi_pb.d);
}

int linearComplementarity_avi_caoferris_setDefaultSolverOptions(SolverOptions* options)
{
  return avi_caoferris_setDefaultSolverOptions(options);
}

