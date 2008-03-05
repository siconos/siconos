/* Siconos-Numerics version 3.0.0, Copyright INRIA 2005-2008.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
#include "LA.h"
#include "MLCP_Solvers.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*#define PATH_SOLVER*/

#ifdef PATH_SOLVER
#include "path/SimpleLCP.h"
#endif /*PATH_SOLVER*/

/*
Warning: this function requires MLCP with M and q, not (A,B,C,D).
The input structure MixedLinearComplementarity_Problem is supposed to fit with this form.
*/

void mlcp_path(MixedLinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options)
{
  double tol = options->dparam[0];
  *info = 1;
#ifdef PATH_SOLVER
  MCP_Termination termination;

  M = problem->M->matrix0;
  q = problem->q;

  nnz = nbNonNulElems(dim, M, 1.0e-18);
  m_i = (int *)calloc(nnz + 1, sizeof(int));
  m_j = (int *)calloc(nnz + 1, sizeof(int));
  m_ij = (double *)calloc(nnz + 1, sizeof(double));
  lb = (double *)calloc(dim + 1, sizeof(double));
  ub = (double *)calloc(dim + 1, sizeof(double));


  FortranToPathSparse(dim, M, 1.0e-18, m_i, m_j, m_ij);
  for (i = 0; i < n; i++)
  {
    lb[i] = -1e20;
    ub[i] = 1e20;
  }
  for (i = n; i < n + m; i++)
  {
    lb[i] = 0;
    ub[i] = 1e20;
  }
  printLCP(dim, nnz, m_i, m_j, m_ij, q, lb, ub);
  SimpleLCP(dim, nnz, m_i, m_j, m_ij, q, lb, ub,
            &termination, z);

  if (termination == MCP_Error)
  {
    *info = 1;
    if (verbose > 0)
      printf("PATH : Error in the solution.\n");
  }
  else if (termination == MCP_Solved)
  {
    for (i = 0; i < n; i++)
    {
      u[i] = z[i];
    }
    for (i = 0; i < m; i++)
    {
      v[i] = z[n + i];
      w[i] = q[n + i];
      for (j = 0; j < dim; j++)
      {
        w[i] += M[i + dim * j] * z[j];
      }
    }

    *info = 0;
    mlcp_compute_error(problem, z, w, tol, &err);

    if (verbose > 0)
      printf("PATH : MLCP Solved, error %10.7f.\n", err);
  }
  else
  {
    if (verbose > 0)
      printf("PATH : MLCP Other error: %d\n", termination);
  }




  free(m_i);
  free(m_j);
  free(m_ij);
  free(lb);
  free(ub);

#endif /*PATH_SOLVER*/


  return;
}
