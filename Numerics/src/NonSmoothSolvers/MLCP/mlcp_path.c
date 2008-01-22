/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2007.
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "LA.h"
#include <math.h>
#include "MLCP_Solvers.h"

/*#define PATH_SOLVER*/

#ifdef PATH_SOLVER
#include "path/SimpleLCP.h"
#endif /*PATH_SOLVER*/

void mlcp_path(int *nn , int* mm, double *A , double *B , double *C , double *D , double *a, double *b, double *u, double *v, double *w , int *info ,  int *iparamMLCP , double *dparamMLCP)
{
  double tol ;
  *info = 1;
  tol   = dparamMLCP[0];
#ifdef PATH_SOLVER
  MCP_Termination termination;

  M = (double*) malloc((dim * dim) * sizeof(double));
  q = (double*) malloc(dim * sizeof(double));
  z = (double*) calloc(dim, sizeof(double));
  ABCDtoM(n , m, A , B , C , D , a, b, M, q);


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
    mlcp_compute_error(nn, mm,  A , B , C , D , a , b, u, v, verbose, w,  &err);

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
  free(M);
  free(q);
  free(z);

#endif /*PATH_SOLVER*/


  return;
}
