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
/*!\file lcp_path_solver.c

  This subroutine allows the resolution of LCP (Linear Complementary Problem).\n
  Try \f$(z,w)\f$ such that:\n
  \f$
   \left\lbrace
    \begin{array}{l}
    A u + Cv +a =0\\
    D u + Bv +b = w
    0 \le v \perp  w \ge 0\\
    \end{array}
   \right.
  \f$

  where M is an (\f$nn \times nn\f$)-matrix, q , w and z nn-vectors.
*/
/*!\fn  void lcp_path_solver( int *nn , double *M , double *q , double *z , double *w , int *info , int *iparamLCP , double *dparamLCP )

  lcp_pgs (Projected Gauss-Seidel) is a basic Projected Gauss-Seidel solver for LCP.\n

  \param nn      On enter, an integer which represents the dimension of the system.
  \param M     On enter, a (\f$nn \times nn\f$)-vector of doubles which contains the components of the matrix with a fortran storage.
  \param q       On enter, a nn-vector of doubles which contains the components of the right hand side vector.
  \param z       On return, a nn-vector of doubles which contains the solution of the problem.
  \param w       On return, a nn-vector of doubles which contains the solution of the problem.
  \param info    On return, an integer which returns the termination value:\n
                 0 : convergence\n
                 1 : iter = itermax\n
                 2 : negative diagonal term

  \param iparamLCP  On enter/return a vector of integers:\n
                - iparamLCP[0] = itermax On enter, the maximum number of iterations allowed.
                - iparamLCP[1] = verbose  On enter, the output log identifiant:\n
                        0 : no output\n
                        >0: active screen output\n
                - iparamLCP[2] = it_end  On enter, the number of iterations performed by the algorithm.

  \param dparamLCP  On enter/return a vector of doubles:\n
                - dparamLCP[0] = tol     On enter, the tolerance required.
                - dparamLCP[1] = omega   On enter, the relaxation parameter (not yet available).
                - dparamLCP[2] = res     On return, the final error value.

  \author Olivier Bonnefon

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "LA.h"
#include <math.h>
#include "mlcp_solvers.h"

/*#define PATH_SOLVER*/

#ifdef PATH_SOLVER
#include "path/SimpleLCP.h"
#endif /*PATH_SOLVER*/

void mlcp_path(int *nn , int* mm, double *A , double *B , double *C , double *D , double *a, double *b, double *u, double *v, double *w , int *info ,  int *iparamMLCP , double *dparamMLCP)
{
  int n = *nn;
  int m = *mm;
  int dim = (n + m);
  int i = 0, j = 0;
  int *m_i = 0;
  int *m_j = 0;
  double *m_ij = 0, *lb = 0, *ub = 0, *z;
  double *M, *q;
  double tol ;
  int nnz = 0;
  double err;
  int verbose = 1;
  double val = 0;
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
