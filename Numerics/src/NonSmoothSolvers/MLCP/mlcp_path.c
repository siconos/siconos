/* Siconos-Numerics, Copyright INRIA 2005-2011.
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
#include "LA.h"
#include "NumericsConfig.h"
#include "MLCP_Solvers.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef HAVE_PATHFERRIS
#include "InterfaceToPathFerris/SimpleLCP.h"
#endif
/*
Warning: this function requires MLCP with M and q, not (A,B,C,D).
The input structure MixedLinearComplementarityProblem is supposed to fit with this form.
*/
int mixedLinearComplementarity_path_setDefaultSolverOptions(MixedLinearComplementarityProblem* problem, SolverOptions* pSolver)
{
#ifdef HAVE_PATHFERRIS
  mixedLinearComplementarity_default_setDefaultSolverOptions(problem, pSolver);
#endif
  return 0;
}


void mlcp_path(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options)
{
  *info = 1;
#ifdef HAVE_PATHFERRIS
  *info = 0;
  MCP_Termination termination;
  double tol = options->dparam[0];

  double * M = problem->M->matrix0;
  double * q = problem->q;
  int nnz, i, j, n, m, dim, numLine;
  n = problem->n;
  m = problem->m;
  dim = m + n;
  /*  if (verbose){
    printf("initial values for z:\n");
    for (int i=0;i<dim;i++)
      printf("%.15e\n", z[i]);
    printf("initial values for w:\n");
    for (int i=0;i<dim;i++)
      printf("%.15e\n", w[i]);

      }*/
  nnz = nbNonNulElems(dim, M, 1.0e-18);
  int * m_i = (int *)calloc(nnz + 1, sizeof(int));
  int * m_j = (int *)calloc(nnz + 1, sizeof(int));
  double * m_ij = (double *)calloc(nnz + 1, sizeof(double));
  double * lb = (double *)calloc(dim + 1, sizeof(double));
  double * ub = (double *)calloc(dim + 1, sizeof(double));
  //  double * u = z;
  //  double * v = z+n;
  double err;



  FortranToPathSparse(dim, M, 1.0e-18, m_i, m_j, m_ij);
  if (problem->blocksRows)
  {
    int numBlock = 0;
    while (problem->blocksRows[numBlock] < n + m)
    {
      if (!problem->blocksIsComp[numBlock])
      {
        for (numLine = problem->blocksRows[numBlock] ; numLine < problem->blocksRows[numBlock + 1]; numLine++)
        {
          lb[numLine] = -1e20;
          ub[numLine] = 1e20;
        }
      }
      else
      {
        for (numLine = problem->blocksRows[numBlock] ; numLine < problem->blocksRows[numBlock + 1]; numLine++)
        {
          lb[numLine] = 0;
          ub[numLine] = 1e20;
        }
      }
      numBlock++;
    }
  }
  else
  {
    printf("DEPRECED MLCP INTERFACE\n");
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
  }
  if (verbose)
    printLCP(dim, nnz, m_i, m_j, m_ij, q, lb, ub);
  SimpleLCP(dim, nnz, m_i, m_j, m_ij, q, lb, ub,
            &termination, z);

  if (termination == MCP_Error)
  {
    *info = 1;
    if (verbose)
      printf("PATH : Error in the solution.\n");
  }
  else if (termination == MCP_Solved)
  {
    /*     for (i=0;i<n;i++){ */
    /*       u[i]=z[i]; */
    /*     } */
    if (problem->blocksRows)
    {
      int numBlock = 0;
      while (problem->blocksRows[numBlock] < n + m)
      {
        if (!problem->blocksIsComp[numBlock])
        {
          for (numLine = problem->blocksRows[numBlock] ; numLine < problem->blocksRows[numBlock + 1]; numLine++)
          {
            w[numLine] = 0;
          }
        }
        else
        {
          for (numLine = problem->blocksRows[numBlock] ; numLine < problem->blocksRows[numBlock + 1]; numLine++)
          {
            w[numLine] = -q[numLine];
            for (int jj = 0; jj < n + m; jj++)
            {
              w[numLine] += M[numLine + dim * jj] * z[jj];
            }
          }
        }
        numBlock++;
      }
    }
    else
    {
      for (i = 0; i < n; i++)
        w[i] = 0;
      for (i = n; i < n + m; i++)
      {
        w[i] = -q[i];
        for (j = 0; j < n + m; j++)
        {
          w[i] += M[i + dim * j] * z[j];
        }
      }
    }



    /*1e-7 because it is the default tol of path.*/
    mlcp_compute_error(problem, z, w, tol, &err);
    if (err > 1e-7)
    {
      printf("PATH : MLCP Solved, error %10.7f.\n", err);
      //*info = 1;
    }
    if (problem->blocksRows)
    {
      int numBlock = 0;
      while (problem->blocksRows[numBlock] < n + m)
      {
        if (problem->blocksIsComp[numBlock])
        {
          for (numLine = problem->blocksRows[numBlock] ; numLine < problem->blocksRows[numBlock + 1]; numLine++)
          {
            if (z[numLine] > w[numLine])
              w[numLine] = 0;
          }
        }
        numBlock++;
      }
    }
    else
    {
      for (i = 0; i < m; i++)
      {
        if (z[n + i] > w[n + i])
          w[n + i] = 0;
      }
    }

    if (verbose)
      printf("PATH : MLCP Solved, error %10.7f.\n", err);
  }
  else
  {
    if (verbose)
      printf("PATH : MLCP Other error: %d\n", termination);
  }




  free(m_i);
  free(m_j);
  free(m_ij);
  free(lb);
  free(ub);

#endif /*HAVE_PATHFERRIS*/


  return;
}
