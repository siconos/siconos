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
#include "SiconosConfig.h"
#include "MLCP_Solvers.h"
#include "SiconosCompat.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "NumericsMatrix.h"
#include "numerics_verbose.h"
#ifdef HAVE_PATHFERRIS
#include "tools/InterfaceToPathFerris/SimpleLCP.h"
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
