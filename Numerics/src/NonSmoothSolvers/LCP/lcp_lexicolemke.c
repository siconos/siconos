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

#include <math.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "LCP_Solvers.h"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"


void lcp_lexicolemke(LinearComplementarityProblem* problem, double *zlem , double *wlem , int *info , SolverOptions* options)
{
  /* matrix M of the lcp */
  double * M = problem->M->matrix0;
  assert(M);
  /* size of the LCP */
  int dim = problem->size;
  assert(dim>0);
  int dim2 = 2 * (dim + 1);

  int i, drive, block, Ifound;
  int ic, jc;
  int ITER;
  int nobasis;
  int itermax = options->iparam[0];

  i=0;
  int n = problem->size;
  double *q = problem->q;
  
  while ((i < (n - 1)) && (q[i] >= 0.)) 
    i++;
  
  if ((i == (n - 1)) && (q[n - 1] >= 0.))
  {

    /* TRIVIAL CASE : q >= 0
     * z = 0 and w = q is solution of LCP(q,M)
     */
    for (int j = 0 ; j < n; j++)
    {
      zlem[j] = 0.0;
      wlem[j] = q[j];
    }
    *info = 0;
    options->iparam[1] = 0;   /* Number of iterations done */
    options->dparam[1] = 0.0; /* Error */
    if (verbose > 0)
      printf("lcp_lexicolemke: found trivial solution for the LCP (positive vector q => z = 0 and w = q). \n");
    return ;
  }
  
  double z0, zb, dblock;
  double pivot, tovip;
  double tmp;
  int *basis;
  double** A;

  /*output*/

  options->iparam[1] = 0;

  /* Allocation */

  basis = (int *)malloc(dim * sizeof(int));
  A = (double **)malloc(dim * sizeof(double*));

  for (ic = 0 ; ic < dim; ++ic)
    A[ic] = (double *)malloc(dim2 * sizeof(double));

  /* construction of A matrix such that
   * A = [ q | Id | -d | -M ] with d = (1,...1)
   */

  /* We need to init only the part corresponding to Id */
  for (ic = 0 ; ic < dim; ++ic)
    for (jc = 1 ; jc <= dim; ++jc)
      A[ic][jc] = 0.0;

  for (ic = 0 ; ic < dim; ++ic)
    for (jc = 0 ; jc < dim; ++jc)
      A[ic][jc + dim + 2] = -M[dim * jc + ic];

  assert(problem->q);

  for (ic = 0 ; ic < dim; ++ic) A[ic][0] = problem->q[ic];

  for (ic = 0 ; ic < dim; ++ic) A[ic][ic + 1 ] =  1.0;
  for (ic = 0 ; ic < dim; ++ic) A[ic][dim + 1] = -1.0;

  DEBUG_PRINT("total matrix\n");
  DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { for(unsigned int j = 0 ; j < dim2; ++j)
      { DEBUG_PRINTF("%1.2e ", A[i][j]) }
      DEBUG_PRINT("\n")});
  /* End of construction of A */

  Ifound = 0;


  for (ic = 0 ; ic < dim  ; ++ic) basis[ic] = ic + 1;

  drive = dim + 1;
  block = 0;
  z0 = A[block][0];
  ITER = 0;

  /* Start research of argmin lexico */
  /* With this first step the covering vector enter in the basis */

  for (ic = 1 ; ic < dim ; ++ic)
  {
    zb = A[ic][0];
    if (zb < z0)
    {
      z0    = zb;
      block = ic;
    }
    else if (zb == z0)
    {
      for (jc = 0 ; jc < dim ; ++jc)
      {
        dblock = A[block][1 + jc] - A[ic][1 + jc];
        if (dblock < 0)
        {
          break;
        }
        else if (dblock > 0)
        {
          block = ic;
          break;
        }
      }
    }
  }

  /* Stop research of argmin lexico */
  DEBUG_PRINTF("Pivoting %i and %i\n", block, drive);

  pivot = A[block][drive];
  tovip = 1.0 / pivot;

  /* Pivot < block , drive > */

  A[block][drive] = 1;
  for (ic = 0       ; ic < drive ; ++ic) A[block][ic] = A[block][ic] * tovip;
  for (ic = drive + 1 ; ic < dim2  ; ++ic) A[block][ic] = A[block][ic] * tovip;

  /* */

  for (ic = 0 ; ic < block ; ++ic)
  {
    tmp = A[ic][drive];
    for (jc = 0 ; jc < dim2 ; ++jc) A[ic][jc] -=  tmp * A[block][jc];
  }
  for (ic = block + 1 ; ic < dim ; ++ic)
  {
    tmp = A[ic][drive];
    for (jc = 0 ; jc < dim2 ; ++jc) A[ic][jc] -=  tmp * A[block][jc];
  }

   nobasis = basis[block];
  basis[block] = drive;

  DEBUG_EXPR_WE( DEBUG_PRINT("new basis: ")
      for (unsigned int i = 0; i < dim; ++i)
      { DEBUG_PRINTF("%i ", basis[i])}
      DEBUG_PRINT("\n"));
  DEBUG_PRINT("total matrix\n");
  DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { for(unsigned int j = 0 ; j < dim2; ++j)
      { DEBUG_PRINTF("%1.2e ", A[i][j]) }
      DEBUG_PRINT("\n")});

  while (ITER < itermax && !Ifound)
  {

    ++ITER;

    if (nobasis < dim + 1)      drive = nobasis + (dim + 1);
    else if (nobasis > dim + 1) drive = nobasis - (dim + 1);

    DEBUG_PRINTF("driving variable %i \n", drive);

    /* Start research of argmin lexico for minimum ratio test */
    pivot = 1e20;
    block = -1;

    for (ic = 0 ; ic < dim ; ++ic)
    {
      zb = A[ic][drive];
      if (zb > 0.0)
      {
        z0 = A[ic][0] / zb;
        if (z0 > pivot) continue;
        if (z0 < pivot)
        {
          pivot = z0;
          block = ic;
        }
        else
        {
          for (jc = 1 ; jc < dim + 1 ; ++jc)
          {
            assert(block >=0 && "lcp_lexicolemke: block <0");
            dblock = A[block][jc] / pivot - A[ic][jc] / zb;
            if (dblock < 0) break;
            else if (dblock > 0)
            {
              block = ic;
              break;
            }
          }
        }
      }
    }
    if (block == -1)
    {
      Ifound = 1;
      DEBUG_PRINT("The pivot column is nonpositive !\n"
          "It either means that the algorithm failed or that the LCP is infeasible\n"
          "Check the class of the M matrix to find out the meaning of this\n");
      break;
    }

    if (basis[block] == dim + 1) Ifound = 1;

    /* Pivot < block , drive > */

    pivot = A[block][drive];
    tovip = 1.0 / pivot;
    A[block][drive] = 1;

    for (ic = 0       ; ic < drive ; ++ic) A[block][ic] = A[block][ic] * tovip;
    for (ic = drive + 1 ; ic < dim2  ; ++ic) A[block][ic] = A[block][ic] * tovip;

    /* */

    for (ic = 0 ; ic < block ; ++ic)
    {
      tmp = A[ic][drive];
      for (jc = 0 ; jc < dim2 ; ++jc) A[ic][jc] -=  tmp * A[block][jc];
    }
    for (ic = block + 1 ; ic < dim ; ++ic)
    {
      tmp = A[ic][drive];
      for (jc = 0 ; jc < dim2 ; ++jc) A[ic][jc] -=  tmp * A[block][jc];
    }

    nobasis = basis[block];
    basis[block] = drive;

    DEBUG_EXPR_WE( DEBUG_PRINT("new basis: ")
      for (unsigned int i = 0; i < dim; ++i)
      { DEBUG_PRINTF("%i ", basis[i])}
      DEBUG_PRINT("\n"));

    DEBUG_PRINT("total matrix\n");
    DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { for(unsigned int j = 0 ; j < dim2; ++j)
      { DEBUG_PRINTF("%1.2e ", A[i][j]) }
      DEBUG_PRINT("\n")});

  } /* end while*/

  DEBUG_EXPR_WE( DEBUG_PRINT("new basis: ")
      for (unsigned int i = 0; i < dim; ++i)
      { DEBUG_PRINTF("%i ", basis[i])}
      DEBUG_PRINT("\n"));

  DEBUG_PRINT("total matrix\n");
  DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { for(unsigned int j = 0 ; j < dim2; ++j)
      { DEBUG_PRINTF("%1.2e ", A[i][j]) }
      DEBUG_PRINT("\n")});

  for (ic = 0 ; ic < dim; ++ic)
  {
    drive = basis[ic];
    if (drive < dim + 1)
    {
      zlem[drive - 1] = 0.0;
      wlem[drive - 1] = A[ic][0];
    }
    else if (drive > dim + 1)
    {
      zlem[drive - dim - 2] = A[ic][0];
      wlem[drive - dim - 2] = 0.0;
    }
  }

  options->iparam[1] = ITER;

  if (Ifound) *info = 0;
  else *info = 1;

  free(basis);

  for (i = 0 ; i < dim ; ++i) free(A[i]);
  free(A);
}


int linearComplementarity_lexicolemke_setDefaultSolverOptions(SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the Lemke Solver\n");
  }

  options->solverId = SICONOS_LCP_LEMKE;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 5;
  options->dSize = 5;
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  options->dWork = NULL;
  options->iWork = NULL;   options->callback = NULL; options->numericsOptions = NULL;
  options->dparam[0] = 1e-6;
  options->iparam[0] = 10000;
  return 0;
}
