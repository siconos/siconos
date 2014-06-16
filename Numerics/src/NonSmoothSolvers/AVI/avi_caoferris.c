/* Siconos-Numerics, Copyright INRIA 2005-2013.
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
#include "AVI_Solvers.h"
#include "pivot-utils.h"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"

void avi_caoferris_stage3(AffineVariationalInequalities* problem, double* restrict u , double* restrict s, unsigned int size_x, unsigned int* restrict A, int *info , SolverOptions* options)
{
  assert(size_x > 0);
  /* matrix M of the avi */
  assert(problem);
  assert(problem->M);
  double * M = problem->M->matrix0;
  assert(M);
  /* size of the AVI */

  unsigned int dim = problem->size;
  assert(dim>0);
  unsigned int dim2 = 2 * (dim + 1);

  unsigned int drive = dim+1;
  unsigned int drive_number = dim+1;
  int block = -1;
  unsigned int has_sol = 0;
  unsigned int nb_iter = 0;
  unsigned int leaving;
  unsigned int itermax = options->iparam[0];


  double z0, zb, dblock;
  double pivot, pivot_inv;
  double tmp;
  unsigned int* basis;
  unsigned int* u_indx;
  unsigned int* s_indx;
  double* mat;

  /*output*/

  options->iparam[1] = 0;

  /* Allocation */
  basis = (unsigned int *)malloc(dim * sizeof(unsigned int));
  mat = (double *)malloc(dim * dim2 * sizeof(double));

  u_indx = (unsigned int *)malloc(dim * sizeof(unsigned int));
  s_indx = (unsigned int *)malloc(dim * sizeof(unsigned int));


  /* construction of mat matrix such that
   * mat = [ q | Id | -d | -M ] with d_i = i in A ? 1 : 0
   */

  /* We need to init only the part corresponding to Id */
  for (unsigned int i = 0 ; i < dim; ++i)
    for (unsigned int j = 1 ; j <= dim; ++j)
      mat[i + j*dim] = 0.0;

  /*  Copy M but mat[dim+2:, :] = -M */
  for (unsigned int i = 0 ; i < dim; ++i)
    for (unsigned int j = 0 ; j < dim; ++j)
      mat[i + dim*(j + dim + 2)] = -M[dim*j + i]; // Siconos is in column major

  assert(problem->q);

  for (unsigned int i = 0 ; i < dim; ++i)
  {
    mat[i] = problem->q[i];
    mat[i + dim*(i + 1)] =  1.0;
  }

  /** Add covering vector */
  assert(problem->d != NULL);
  for (unsigned int i = 0; i < size_x  ; ++i) mat[i + dim*(dim + 1)] = problem->d[i];
  for (unsigned int i = size_x; i < dim; ++i) mat[i + dim*(dim + 1)] = 0.0;


  DEBUG_PRINT("total matrix\n");
  DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { for(unsigned int j = 0 ; j < dim2; ++j)
      { DEBUG_PRINTF("% 2.2e ", mat[i + j*dim]) }
      DEBUG_PRINT("\n")});
  /* End of construction of mat */


  unsigned int val_A = A[0];

  /** Contruct the basis and the index maps for u and s
   * basis = (u_A, s_I) and nonbasic variables are (u_I, s_A)*/
  for (unsigned int i = 0, indx_A = 0, indx_I = 0; i < dim; ++i)
  {
    if (i == val_A-1) // i is in A, u_i is basic and s_i is nonbasic
    {
      basis[indx_A] = val_A;
      u_indx[i] = indx_A + 1;
      s_indx[i] = dim2 - size_x + i;
      if (++indx_A < size_x)
        val_A = A[indx_A];
    }
    else // i is not in A (therefore in I), s_i is basic, u_i is nonbasic
    {
      basis[size_x+indx_I] = dim + 2 + i;
      u_indx[i] = dim + 2 + indx_I;
      s_indx[i] = size_x + indx_I + 1;
      ++indx_I;
    }
  }
  DEBUG_PRINT("basis u_indx s_indx\n");
  DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { DEBUG_PRINTF("%i %i %i\n", basis[i], u_indx[i], s_indx[i]) });

  /* Start research of argmin lexico
   * lexicographic order is simple: we just look for the min of index in case
   * of tie */
  /* With this first step the covering vector enter in the basis */

  block = pivot_init_lemke(mat, size_x);

  /* Stop research of argmin lexico */


  if (block == -1)
  {
    /** exit, the solution is at hand with the current basis */
    DEBUG_PRINT("Trivial solution\n");
    goto exit_caoferris;
  }

  /* Pivot < mu , driver > */

  DEBUG_PRINTF("Pivoting %i and %i\n", block, drive);
  pivot = mat[block + drive*dim];
  pivot_inv = 1.0/pivot;

  /* Update column mat[block, :] */
  mat[block + drive*dim] = 1;
  for (unsigned int i = 0        ; i < drive; ++i) mat[block + i*dim] *= pivot_inv;
  for (unsigned int i = drive + 1; i < dim2 ; ++i) mat[block + i*dim] *= pivot_inv;

  /* Update other columns*/
  for (unsigned int i = 0; i < block; ++i)
  {
    tmp = mat[i + drive*dim];
    for (unsigned int j = 0; j < dim2; ++j) mat[i + j*dim] -= tmp*mat[block + j*dim];
  }
  for (unsigned int i = block + 1; i < dim; ++i)
  {
    tmp = mat[i + drive*dim];
    for (unsigned int j = 0; j < dim2; ++j) mat[i + j*dim] -= tmp*mat[block + j*dim];
  }

  /** one basic u is leaving and mu enters the basis */
  leaving = basis[block];
  basis[block] = drive;

  DEBUG_EXPR_WE( DEBUG_PRINT("new basis: ")
      for (unsigned int i = 0; i < dim; ++i)
      { DEBUG_PRINTF("%i ", basis[i])}
      DEBUG_PRINT("\n"));

  DEBUG_PRINT("total matrix\n");
  DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { for(unsigned int j = 0 ; j < dim2; ++j)
      { DEBUG_PRINTF("% 2.2e ", mat[i + j*dim]) }
      DEBUG_PRINT("\n")});

  while (nb_iter < itermax && !has_sol)
  {

    ++nb_iter;

    if (leaving < dim + 1)
    {
      drive_number = leaving + dim + 1;
      drive = s_indx[leaving-1];
    }
    else if (leaving > dim + 1)
    {
      drive_number = leaving - (dim + 1);
      drive = u_indx[leaving - (dim + 2)];
    }

    DEBUG_PRINTF("driving variable %i \n", drive_number);
    /* Start research of argmin lexico for minimum ratio test */

    pivot = 1e20;
    block = -1;

    for (unsigned int i = 0 ; i < dim ; ++i)
    {
      zb = mat[i + drive*dim];
      if (zb > 0.0)
      {
        z0 = mat[i] / zb;
        if (z0 > pivot) continue;
        if (z0 < pivot)
        {
          pivot = z0;
          block = i;
        }
        else
        {
          for (unsigned int j = 1; j <= dim; ++j)
          {
            assert(block >=0 && "avi_caoferris: block <0");
            dblock = mat[block + j*dim] / pivot - mat[i + j*dim] / zb;
            if (dblock < 0) break;
            else if (dblock > 0)
            {
              block = i;
              break;
            }
          }
        }
      }
    }
    if (block == -1) break;

    if (basis[block] == dim + 1) has_sol = 1;

    /* Pivot < block , drive > */
    DEBUG_PRINTF("Pivoting %i and %i\n", block, drive);

    pivot = mat[block + drive*dim];
    pivot_inv = 1.0/pivot;

    /* Update column mat[block, :] */
    mat[block + drive*dim] = 1;
    for (unsigned int i = 0        ; i < drive; ++i) mat[block + i*dim] *= pivot_inv;
    for (unsigned int i = drive + 1; i < dim2 ; ++i) mat[block + i*dim] *= pivot_inv;

    /* Update other columns*/
    for (unsigned int i = 0; i < block; ++i)
    {
      tmp = mat[i + drive*dim];
      for (unsigned int j = 0; j < dim2; ++j) mat[i + j*dim] -= tmp*mat[block + j*dim];
    }
    for (unsigned int i = block + 1; i < dim; ++i)
    {
      tmp = mat[i + drive*dim];
      for (unsigned int j = 0; j < dim2; ++j) mat[i + j*dim] -= tmp*mat[block + j*dim];
    }

    /** one basic variable is leaving and driver enters the basis */
    leaving = basis[block];
    basis[block] = drive_number;

  DEBUG_EXPR_WE( DEBUG_PRINT("new basis: ")
      for (unsigned int i = 0; i < dim; ++i)
      { DEBUG_PRINTF("%i ", basis[i])}
      DEBUG_PRINT("\n"));

  DEBUG_PRINT("total matrix\n");
  DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { for(unsigned int j = 0 ; j < dim2; ++j)
      { DEBUG_PRINTF("% 2.2e ", mat[i + j*dim]) }
      DEBUG_PRINT("\n")});
  } /* end while*/

  DEBUG_EXPR_WE( DEBUG_PRINT("new basis: ")
      for (unsigned int i = 0; i < dim; ++i)
      { DEBUG_PRINTF("%i ", basis[i])}
      DEBUG_PRINT("\n"));

  DEBUG_PRINT("total matrix\n");
  DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { for(unsigned int j = 0 ; j < dim2; ++j)
      { DEBUG_PRINTF("% 2.2e ", mat[i + j*dim]) }
      DEBUG_PRINT("\n")});

exit_caoferris:

  /* Recover solution */
  for (unsigned int i = 0 ; i < dim; ++i)
  {
    drive = basis[i];
    if (drive < dim + 1)
    {
      s[drive - 1] = 0.0;
      u[drive - 1] = mat[i];
    }
    else if (drive > dim + 1)
    {
      s[drive - dim - 2] = mat[i];
      u[drive - dim - 2] = 0.0;
    }
  }

  DEBUG_PRINT("u s\n");
  DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { DEBUG_PRINTF("%e %e\n", u[i], s[i]) });

  options->iparam[1] = nb_iter;

  if (has_sol) *info = 0;
  else *info = 1;

  free(basis);
  free(u_indx);
  free(s_indx);

  free(mat);
}


int avi_caoferris_setDefaultSolverOptions(SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the Cao-Ferris Solver\n");
  }

  options->solverId = SICONOS_AVI_CAOFERRIS;

  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 5;
  options->dSize = 5;
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  options->dWork = NULL;
  options->iWork = NULL;
  options->callback = NULL;
  options->numericsOptions = NULL;
  options->dparam[0] = 1e-6;
  options->iparam[0] = 10000;


  return 0;
}
