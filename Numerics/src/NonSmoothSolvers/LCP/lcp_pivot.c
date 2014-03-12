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
#include <float.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "LCP_Solvers.h"
#include "pivot-utils.h"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"

int pivot_selection_bard(double** mat, unsigned int dim)
{
  int block = -1;
  double zb, z0, dblock;

  for (unsigned int i = 0; i < dim; ++i)
  {
    zb = mat[i][0];
    if (zb < 0.0)
    {
      if (block == -1)
      {
        z0    = zb;
        block = i;
      }
      else
      {
        for (unsigned int j = 1; j <= dim; ++j)
        {
          dblock = mat[i][j]/zb - mat[block][j]/z0;
          if (dblock < 0.0)
          {
            break;
          }
          else if (dblock > 0.0)
          {
            block = i;
            break;
          }
        }
      }
    }
  }
  return block;
}

int pivot_selection_least_index(double** mat, unsigned int dim)
{
  int block = -1;

  for (unsigned int i = 0; i < dim; ++i)
  {
    if (mat[i][0] < 0.0)
    {
      block = i;
      break;
    }
  }
  return block;
}

void init_M_bard(double** mat, double* M, unsigned int dim, double* q)
{
  /* construction of mat matrix such that
   * mat = [ q | Id | -M ]
   */

  /*  Copy M but mat[dim+1:, :] = -M */
  for (unsigned int i = 0; i < dim; ++i)
  {
    for (unsigned int j = 1; j <= dim; ++j)
    {
      mat[i][j] = 0.0; /* We need to init only the part corresponding to Id */
      mat[i][j + dim] = -M[dim*(j-1) + i]; /* Siconos is in column major */
    }
  }

  for (unsigned int i = 0; i < dim; ++i)
  {
    mat[i][0] = q[i];
    mat[i][i + 1] =  1.0;
  }
}

void init_M_least_index(double** mat, double* M, unsigned int dim, double* q)
{
  /* construction of mat matrix such that
   * mat = [ q | M ]
   */

  /* We need to init only the part corresponding to Id */
  for (unsigned int i = 0 ; i < dim; ++i)
  {
    mat[i][0] = q[i];
    for (unsigned int j = 0 ; j < dim; ++j)
      mat[i][j+1] = M[dim*j + i];
  }
}

void lcp_pivot(LinearComplementarityProblem* problem, double* u , double* s, int *info , SolverOptions* options)
{
  /* matrix M of the LCP */
  assert(problem);
  assert(problem->M);
  double * M = problem->M->matrix0;
  assert(M);
  /* size of the LCP */

  unsigned int dim = problem->size;
  assert(dim>0);
  unsigned int dim2;

  unsigned int drive = dim+1;
  int block = -1;
  unsigned int has_sol = 0;
  unsigned int nb_iter = 0;
  unsigned int leaving = 0;
  unsigned int itermax = options->iparam[0];
  unsigned int preAlloc = options->iparam[2];
  unsigned int pivot_selection_rule = options->iparam[3];

  double pivot;
  double tmp;
  unsigned int* basis;
  double** mat;

  /*output*/

  options->iparam[1] = 0;

  /* Allocation */
  basis = (unsigned int *)malloc(dim * sizeof(unsigned int));
  mat = (double **)malloc(dim * sizeof(double*));

  switch (pivot_selection_rule)
  {
    case SICONOS_LCP_PIVOT_BARD:
      dim2 = 2*dim + 1;
      break;
    case SICONOS_LCP_PIVOT_LEAST_INDEX:
      dim2 = dim + 1;
      break;
    case SICONOS_LCP_PIVOT_LEMKE:
    default:
      dim2 = 2 * (dim + 1);
  }

  for (unsigned int i = 0 ; i < dim; ++i)
    mat[i] = (double *)malloc(dim2 * sizeof(double));

  assert(problem->q);

  switch (pivot_selection_rule)
  {
    case SICONOS_LCP_PIVOT_BARD:
      init_M_bard(mat, M, dim, problem->q);
      break;
    case SICONOS_LCP_PIVOT_LEAST_INDEX:
      init_M_least_index(mat, M, dim, problem->q);
      break;
    case SICONOS_LCP_PIVOT_LEMKE:
    default:
      init_M_lemke(mat, M, dim, dim2, dim, problem->q, NULL); /* TODO support custom covering vector */
  }

  DEBUG_PRINT("total matrix\n");
  DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { for(unsigned int j = 0 ; j < dim2; ++j)
      { DEBUG_PRINTF("% 2.2e ", mat[i][j]) }
      DEBUG_PRINT("\n")});
  /* End of construction of mat */
  for (unsigned int i = 0 ; i < dim  ; ++i) basis[i] = i + 1;
  /* Looking for pivot */
  switch (pivot_selection_rule)
  {
    case SICONOS_LCP_PIVOT_BARD:
      block = pivot_selection_bard(mat, dim);
      drive = block + dim + 1;
      break;
    case SICONOS_LCP_PIVOT_LEAST_INDEX:
      block = pivot_selection_least_index(mat, dim);
      drive = block + 1;
      break;
    case SICONOS_LCP_PIVOT_LEMKE:
    default:
      block = pivot_init_lemke(mat, dim);
  }

  if (block == -1)
  {
    /** exit, the solution is at hand with the current basis */
    DEBUG_PRINT("Trivial solution\n");
    goto exit_lcp_pivot;
  }

  /* Pivot < mu , drive >  or < drive, drive > */

  DEBUG_PRINTF("Pivoting %i and %i\n", block, drive);
  pivot = mat[block][drive];

  if (fabs(pivot) < DBL_EPSILON)
  {
    *info = 3;
    if (verbose > 0)
      printf("the pivot is nul, the algorithm cannot be used !\n");
    goto exit_lcp_pivot;
  }

  switch (pivot_selection_rule)
  {
    case SICONOS_LCP_PIVOT_BARD:
      do_pivot_driftless(mat, dim, dim2, block, drive);
      break;
    case SICONOS_LCP_PIVOT_LEAST_INDEX:
      do_pivot(mat, dim, dim2, block, drive);
      break;
    case SICONOS_LCP_PIVOT_LEMKE:
    default:
      do_pivot_driftless(mat, dim, dim2, block, drive);
  }

  /* Update the basis */
  switch (pivot_selection_rule)
  {
    /* Principal Pivoting Methods  */
    case SICONOS_LCP_PIVOT_BARD:
      basis[block] = basis[block] <= dim ? block + dim + 2 : block + 1;
      break;
    case SICONOS_LCP_PIVOT_LEAST_INDEX:
      basis[block] = basis[block] <= dim ? block + dim + 2 : block + 1;
      break;
    case SICONOS_LCP_PIVOT_LEMKE:
    default:
      /** one basic u is leaving and mu enters the basis */
      leaving = basis[block];
      basis[block] = drive;
  }

  DEBUG_EXPR_WE( DEBUG_PRINT("new basis: ")
      for (unsigned int i = 0; i < dim; ++i)
      { DEBUG_PRINTF("%i ", basis[i])}
      DEBUG_PRINT("\n"));

  DEBUG_PRINT("total matrix\n");
  DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { for(unsigned int j = 0 ; j < dim2; ++j)
      { DEBUG_PRINTF("% 2.2e ", mat[i][j]) }
      DEBUG_PRINT("\n")});

  while (nb_iter < itermax && !has_sol)
  {

    ++nb_iter;

    /* Start research of argmin lexico for minimum ratio test */

    /* Looking for pivot */
    switch (pivot_selection_rule)
    {
      case SICONOS_LCP_PIVOT_BARD:
        block = pivot_selection_bard(mat, dim);
        drive = block + dim + 1;
        break;
      case SICONOS_LCP_PIVOT_LEAST_INDEX:
        block = pivot_selection_least_index(mat, dim);
        drive = block + 1;
        break;
      case SICONOS_LCP_PIVOT_LEMKE:
      default:
        if (leaving < dim + 1)
        {
          drive = leaving + dim + 1;
        }
        else if (leaving > dim + 1)
        {
          drive = leaving - (dim + 1);
        }
        block = pivot_selection_lemke(mat, dim, drive);
    }

    /* We stop here: it either mean that the algorithm stops here or that there
     * is an issue with the LCP */
    if (block == -1)
    {
      if (pivot_selection_rule == SICONOS_LCP_PIVOT_LEMKE || pivot_selection_rule == 0)
      {
        DEBUG_PRINT("The pivot column is nonpositive !\n"
          "It either means that the algorithm failed or that the LCP is infeasible\n"
          "Check the class of the M matrix to find out the meaning of this\n");
      }
      break;
    }

    DEBUG_PRINTF("driving variable %i \n", drive);
    if (basis[block] == dim + 1) has_sol = 1;

    /* Pivot < block , drive > */
    DEBUG_PRINTF("Pivoting %i and %i\n", block, drive);

    pivot = mat[block][drive];
    if (fabs(pivot) < DBL_EPSILON)
    {
      *info = 3;
      if (verbose > 0)
        printf("the pivot is nul, the algorithm cannot be used !\n");
      goto exit_lcp_pivot;
    }

    switch (pivot_selection_rule)
    {
      case SICONOS_LCP_PIVOT_BARD:
        do_pivot_driftless(mat, dim, dim2, block, drive);
        break;
      case SICONOS_LCP_PIVOT_LEAST_INDEX:
        do_pivot(mat, dim, dim2, block, drive);
        break;
      case SICONOS_LCP_PIVOT_LEMKE:
      default:
        do_pivot_driftless(mat, dim, dim2, block, drive);
    }

    switch (pivot_selection_rule)
    {
      /* Principal Pivoting Methods  */
      case SICONOS_LCP_PIVOT_BARD:
        basis[block] = basis[block] <= dim ? block + dim + 2 : block + 1;
        break;
      case SICONOS_LCP_PIVOT_LEAST_INDEX:
        basis[block] = basis[block] <= dim ? block + dim + 2 : block + 1;
        break;
      case SICONOS_LCP_PIVOT_LEMKE:
      default:
        /** one basic variable is leaving and the driving one enters the basis */
        leaving = basis[block];
        basis[block] = drive;
    }

    DEBUG_EXPR_WE( DEBUG_PRINT("new basis: ")
      for (unsigned int i = 0; i < dim; ++i)
      { DEBUG_PRINTF("%i ", basis[i])}
      DEBUG_PRINT("\n"));

    DEBUG_PRINT("total matrix\n");
    DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { for(unsigned int j = 0 ; j < dim2; ++j)
      { DEBUG_PRINTF("% 2.2e ", mat[i][j]) }
      DEBUG_PRINT("\n")});
  } /* end while*/

  DEBUG_EXPR_WE( DEBUG_PRINT("final basis: ")
      for (unsigned int i = 0; i < dim; ++i)
      { DEBUG_PRINTF("%i ", basis[i])}
      DEBUG_PRINT("\n"));

  DEBUG_PRINT("total matrix\n");
  DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { for(unsigned int j = 0 ; j < dim2; ++j)
      { DEBUG_PRINTF("% 2.2e ", mat[i][j]) }
      DEBUG_PRINT("\n")});

exit_lcp_pivot:

  /* Recover solution */
  for (unsigned int i = 0 ; i < dim; ++i)
  {
    drive = basis[i];
//    assert(drive >0);
    if (drive < dim + 1)
    {
      u[drive - 1] = 0.0;
      s[drive - 1] = mat[i][0];
    }
    else if (drive > dim + 1)
    {
      u[drive - dim - 2] = mat[i][0];
      s[drive - dim - 2] = 0.0;
    }
  }

  DEBUG_PRINT("u s\n");
  DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { DEBUG_PRINTF("%e %e\n", u[i], s[i]) });

  options->iparam[1] = nb_iter;

  switch (pivot_selection_rule)
  {
    /* Principal Pivoting Methods  */
    case SICONOS_LCP_PIVOT_BARD:
    case SICONOS_LCP_PIVOT_LEAST_INDEX:
      *info = lcp_compute_error(problem, u, s, options->dparam[0], &tmp);
      break;
    case SICONOS_LCP_PIVOT_LEMKE:
    default:
      /** one basic variable is leaving and the driving one enters the basis */
      if (has_sol) *info = 0;
      else *info = 1;
  }

  if (info > 0)
  {
    DEBUG_PRINT("No solution found !\n");
  }

  free(basis);

  for (unsigned int i = 0 ; i < dim ; ++i) free(mat[i]);
  free(mat);
}

int linearComplementarity_pivot_setDefaultSolverOptions(SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the generic pivot Solver\n");
  }

  options->solverId = SICONOS_LCP_PIVOT;
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
