/* Siconos-Numerics, Copyright INRIA 2005-2014
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

#include "SiconosLapack.h"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"

#include "lcp_pivot.h"

int pivot_selection_bard(double* mat, unsigned int dim)
{
  int block = -1;
  double zb, z0, dblock;

  for (unsigned int i = 0; i < dim; ++i)
  {
    zb = mat[i];
    if (zb < 0.0)
    {
      if (block == -1)
      {
        z0    = zb;
        block = i;
      }
      else
      {
        for (unsigned int j = dim; j <= dim*dim; j += dim)
        {
          dblock = mat[i + j]/zb - mat[block + j]/z0;
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

int pivot_selection_least_index(double* mat, unsigned int dim)
{
  int block = -1;

  for (unsigned int i = 0; i < dim; ++i)
  {
    if (mat[i] < 0.0)
    {
      block = i;
      break;
    }
  }
  return block;
}

void init_M_bard(double* restrict mat, double* restrict M, unsigned int dim, double* restrict q)
{
  /* construction of mat matrix such that
   * mat = [ q | Id | -M ]
   */

  /*  Copy M but mat[dim+1:, :] = -M */
  for (unsigned int i = 0; i < dim; ++i)
  {
    for (unsigned int j = 0; j < dim*dim; j += dim)
    {
      mat[i + j + dim*(1 + dim)] = -M[j + i]; /* Siconos is in column major */
    }
  }

  memset(&mat[dim], 0, sizeof(double) * dim * dim);

  for (unsigned int i = 0, j = dim; i < dim; ++i, j+= dim)
  {
    mat[i] = q[i];
    mat[i + j] =  1.0;
  }
}

void init_M_least_index(double* restrict mat, double* restrict M, unsigned int dim, double* restrict q)
{
  /* construction of mat matrix such that
   * mat = [ q | M ]
   */

  /* We need to init only the part corresponding to Id */
  for (unsigned int i = 0 ; i < dim; ++i)
  {
    mat[i] = q[i];
    for (unsigned int j = 0 ; j < dim*dim; j += dim)
      mat[i + j + dim] = M[j + i];
  }
}

int init_M_lemke_warm_start(int n, double* restrict u, double* restrict mat, double* restrict M, double* restrict q, int* restrict basis, double* restrict cov_vec)
{
  /* crash the basis and form the matrix for Lemke's algorithm
   * mat = [ q | Id | -d | -M ]
   * TODO: implement a version with no memory allocation
   * TODO: think if it is necessary to compute the new covering vector if
   * cov_vec == NULL
   */

  /* q vector */
  double* q_bar = mat;
  cblas_dcopy(n, q, 1, q_bar, 1);

  /* covering vector for the auxiliary variable */
  double* d = &mat[(n+1)*n];
  if (cov_vec)  cblas_dcopy(n, cov_vec, 1, d, 1);
  else for (unsigned int i = 0; i < n; ++i) d[i] = 1.0;

  /* take care of M */
  double* mat_basic = &mat[n];
  double* mat_nonbasic = &mat[(n+2)*n];
  for (unsigned int i = 0; i < n; ++i)
  {
    if (u[i] > DBL_EPSILON) // M_bas[:, i] = M[:, i]
    {
      basis[i] = i + 2 + n;
      cblas_dcopy(n, &M[i*n], 1, &mat_basic[i*n], 1);
      memset(&mat_nonbasic[i*n], 0, sizeof(double) * n);
      mat_nonbasic[i*n + i] = -1.0;
    }
    else /* s[i] > 0.0 and if both s and u are nonbasic, we insert a column from the identity matrix
          * this choice could ne different
          * M_bas[:, i] = -I[:, i] */
    {
      basis[i] = i + 1;
      cblas_dcopy(n, &M[i*n], 1, &mat_nonbasic[i*n], 1);
      memset(&mat_basic[i*n], 0, sizeof(double) * n);
      mat_basic[i*n + i] = -1.0;
    }
  }

  /* data for LAPACKE */
  int *ipiv = (int*)malloc((n+1)*sizeof(int));
  int info = 0;

  /* Compute LU factorisation of basis */
  DGETRF(n, n, mat_basic, n, ipiv, &info);

  assert(info <= 0 && "crash_pivot_basis :: pivot form DGETRF > 0, this should not append !\n");
  if (info < 0)
  {
    printf("crash_pivot_basis :: the crash basis is singular, cannot inverse the matrix.\n\
            The (first) diagonal element of U to be 0.0 is %d\n\
            A remedy remains to be implemented !\n", info);
    return info;
  }

  /* Compute new values Mbar = (M_bas)^{-1} M_nonbas ; q_bar = (M_bas)^{-1} q ;
   * d_bar = (M_bas)^{-1} d */
  DGETRS(LA_NOTRANS, n, n, mat_basic, n, ipiv, mat_nonbasic, n, &info);
  DGETRS(LA_NOTRANS, n, 1, mat_basic, n, ipiv, q_bar, n, &info);
  DGETRS(LA_NOTRANS, n, 1, mat_basic, n, ipiv, d, n, &info);

  /* take the opposite of d and M, see matrix definition */
  cblas_dscal(n, -1.0, q_bar, 1);
//  cblas_dscal(n, -1.0, d, 1);
//  cblas_dscal(n*n, -1.0, mat_nonbasic, 1);

  /* set the identity part in the matrix (for the lexicographic ordering) */
  memset(mat_basic, 0, sizeof(double) * n * n);
  for (unsigned int i = 0; i < n; ++i) mat[i + n*(i + 1)] =  1.0;

  free(ipiv);
  return info;
}

void lcp_pivot(LinearComplementarityProblem* problem, double* u , double* s, int *info , SolverOptions* options)
{
  lcp_pivot_covering_vector(problem, u, s, info, options, NULL);
}


void lcp_pivot_covering_vector(LinearComplementarityProblem* problem, double* restrict u , double* restrict s, int *info , SolverOptions* options, double* restrict cov_vec)
{
  /* matrix M of the LCP */
  assert(problem);
  assert(problem->M);
  double* M = problem->M->matrix0;
  assert(M);
  assert(problem->q);


  unsigned int dim = problem->size;
  assert(dim>0);
  unsigned int dim2;
  /* size of the LCP */
  DEBUG_EXPR_WE( DEBUG_PRINT("matrix M: ") display(problem->M); DEBUG_PRINT("vector q: ")
      for(unsigned i = 0; i < dim; ++i) {printf("%e ", problem->q[i]);} printf("\n");
      if (cov_vec) { DEBUG_PRINT("covering vector: ") for(unsigned i = 0; i < dim; ++i) {printf("%e ", cov_vec[i]);}printf("\n");});

  unsigned int drive = dim+1;
  int bck_drive = -1;
  int block = -1;
  unsigned int has_sol = 0;
  unsigned int nb_iter = 0;
  unsigned int leaving = 0;
  unsigned int itermax = options->iparam[0];
  unsigned int preAlloc = options->iparam[SICONOS_IPARAM_PREALLOC];
  unsigned int pivot_selection_rule = options->iparam[SICONOS_IPARAM_PIVOT_RULE];

  double pivot;
  double tmp;
  int* basis;
  int basis_init = 0; /* 0 is basis was not initialized, 1 otherwise*/
  unsigned t_indx = 0;
  double* t_stack = NULL;
  double* mat;

  *info = 0;

  /*output*/

  options->iparam[1] = 0;

  /* Allocation */
  switch (pivot_selection_rule)
  {
    case SICONOS_LCP_PIVOT_BARD:
      dim2 = 2*dim + 1;
      break;
    case SICONOS_LCP_PIVOT_LEAST_INDEX:
      dim2 = dim + 1;
      break;
    case SICONOS_LCP_PIVOT_LEMKE:
    case SICONOS_LCP_PIVOT_PATHSEARCH:
    default:
      dim2 = 2 * (dim + 1);
  }

  int stack_size = 0;
  // with pathsearch we need a stack of the basis
  if (pivot_selection_rule == SICONOS_LCP_PIVOT_PATHSEARCH)
  {
    stack_size = options->iparam[SICONOS_IPARAM_PATHSEARCH_STACKSIZE];
    assert(stack_size >= 1);
  }
  if (preAlloc)
  {
    if (!options->iWork)
    {
      options->iWork = (int *)malloc(dim * (1 + stack_size) * sizeof(int));
    }
    if(!options->dWork)
      options->dWork = (double *)malloc((stack_size + dim * dim2 )*sizeof(double));

    basis = options->iWork;
    mat = &options->dWork[stack_size];
    t_stack = options->dWork;
  }
  else
  {
    // with pathsearch we need a stack of the basis
    if (pivot_selection_rule == SICONOS_LCP_PIVOT_PATHSEARCH)
    {
      stack_size = options->iparam[SICONOS_IPARAM_PATHSEARCH_STACKSIZE];
      assert(stack_size >= 1);
      basis = (int *)malloc(dim * stack_size * sizeof(int));
    }
    else
    {
      basis = (int *)malloc(dim * sizeof(int));
    }
    mat = (double *)malloc((stack_size + dim * dim2 )* sizeof(double));
    t_stack = &mat[dim * dim2];
  }

  assert((pivot_selection_rule != SICONOS_LCP_PIVOT_PATHSEARCH) ||
      ((pivot_selection_rule == SICONOS_LCP_PIVOT_PATHSEARCH) && t_stack != mat));

  assert(problem->q);

  switch (pivot_selection_rule)
  {
    case SICONOS_LCP_PIVOT_BARD:
      init_M_bard(mat, M, dim, problem->q);
      break;
    case SICONOS_LCP_PIVOT_LEAST_INDEX:
      init_M_least_index(mat, M, dim, problem->q);
      break;
    case SICONOS_LCP_PIVOT_PATHSEARCH:
      init_M_lemke_warm_start(dim, u, mat, M, problem->q, basis, cov_vec);
      basis_init = 1;
      DEBUG_EXPR_WE( DEBUG_PRINT("basis after hot start: ")
          for (unsigned int i = 0; i < dim; ++i)
          { DEBUG_PRINTF("%i ", basis[i])}
          DEBUG_PRINT("\n"));

      break;
    case SICONOS_LCP_PIVOT_LEMKE:
    default:
        init_M_lemke(mat, M, dim, dim, problem->q, cov_vec);
  }
  DEBUG_PRINT("lcp_pivot: init done, starting resolution\n");

  DEBUG_PRINT("total matrix\n");
  DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { for(unsigned int j = 0 ; j < dim2; ++j)
      { DEBUG_PRINTF("% 2.2e ", mat[i + j*dim]) }
      DEBUG_PRINT("\n")});
  /* End of construction of mat */

  /* Init basis if necessary */
  if (!basis_init)
  {
    for (unsigned int i = 0 ; i < dim  ; ++i) basis[i] = i + 1;
  }

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
    case SICONOS_LCP_PIVOT_PATHSEARCH:
      block = pivot_init_pathsearch(dim, mat, &t_indx);
      break;
    case SICONOS_LCP_PIVOT_LEMKE:
    default:
      block = pivot_init_lemke(mat, dim);
  }

  if (block < 0)
  {
    if (block == -1)
    {
      /** exit, the solution is at hand with the current basis */
      DEBUG_PRINT("Trivial solution\n");
      goto exit_lcp_pivot;
    }
    else if (block == PIVOT_PATHSEARCH_SUCCESS)
    {
      DEBUG_PRINTF("lcp_pivot :: path search successful ! t_indx = %d\n", t_indx);
      bck_drive = t_indx; /* XXX correct ? */
      t_stack[nb_iter%stack_size] = 1.0;
      double pivot = 1.0; /* force value of pivot to avoid numerical issues */
      for (unsigned int i = 0; i < dim; ++i) mat[i] -= mat[i + drive*dim]*pivot;
      *info = 0;
      goto exit_lcp_pivot;
    }
    else if (block == -LCP_PATHSEARCH_NON_ENTERING_T)
    {
      /* exit, t could not become basic */
      assert(pivot_selection_rule == SICONOS_LCP_PIVOT_PATHSEARCH);
      DEBUG_PRINT("lcp_pivot :: t could not become basic, exiting\n");
      *info = LCP_PATHSEARCH_NON_ENTERING_T;
      goto exit_lcp_pivot;
    }
  }
  /* Pivot < mu , drive >  or < drive, drive > */

  DEBUG_PRINTF("Pivoting %i and %i\n", block, drive);
  pivot = mat[block + drive*dim];

  if (fabs(pivot) < DBL_EPSILON)
  {
    *info = LCP_PIVOT_NUL;
    if (verbose > 0)
      printf("the pivot is nul, the algorithm cannot be used !\n");
    goto exit_lcp_pivot;
  }

  /* update matrix */
  switch (pivot_selection_rule)
  {
    case SICONOS_LCP_PIVOT_BARD:
      do_pivot_driftless(mat, dim, dim2, block, drive);
      break;
    case SICONOS_LCP_PIVOT_LEAST_INDEX:
      do_pivot(mat, dim, dim2, block, drive);
      break;
    case SICONOS_LCP_PIVOT_LEMKE:
    case SICONOS_LCP_PIVOT_PATHSEARCH:
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
    case SICONOS_LCP_PIVOT_PATHSEARCH:
      DEBUG_PRINTF("t value : %le\n", mat[t_indx]);
    case SICONOS_LCP_PIVOT_LEMKE:
    default:
      /** one basic u is leaving and mu enters the basis */
      leaving = basis[block];
      basis[block] = drive;
  }

  DEBUG_PRINTF("leaving variable = %d\n", leaving);

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
      case SICONOS_LCP_PIVOT_PATHSEARCH:
        if (leaving < dim + 1)
        {
          drive = leaving + dim + 1;
        }
        else if (leaving > dim + 1)
        {
          drive = leaving - (dim + 1);
        }
        else /* XXX oulalla */
        {
          assert(0 && "leaving variable is t");
          drive = dim + 1;
        }
        block = pivot_selection_pathsearch(mat, dim, drive, t_indx);
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

    DEBUG_PRINTF("Blocking variable: %d\tDriving variable: %d\n", block, drive);

    if(block < 0)
    {

      /* We stop here: it either mean that the algorithm stops here or that there
       * is an issue with the LCP */
      if (block == -1)
      {
        switch (pivot_selection_rule)
        {
          case SICONOS_LCP_PIVOT_LEMKE:
          case SICONOS_LCP_PIVOT_PATHSEARCH:
            *info = LCP_PIVOT_RAY_TERMINATION;
            DEBUG_PRINT("The pivot column is nonpositive ! We are on ray !\n"
                "It either means that the algorithm failed or that the LCP is infeasible\n"
                "Check the class of the M matrix to find out the meaning of this\n");
          default:
            bck_drive = drive < dim + 1 ? drive - 1 : drive - dim - 2;
        }
        break;
      }
      /* path search was successful, t = 1, we need to update the value of the
       * basic variable, but we are done here :) */
      else if (block == PIVOT_PATHSEARCH_SUCCESS)
      {
        DEBUG_PRINTF("lcp_pivot :: path search successful ! t_indx = %d\n", t_indx);
        basis[t_indx] = drive;
        t_stack[nb_iter%stack_size] = 1.0;
        double pivot = (mat[t_indx] - 1.0)/mat[t_indx + drive*dim];
        for (unsigned int i = 0; i < dim; ++i) mat[i] -= mat[i + drive*dim]*pivot;
        mat[t_indx] = pivot;
        *info = 0;
        break;
      }

    }

    DEBUG_PRINTF("driving variable %i \n", drive);
    if (basis[block] == dim + 1)
    {
      if (pivot_selection_rule != SICONOS_LCP_PIVOT_PATHSEARCH)
      {
        has_sol = 1;
      }
      else
      {
        DEBUG_PRINT("t variable leaving !\n");
        *info = LCP_PATHSEARCH_LEAVING_T;
        bck_drive = drive < dim + 1 ? drive - 1 : drive - dim - 2;
        options->dparam[2] = mat[t_indx];
        goto exit_lcp_pivot;
      }
    }


    /* Pivot < block , drive > */
    DEBUG_PRINTF("Pivoting %i and %i\n", block, drive);

    pivot = mat[block + drive*dim];
    if (fabs(pivot) < DBL_EPSILON)
    {
      *info = LCP_PIVOT_NUL;
      if (verbose > 0)
        printf("the pivot is nul, the algorithm cannot be used !\n");
      goto exit_lcp_pivot;
    }

    /* update matrix */
    switch (pivot_selection_rule)
    {
      case SICONOS_LCP_PIVOT_BARD:
        do_pivot_driftless(mat, dim, dim2, block, drive);
        break;
      case SICONOS_LCP_PIVOT_LEAST_INDEX:
        do_pivot(mat, dim, dim2, block, drive);
        break;
      case SICONOS_LCP_PIVOT_LEMKE:
      case SICONOS_LCP_PIVOT_PATHSEARCH:
      default:
        do_pivot_driftless(mat, dim, dim2, block, drive);
    }

    /* determine leaving variable and update basis */
    switch (pivot_selection_rule)
    {
      /* Principal Pivoting Methods  */
      case SICONOS_LCP_PIVOT_BARD:
        basis[block] = basis[block] <= dim ? block + dim + 2 : block + 1;
        break;
      case SICONOS_LCP_PIVOT_LEAST_INDEX:
        basis[block] = basis[block] <= dim ? block + dim + 2 : block + 1;
        break;
      case SICONOS_LCP_PIVOT_PATHSEARCH:
        leaving = basis[block];
        //memcpy(basis, basis+dim, dim*sizeof(int));
        //basis += dim;
        basis[block] = drive;
        t_stack[nb_iter%stack_size] = mat[t_indx];
        DEBUG_PRINTF("t value : %2.2e\n", mat[t_indx]);
        DEBUG_PRINTF("t-1.0 value : %2.2e\n", mat[t_indx]-1.0);
        /* XXX to test */
//        if (fabs(mat[t_indx] -1.0) < 1e-8)
//        {
//          double pivot = (mat[t_indx] - 1.0)/mat[t_indx + drive*dim];
//          for (unsigned int i = 0; i < dim; ++i) mat[i] -= mat[i + drive*dim]*pivot;
//          mat[t_indx] = 0;
//          *info = 0;
//          has_sol = 1;
//        }
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
      { DEBUG_PRINTF("% 2.2e ", mat[i+ j*dim]) }
      DEBUG_PRINT("\n")});
  } /* end while*/

exit_lcp_pivot:

  DEBUG_EXPR_WE( DEBUG_PRINT("final basis: ")
      for (unsigned int i = 0; i < dim; ++i)
      { DEBUG_PRINTF("%i ", basis[i])}
      DEBUG_PRINT("\n"));

  DEBUG_PRINT("total matrix\n");
  DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { for(unsigned int j = 0 ; j < dim2; ++j)
      { DEBUG_PRINTF("% 2.2e ", mat[i + j*dim]) }
      DEBUG_PRINT("\n")});

  /* Recover solution */
  for (unsigned int i = 0 ; i < dim; ++i)
  {
    drive = basis[i];
    assert(drive > 0);
    //assert(drive != dim + 1);
    if (drive < dim + 1)
    {
      u[drive - 1] = 0.0;
      s[drive - 1] = mat[i];
    }
    else if (drive > dim + 1)
    {
      u[drive - dim - 2] = mat[i];
      s[drive - dim - 2] = 0.0;
    }
    else
    {
      assert(bck_drive >= 0);
      u[bck_drive] = 0.0;
      s[bck_drive] = 0.0;
    }

  }

  DEBUG_PRINT("u s\n");
  DEBUG_EXPR_WE(for (unsigned int i = 0; i < dim; ++i)
      { DEBUG_PRINTF("%e %e\n", u[i], s[i]) });

  options->iparam[1] = nb_iter;

  /* update info */
  switch (pivot_selection_rule)
  {
    /* Principal Pivoting Methods  */
    case SICONOS_LCP_PIVOT_BARD:
    case SICONOS_LCP_PIVOT_LEAST_INDEX:
      *info = lcp_compute_error(problem, u, s, options->dparam[0], &tmp);
      break;
    case SICONOS_LCP_PIVOT_PATHSEARCH:
      break; /* info should already be set */
    case SICONOS_LCP_PIVOT_LEMKE:
    default:
      /* info may already be set*/
      if (*info == 0)
      {
        if (has_sol) *info = 0;
        else *info = 1;
      }
  }

  if (*info > 0)
  {
    DEBUG_PRINT("No solution found !\n");
  }

  if(!preAlloc)
  {
    free(basis);
    free(mat);
  }
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
  options->iSize = 6;
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
