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


#include <math.h>
#include <float.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "LinearComplementarityProblem.h"
#include "LCP_Solvers.h"
#include "lcp_cst.h"
#include "SolverOptions.h"
#include "NumericsMatrix.h"

#include "pivot-utils.h"
#include "lumod_wrapper.h"
#include "numerics_verbose.h"
#include "SiconosLapack.h"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
//#define DEBUG_NO_MATRIX
#include "debug.h"

#define NO_LEXICO_MAT

#define WARN_ONLY_SMALL_PIVOT
#include "lcp_pivot.h"
#include "pivot-utils.h"
#define LEXICO_TOL 1e3*DBL_EPSILON

DEBUG_GLOBAL_VAR_DECL(unsigned * basis_global;);

inline static double* get_q_tilde(double* mat, unsigned n)
{
  return mat;
}

inline static double* get_driving_col(double* mat, unsigned n)
{
  return &mat[n];
}

inline static double* get_lexico_mat(double* mat, unsigned n)
{
  return &mat[2*n];
}

inline static double* get_col_tilde(double* mat, unsigned n)
{
  return &mat[(n+2)*n];
}

inline static double* get_cov_vec(double* mat, unsigned n)
{
  return &mat[(n+3)*n];
}

void lcp_pivot_lumod(LinearComplementarityProblem* problem, double* u , double* s, int *info , SolverOptions* options)
{
  lcp_pivot_lumod_covering_vector(problem, u, s, info, options, NULL);
}


void lcp_pivot_lumod_covering_vector(LinearComplementarityProblem* problem, double* restrict u , double* restrict s, int *info , SolverOptions* options, double* restrict cov_vec)
{
  /* matrix M of the LCP */
  assert(problem);
  assert(problem->M);
  double* M = problem->M->matrix0;
  assert(M);
  assert(problem->q);


  unsigned int dim = problem->size;
  assert(dim>0);
  /* unsigned int dim2; */
  /* size of the LCP */
  DEBUG_EXPR_WE( DEBUG_PRINT("matrix M: ") NM_display(problem->M); DEBUG_PRINT("vector q: ")
      for(unsigned i = 0; i < dim; ++i) {printf("%e ", problem->q[i]);} printf("\n");
      if (cov_vec) { DEBUG_PRINT("covering vector: ") for(unsigned i = 0; i < dim; ++i) {printf("%e ", cov_vec[i]);}printf("\n");});

  unsigned drive = dim+1;
  int bck_drive = -1;
  int block = -1;
  unsigned has_sol = 0;
  unsigned nb_iter = 0;
  unsigned leaving = 0;
  unsigned itermax = options->iparam[0];
  unsigned preAlloc = options->iparam[SICONOS_IPARAM_PREALLOC];
  unsigned pivot_selection_rule = options->iparam[SICONOS_IPARAM_PIVOT_RULE];

  assert(itermax > 0 && "lcp_pivot_lumod_covering_vector itermax == 0, the algorithm will not run");
  double pivot;
  double tmp;
  unsigned* basis = (unsigned*) malloc(dim*sizeof(unsigned));
  DEBUG_EXPR_WE(basis_global = basis;);
  unsigned* candidate_indx = (unsigned*) malloc(dim*sizeof(unsigned));
  int basis_init = 0; /* 0 if basis was not initialized, 1 otherwise*/
  unsigned t_indx = 0;
  unsigned aux_indx = 0;
#if 0
  double* t_stack = NULL;
#endif
  /* This matrix contains q, the solution to the linear system Hk x = driving_col,
   * the matrix for the lexicographic ordering and the solution to the linear
   * system H x = driving_col. */
  double* mat = (double*) calloc((dim+4)*dim, sizeof(double)); /* XXX memory save */
  assert(problem->q);
  cblas_dcopy(dim, problem->q, 1, get_q_tilde(mat, dim), 1);

  if (cov_vec)
  {
    cblas_dcopy(dim, cov_vec, 1, get_cov_vec(mat, dim), 1);
  }
  else
  {
    double* d = get_cov_vec(mat, dim);
    for (unsigned i = 0; i < dim; ++i) d[i] = 1.;
  }

  /* Init the lexicographic mat */
  double* lexico_mat = get_lexico_mat(mat, dim);
  for (unsigned i = 0; i < dim*dim; i += dim+1) lexico_mat[i] = 1.;
  DEBUG_PRINT_MAT_ROW_MAJOR_NCOLS_SMALL_STR("lexico_mat", lexico_mat, dim, dim, dim);

  /* Maximum number of columns changed in the matrix */
  /* TODO: user settable and should not be bigger than the size of the matrix? */
  unsigned maxmod = 50;

  *info = 0;

  /*output*/
  options->iparam[1] = 0;

  /* Allocation */
  SN_lumod_dense_data* lumod_data = SN_lumod_dense_allocate(dim, maxmod);

/*   switch (pivot_selection_rule) */
/*   { */
/* /\*     case SICONOS_LCP_PIVOT_BARD: */
/*       dim2 = 2*dim + 1; */
/*       break; */
/*     case SICONOS_LCP_PIVOT_LEAST_INDEX: */
/*       dim2 = dim + 1; */
/*       break;*\/ */
/*     case SICONOS_LCP_PIVOT_LEMKE: */
/*     case SICONOS_LCP_PIVOT_PATHSEARCH: */
/*     default: */
/*       dim2 = 2 * (dim + 1); */
/*   } */


  /* Init basis if necessary */
  if (!basis_init)
  {
    for (unsigned i = 0 ; i < dim ; ++i) basis[i] = i + 1;
  }

  /* Looking for pivot */
  switch (pivot_selection_rule)
  {
/*     case SICONOS_LCP_PIVOT_BARD:
      block = pivot_selection_bard(mat, dim);
      drive = block + dim + 1;
      break;
    case SICONOS_LCP_PIVOT_LEAST_INDEX:
      block = pivot_selection_least_index(mat, dim);
      drive = block + 1;
      break;
    case SICONOS_LCP_PIVOT_PATHSEARCH:
      block = pivot_init_pathsearch(dim, mat, &t_indx);
      break;*/
    case SICONOS_LCP_PIVOT_LEMKE:
    default:
//      block = pivot_init_lemke(get_q_tilde(mat, dim), dim);
      block = pivot_selection_lemke2(dim, get_cov_vec(mat, dim), get_q_tilde(mat, dim), get_lexico_mat(mat, dim), INT_MAX, LEXICO_TOL);
  }

  if (block < 0)
  {
    if (block == -1)
    {
      /** exit, the solution is at hand with the current basis */
      DEBUG_PRINT("Trivial solution\n");
      goto exit_lcp_pivot;
    }
#if 0
    else if (block == PIVOT_PATHSEARCH_SUCCESS)
    {
      assert(pivot_selection_rule == SICONOS_LCP_PIVOT_PATHSEARCH);
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
#endif
  }

  /* save the position of the auxiliary variable */
  assert(block >= 0);
  aux_indx = block;

  /* Update the basis */
  switch (pivot_selection_rule)
  {
#if 0
    /* Principal Pivoting Methods  */
    case SICONOS_LCP_PIVOT_BARD:
      basis[block] = basis[block] <= dim ? block + dim + 2 : block + 1;
      break;
    case SICONOS_LCP_PIVOT_LEAST_INDEX:
      basis[block] = basis[block] <= dim ? block + dim + 2 : block + 1;
      break;
    case SICONOS_LCP_PIVOT_PATHSEARCH:
      DEBUG_PRINTF("t value : %le\n", mat[t_indx]);
#endif
    case SICONOS_LCP_PIVOT_LEMKE:
    default:
      /** one basic u is leaving and mu enters the basis */
      leaving = basis[block];
      basis[block] = drive;
  }

  /* Init the LUMOD data and perform the pivot < aux_variable , drive > */
  switch (pivot_selection_rule)
  {
/*     case SICONOS_LCP_PIVOT_BARD:
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

      break;*/
    case SICONOS_LCP_PIVOT_LEMKE:
    default:
      SN_lumod_factorize(lumod_data, basis, problem->M, get_cov_vec(mat, dim));
      pivot = get_cov_vec(mat, dim)[block];
  }
  DEBUG_PRINT("lcp_pivot: init done, starting resolution\n");

  /* XXX Maybe we should compute theta = q_i/pivot */
  if (fabs(pivot) < DBL_EPSILON)
  {
    if (verbose > 0)
      printf("the pivot is quasi-nul %e, the algorithm cannot be used !\n", pivot);
#ifndef WARN_ONLY_SMALL_PIVOT
    *info = LCP_PIVOT_NUL;
    goto exit_lcp_pivot;
#endif
  }

  /* Update q 
   * XXX maybe this code should go to pivot-utils*/
  {
    double* q = get_q_tilde(mat, dim);
    double theta = -q[block]/pivot;
    cblas_daxpy(dim, theta, get_cov_vec(mat, dim), 1, q, 1);
    q[block] = theta;

    unsigned block_row_indx = block*dim;
    for (unsigned i = 0, j = 0; i < dim; ++i, j += dim)
    {
      if (j == block_row_indx) continue;
      cblas_daxpy(dim, -get_cov_vec(mat, dim)[i]/pivot, &lexico_mat[block_row_indx], 1, &lexico_mat[j], 1);
    }
    cblas_dscal(dim, -1./pivot, &lexico_mat[block_row_indx], 1);
    DEBUG_PRINT_MAT_ROW_MAJOR_NCOLS_SMALL2_STR("lexico_mat", lexico_mat, dim, dim, dim, get_cov_vec(mat, dim));
  }
  DEBUG_PRINT_VEC(get_q_tilde(mat, dim), dim);


  DEBUG_EXPR_WE( DEBUG_PRINT("new basis: ")
      for (unsigned int i = 0; i < dim; ++i)
      { DEBUG_PRINTF("%i ", basis[i])}
      DEBUG_PRINT("\n"));

  while (nb_iter < itermax && !has_sol)
  {

    ++nb_iter;
    /*  Prepare the search for leaving variable */
    double* driving_col = get_driving_col(mat, dim);
    if (leaving < dim + BASIS_OFFSET) /* the leaving variable is w_i -> the driving variable is z_i */
    {
      drive = leaving + dim + BASIS_OFFSET;
      cblas_dcopy(dim, &M[dim*(leaving-BASIS_OFFSET)], 1, driving_col, 1);
    }
    else if (leaving > dim + BASIS_OFFSET) /*  the leaving variable is z_i -> the driving variable is w_i */
    {
      drive = leaving - (dim + BASIS_OFFSET);
      memset(driving_col, 0, sizeof(double) * dim);
      driving_col[drive - BASIS_OFFSET] = -1.;
    }
    else
    {
      printf("lcp_pivot_lumod the leaving variable is the auxiliary variable; we should not execute those lines!\n");
      exit(EXIT_FAILURE);
    }
    DEBUG_EXPR_WE( DEBUG_PRINT("basis= "); for (unsigned i = 0; i < dim; ++i) { DEBUG_PRINTF("%s%d ", basis_to_name(basis[i], dim), basis_to_number(basis[i], dim)); } DEBUG_PRINT("\n"));
    int solve_info = SN_lumod_dense_solve(lumod_data, driving_col, get_col_tilde(mat, dim));
    if (SN_lumod_need_refactorization(solve_info))
    {
      DEBUG_PRINT("Refactorizing!\n");
      SN_lumod_factorize(lumod_data, basis, problem->M, get_cov_vec(mat, dim));
      if (leaving < dim + BASIS_OFFSET) /* the leaving variable is w_i -> the driving variable is z_i */
      {
        drive = leaving + dim + BASIS_OFFSET;
        cblas_dcopy(dim, &M[dim*(leaving-BASIS_OFFSET)], 1, driving_col, 1);
      }
      else if (leaving > dim + BASIS_OFFSET) /*  the leaving variable is z_i -> the driving variable is w_i */
      {
        drive = leaving - (dim + BASIS_OFFSET);
        memset(driving_col, 0, sizeof(double) * dim);
        assert(drive >= BASIS_OFFSET);
        driving_col[drive - BASIS_OFFSET] = -1.;
      }
      solve_info = SN_lumod_dense_solve(lumod_data, driving_col, get_col_tilde(mat, dim));
    }
    if (solve_info != 0)
    {
      printf("lcp_pivot_lumod :: SN_lumod_dense_solve failed!, info = %d", solve_info);
      *info = LCP_PIVOT_LUMOD_FAILED;
      goto exit_lcp_pivot;
    }


    /* Start research of argmin lexico for minimum ratio test */

    /* Looking for pivot */
    switch (pivot_selection_rule)
    {
/*       case SICONOS_LCP_PIVOT_BARD:
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
        else // XXX oulalla
        {
          assert(0 && "leaving variable is t");
          drive = dim + 1;
        }
        block = pivot_selection_pathsearch(mat, dim, drive, t_indx);
        break;*/
      case SICONOS_LCP_PIVOT_LEMKE:
      default:
#ifndef NO_LEXICO_MAT
        block = pivot_selection_lemke2(dim, get_driving_col(mat, dim), get_q_tilde(mat, dim), get_lexico_mat(mat, dim), aux_indx, LEXICO_TOL);
#else
        block = pivot_selection_lemke3(dim, get_driving_col(mat, dim), get_q_tilde(mat, dim), get_lexico_mat(mat, dim), basis, candidate_indx, lumod_data, aux_indx, LEXICO_TOL);
#endif
    }

    DEBUG_PRINTF("leaving variable %s%d entering variable %s%d\n", basis_to_name(basis[block], dim), basis_to_number(basis[block], dim), basis_to_name(drive, dim), basis_to_number(drive, dim));

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
#if 0
      /* path search was successful, t = 1, we need to update the value of the
       * basic variable, but we are done here :) */
      else if (block == PIVOT_PATHSEARCH_SUCCESS)
      {
        assert(pivot_selection_rule == SICONOS_LCP_PIVOT_PATHSEARCH);
        DEBUG_PRINTF("lcp_pivot :: path search successful ! t_indx = %d\n", t_indx);
        basis[t_indx] = drive;
        t_stack[nb_iter%stack_size] = 1.0;
        double pivot = (mat[t_indx] - 1.0)/mat[t_indx + drive*dim];
        for (unsigned int i = 0; i < dim; ++i) mat[i] -= mat[i + drive*dim]*pivot;
        mat[t_indx] = pivot;
        *info = 0;
        break;
      }
#endif

    }

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
    DEBUG_PRINTF("Pivoting variable at pos %d in basis (%s%d) and (%s%d)\n", block, basis_to_name(basis[block], dim), basis_to_number(basis[block], dim), basis_to_name(drive, dim), basis_to_number(drive, dim));

    pivot = get_driving_col(mat, dim)[block];
    if (fabs(pivot) < DBL_EPSILON)
    {
      if (verbose > 0)
        printf("the pivot is quasi-nul %e, danger !\nq[block] = %e; z = %e\n", pivot, get_q_tilde(mat, dim)[block], get_q_tilde(mat, dim)[block]/pivot);
#ifndef WARN_ONLY_SMALL_PIVOT
      *info = LCP_PIVOT_NUL;
      goto exit_lcp_pivot;
#endif
    }

    /* update matrix and q*/
    switch (pivot_selection_rule)
    {
      /*    case SICONOS_LCP_PIVOT_BARD:
            do_pivot_driftless(mat, dim, dim2, block, drive);
            break;
            case SICONOS_LCP_PIVOT_LEAST_INDEX:
            do_pivot(mat, dim, dim2, block, drive);
            break;*/
      case SICONOS_LCP_PIVOT_LEMKE:
      case SICONOS_LCP_PIVOT_PATHSEARCH:
      default:
        do_pivot_lumod(lumod_data, problem->M, get_q_tilde(mat, dim), get_lexico_mat(mat, dim), get_driving_col(mat, dim), get_col_tilde(mat, dim), basis, block, drive);
    }
    DEBUG_PRINT_VEC(get_q_tilde(mat, dim), dim);

    /* determine leaving variable and update basis */
    switch (pivot_selection_rule)
    {
#if 0
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
#endif        /* XXX to test */
//        if (fabs(mat[t_indx] -1.0) < 1e-8)
//        {
//          double pivot = (mat[t_indx] - 1.0)/mat[t_indx + drive*dim];
//          for (unsigned int i = 0; i < dim; ++i) mat[i] -= mat[i + drive*dim]*pivot;
//          mat[t_indx] = 0;
//          *info = 0;
//          has_sol = 1;
//        }
//        break;
      case SICONOS_LCP_PIVOT_LEMKE:
      default:
        /** one basic variable is leaving and the driving one enters the basis */
        leaving = basis[block];
        basis[block] = drive;
    }

    DEBUG_PRINT_VEC_STR("basis value", get_q_tilde(mat, dim), dim);

    DEBUG_EXPR_WE( DEBUG_PRINT("new basis: ")
      for (unsigned int i = 0; i < dim; ++i)
      { DEBUG_PRINTF("%i ", basis[i])}
      DEBUG_PRINT("\n"));

  } /* end while*/

exit_lcp_pivot:

  DEBUG_EXPR_WE( DEBUG_PRINT("final basis: ")
      for (unsigned int i = 0; i < dim; ++i)
      { DEBUG_PRINTF("%i ", basis[i])}
      DEBUG_PRINT("\n"));

  /* Recover solution */
  double* finalq = get_q_tilde(mat, dim);
  for (unsigned int i = 0 ; i < dim; ++i)
  {
    drive = basis[i];
    assert(drive > 0);
    //assert(drive != dim + 1);
    if (drive < dim + 1)
    {
      u[drive - 1] = 0.0;
      s[drive - 1] = finalq[i];
    }
    else if (drive > dim + 1)
    {
      u[drive - dim - 2] = finalq[i];
      s[drive - dim - 2] = 0.0;
    }
    else
    {
      if (nb_iter < itermax)
      {
        assert(bck_drive >= 0);
        u[bck_drive] = 0.0;
        s[bck_drive] = 0.0;
      }
    }

  }
  /* End recover solution  */

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
    SM_lumod_dense_free(lumod_data);
  }

  /*  XXX Clean that --xhub */
  free(candidate_indx);
}


/*
int linearComplementarity_pivot_lumod_setDefaultSolverOptions(SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the generic pivot Solver\n");
  }

  solver_options_set(options, SICONOS_LCP_PIVOT_LUMOD);
  return 0;
}
*/
