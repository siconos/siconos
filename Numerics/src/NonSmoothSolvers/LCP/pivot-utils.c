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
#include <assert.h>
#include <stddef.h>
#include <string.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include "pivot-utils.h"
#include "lcp_pivot.h"

#include "lumod_wrapper.h"
#include "SiconosCompat.h"
#include "SiconosBlas.h"

#define DEBUG_STDOUT
#define DEBUG_MESSAGES
#include "debug.h"

int pivot_init_lemke(double* mat, unsigned int dim)
{
  int block = 0;
  double zb, dblock;
  double z0 = mat[0];

  for (unsigned int i = 1 ; i < dim ; ++i)
  {
    zb = mat[i];
    if (zb < z0)
    {
      z0 = zb;
      block = i;
    }
    else if (zb == z0)
    {
      for (unsigned int j = 1; j <= dim; ++j)
      {
        dblock = mat[block + j*dim] - mat[i + j*dim];
        if (dblock < 0.) break;
        else if (dblock > 0.)
        {
          block = i;
          break;
        }
      }
    }
  }
  /* XXX check that */
  return z0 < 0.0 ? block : -1 ;
}

int pivot_init_pathsearch(unsigned dim, double* mat, unsigned* t_indx)
{
  int block = -LCP_PATHSEARCH_NON_ENTERING_T;
  double* minus_r = &mat[dim*(dim+1)];
  double t_min =  INFINITY;
  for (unsigned i = 0; i < dim; ++i)
  {
    /* -r = mat[dim(dim+1) + ...] */
    if (fabs(minus_r[i]) > DBL_EPSILON) /* XXX tolerance */
    {
      double tt = mat[i]/minus_r[i];
      if (tt >= 0.0)
      {
        if (tt < t_min)
      {
        t_min = tt;
        block = i;
      }
      }
      else /* ratio neg, t could be increased at will */
      {
        if (block < 0) /* accept only if no variable has been detected as blocking */
        {
      DEBUG_PRINTF("pivot_init_pathsearch :: non-blocking variable %d\n", i);
          t_min = 1.0;
          block = i;
        }
      }
    }
  }

  *t_indx = block;

  if (block >= 0)
  {
//    if (t_min >= 1.0 - DBL_EPSILON*fabs(minus_r[block])) /* XXX need a tol here ?*/
    if (t_min >= 1.0)
    {
      DEBUG_PRINTF("pivot_init_pathsearch :: search successful, t = %.*e\n", DECIMAL_DIG, t_min);
      block = PIVOT_PATHSEARCH_SUCCESS;
    }
    else
    {
      DEBUG_PRINTF("pivot_init_pathsearch :: expected t value = %.*le; blocking variable indx = %d\n", DECIMAL_DIG, t_min, block);
    }
  }


  return block;
}

int pivot_selection_lemke(double* mat, unsigned dim, unsigned drive, unsigned aux_indx)
{
  int block = -1;
  double candidate_pivot, candidate_ratio, dblock;
  double ratio = INFINITY;
  for (unsigned i = 0 ; i < dim ; ++i)
  {
    candidate_pivot = mat[i + drive*dim];
    if (candidate_pivot > 0.)
    {
      candidate_ratio = mat[i] / candidate_pivot;
      if (candidate_ratio > ratio) continue;
      else if (candidate_ratio < ratio)
      {
        ratio = candidate_ratio;
        block = i;
      }
      else
      {
        if (block == aux_indx || i == aux_indx)
        {
          /* We want the auxilliary variable to exit before any othe.
           * see CPS p. 279 and example 4.4.16 */
          block = aux_indx;
        }
        else
        {
          double current_pivot = mat[block + drive*dim];
          for (unsigned j = 1; j <= dim; ++j)
          {
            assert(block >= 0 && "ratio_selection_lemke: block < 0");
            dblock = mat[block + j*dim] / current_pivot - mat[i + j*dim] / candidate_pivot;
            if (dblock < 0.) break;
            else if (dblock > 0.)
            {
              block = i;
              break;
            }
          }
        }
      }
    }
  }
  return block;
}

/* This version is with   */
int pivot_selection_lemke2(unsigned n, double* restrict col_drive, double* restrict q_tilde, double* restrict lexico_mat, unsigned drive, unsigned aux_indx)
{
  int block = -1;
  double candidate_pivot, current_pivot, candidate_ratio, dblock;
  double ratio = INFINITY;
  for (unsigned i = 0 ; i < n ; ++i)
  {
    candidate_pivot = col_drive[i];
    if (candidate_pivot > DBL_EPSILON)
    {
      candidate_ratio = q_tilde[i] / candidate_pivot;
      if (candidate_ratio > ratio) continue;
      else if (candidate_ratio < ratio)
      {
        ratio = candidate_ratio;
        block = i;
        current_pivot = col_drive[block];
      }
      else
      {
        if (block == aux_indx || i == aux_indx)
        {
          /* We want the auxilliary variable to exit before any othe.
           * see CPS p. 279 and example 4.4.16 */
          block = aux_indx;
        }
        else
        {
          for (unsigned j = 0; j < n; ++j)
          {
            assert(block >= 0 && "ratio_selection_lemke: block < 0");
            dblock = lexico_mat[block*n + j] / current_pivot - lexico_mat[i*n + j] / candidate_pivot;
            if (dblock < 0.) break;
            else if (dblock > 0.)
            {
              block = i;
              break;
            }
          }
        }
      }
    }
  }
  return block;
}


int pivot_selection_pathsearch(double* mat, unsigned dim, unsigned drive, unsigned t_indx)
{
  int block = pivot_selection_lemke(mat, dim, drive, t_indx);
  double pivot_t = mat[t_indx + drive*dim];
  if (pivot_t <= -DBL_EPSILON)
  {
    /* try to set t to 1 */
    double ratio_t = (mat[t_indx] - 1.0)/pivot_t;
    if (block >= 0)
    {
      double ratio_lemke = mat[block]/mat[block + drive*dim];
      if (ratio_t <= ratio_lemke)
      {
        block = PIVOT_PATHSEARCH_SUCCESS;
      }
    }
    else
    {
      block = PIVOT_PATHSEARCH_SUCCESS;
    }
  }

  return block;
}

void init_M_lemke(double* restrict mat, double* restrict M, unsigned int dim, unsigned int size_x, double* restrict q, double* restrict d)
{
  /* construction of mat matrix such that
   * mat = [ q | Id | -d | -M ] with d_i = 1 if i < size_x
   */

  /* We need to init only the part corresponding to Id */
  memset(&mat[dim], 0, sizeof(double) * dim * dim);

  /*  Copy M but mat[dim+2:, :] = -M */
  for (unsigned int i = 0 ; i < dim; ++i)
    for (unsigned int j = 0 ; j < dim; ++j)
      mat[i + dim*(j + dim + 2)] = -M[dim*j + i]; // Siconos is in column major


  for (unsigned int i = 0 ; i < dim; ++i)
  {
    mat[i] = q[i];
    mat[i + dim*(i + 1)] =  1.0;
  }

  /** Add covering vector */
  if (d != NULL)
    for (unsigned int i = 0; i < size_x  ; ++i) mat[i + dim*(dim + 1)] = d[i];
  else
    for (unsigned int i = 0; i < size_x  ; ++i) mat[i + dim*(dim + 1)] = -1.0;
  for (unsigned int i = size_x; i < dim; ++i) mat[i + dim*(dim + 1)] = 0.0;
}

void do_pivot_driftless(double* mat, unsigned int dim, unsigned int dim2, unsigned int block, unsigned int drive)
{
  if (fabs(mat[block + drive*dim]) < DBL_EPSILON)
  {
    printf("do_pivot_driftless :: pivot value too small %e; q[block] = %e; theta = %e\n", mat[block + drive*dim], mat[block], mat[block]/mat[block + drive*dim]);
  }
  double pivot_inv = 1.0/mat[block + drive*dim];
  unsigned ncols = dim*dim2;

  /* Update column mat[block, :] */
  mat[block + drive*dim] = 1.; /* nm_rs = 1 */
  /* nm_rj = m_rj/m_rs */
  for (unsigned int i = 0        ; i < drive; ++i) mat[block + i*dim] *= pivot_inv;
  for (unsigned int i = drive + 1; i < dim2 ; ++i) mat[block + i*dim] *= pivot_inv;

  /* Update other columns*/
  for (unsigned int i = 0; i < block; ++i)
  {
    double tmp = mat[i + drive*dim];
    /* nm_ij = m_ij + (m_ir/m_rs)m_rj = m_ij - m_is*nm_rj */
    for (unsigned int j = 0; j < ncols; j+=dim) mat[i + j] -= tmp*mat[block + j];
  }
  for (unsigned int i = block + 1; i < dim; ++i)
  {
    double tmp = mat[i + drive*dim];
    /* nm_ij = m_ij + (m_ir/m_rs)m_rj = m_ij - m_is*nm_rj */
    for (unsigned int j = 0; j < ncols; j+=dim) mat[i + j] -= tmp*mat[block + j];
  }
}

void do_pivot_driftless2(double* mat, unsigned int dim, unsigned int dim2, unsigned int block, unsigned int drive)
{
  double tmp;
  double pivot_inv = 1.0/mat[block + drive*dim];

  /* Update column mat[block, :] */
  mat[block + drive*dim] = pivot_inv; /* nm_rs = 1 */
  /*  nm_rj = m_rj/m_rs */
  for (unsigned int i = 0        ; i < drive; ++i) mat[block + i*dim] *= pivot_inv;
  for (unsigned int i = drive + 1; i < dim2 ; ++i) mat[block + i*dim] *= pivot_inv;

  /* Update other columns*/
  for (unsigned int i = 0; i < block; ++i)
  {
    tmp = mat[i + drive*dim];
    /* nm_ij = m_ij + (m_ir/m_rs)m_rj = m_ij - m_is*nm_rj */
    for (unsigned int j = 0        ; j < drive; ++j) mat[i + j*dim] -= tmp*mat[block + j*dim];
    for (unsigned int j = drive + 1; j < dim2 ; ++j) mat[i + j*dim] -= tmp*mat[block + j*dim];
  }
  for (unsigned int i = block + 1; i < dim; ++i)
  {
    tmp = mat[i + drive*dim];
    for (unsigned int j = 0        ; j < drive; ++j) mat[i + j*dim] -= tmp*mat[block + j*dim];
    for (unsigned int j = drive + 1; j < dim2 ; ++j) mat[i + j*dim] -= tmp*mat[block + j*dim];
  }
}

/* Standard pivot <block, drive>  */
void do_pivot(double* mat, unsigned int dim, unsigned int dim2, unsigned int block, unsigned int drive)
{
  double tmp;
  double pivot_inv = 1.0/mat[block + drive*dim];
  double mpivot_inv = -1.0/mat[block + drive*dim];

  /* Update column mat[block, :] */
  mat[block + drive*dim] = pivot_inv;   /* nm_rs = 1/m_rs  */
  /* nm_rj = -m_rj/m_rs */
  for (unsigned int i = 0        ; i < drive; ++i) mat[block + i*dim] *= mpivot_inv;
  for (unsigned int i = drive + 1; i < dim2 ; ++i) mat[block + i*dim] *= mpivot_inv;

  /* Update other lines*/
  for (unsigned int i = 0; i < block; ++i)
  {
    tmp = mat[i + drive*dim];
    /* nm_ij = m_ij - (m_is/m_rs)m_rj = m_ij + m_is*nm_rj */
    for (unsigned int j = 0        ; j < drive; ++j) mat[i + j*dim] += tmp*mat[block + j*dim];
    for (unsigned int j = drive + 1; j < dim2 ; ++j) mat[i + j*dim] += tmp*mat[block + j*dim];
    mat[i + drive*dim] *= pivot_inv; /* nm_is = m_is/m_rs */
  }
  for (unsigned int i = block + 1; i < dim; ++i)
  {
    tmp = mat[i + drive*dim];
    for (unsigned int j = 0        ; j < drive; ++j) mat[i + j*dim] += tmp*mat[block + j*dim];
    for (unsigned int j = drive + 1; j < dim2 ; ++j) mat[i + j*dim] += tmp*mat[block + j*dim];
    mat[i + drive*dim] *= pivot_inv;
  }
}

void do_pivot_lumod(SN_lumod_dense_data* restrict lumod_data, NumericsMatrix* restrict M, double* restrict q_tilde, double* restrict lexico_mat, double* restrict col_drive, double* restrict col_tilde, unsigned* basis, unsigned block, unsigned drive)
{
  unsigned n = lumod_data->n;

  /* Get */
  unsigned leaving_var_indx = basis[block]-1;
  unsigned entering_var_indx = drive-1;
  /* this variable is p+1 in QianLi */
  unsigned pos_leaving_var_in_factorized_basis = lumod_data->factorized_basis[leaving_var_indx];
  /* Determine the type of pivot operation */
  unsigned leaving_type = pos_leaving_var_in_factorized_basis > 0 ? 1 : 0;
  unsigned entering_type = lumod_data->factorized_basis[entering_var_indx] > 0 ? 2 : 0;

  switch (leaving_type + entering_type)
  {
    case 0:
      /* Both variables are not in the factorized basis
       * -> change the column corresponding to the leaving variable */
      SN_lumod_replace_col(lumod_data, -lumod_data->row_col_indx[leaving_var_indx], col_tilde);
      lumod_data->row_col_indx[entering_var_indx] = lumod_data->row_col_indx[leaving_var_indx];
      lumod_data->row_col_indx[leaving_var_indx] = 0;
      break;
    case 1:
      /* Leaving variable is in the factorized basis, entering not
       * -> add one row and col */
      /* Save the index of the column associated with the entering variable that
       * was not in the basis */
        lumod_data->row_col_indx[entering_var_indx] = -lumod_data->k; /* col index is negative */
        lumod_data->row_col_indx[leaving_var_indx] = lumod_data->k;
        SN_lumod_add_row_col(lumod_data, pos_leaving_var_in_factorized_basis-1, col_tilde);
      break;
    case 2:
      {
      /*  Leaving variable is not in the factorized basis, entering is.
       *  -> remove one row and one col */
        DEBUG_PRINT_VEC_INT_STR("old index", lumod_data->row_col_indx, 2*n+1);
      int index_row = lumod_data->row_col_indx[entering_var_indx];
      int index_col = lumod_data->row_col_indx[leaving_var_indx];
      assert(index_row >= 0);
      assert(index_col <= 0);
      SN_lumod_delete_row_col(lumod_data, index_row, -index_col);
      /*  update the index, since a row ans col was deleted */
      unsigned changed_row = SN_lumod_find_arg_var(lumod_data->row_col_indx, lumod_data->k, n);
      unsigned changed_col = SN_lumod_find_arg_var(lumod_data->row_col_indx, -lumod_data->k, n);
      DEBUG_PRINTF("Changing row for variable %d and col for variable %d\n", changed_row, changed_col);
      lumod_data->row_col_indx[changed_row] = index_row;
      lumod_data->row_col_indx[changed_col] = index_col;
      lumod_data->row_col_indx[entering_var_indx] = 0;
      lumod_data->row_col_indx[leaving_var_indx] = 0;
      DEBUG_PRINT_VEC_INT_STR("new index", lumod_data->row_col_indx, 2*n+1);
/*       if ((index_row <= lumod_data->k) || (index_col >= -lumod_data->k))
      {
        for (unsigned i = 0; i < 2*n + 1; ++i)
        {
          if (lumod_data->row_col_indx[i] > index_row)
          {
            --(lumod_data->row_col_indx[i]);
          }
          else if(lumod_data->row_col_indx[i] < index_col)
          {
            ++(lumod_data->row_col_indx[i]);
          }
        }
      }*/
        DEBUG_PRINT_VEC_INT(lumod_data->row_col_indx, 2*n+1);
      }
      break;
    case 3:
      /* Both variables are in the factorized basis 
       * -> replace the row corresponding to the entering variable by the
       *  leaving one*/
      SN_lumod_replace_row(lumod_data, lumod_data->row_col_indx[entering_var_indx], pos_leaving_var_in_factorized_basis-1);
      lumod_data->row_col_indx[leaving_var_indx] = lumod_data->row_col_indx[entering_var_indx];
      lumod_data->row_col_indx[entering_var_indx] = 0;
      break;
    default:
      printf("do_pivot_lumod :: unkown update type occuring\n");
      exit(EXIT_FAILURE);
  }

  /* Update q_tilde */
  double theta = q_tilde[block]/col_drive[block];
  DEBUG_PRINTF("theta = %e\n", theta);
  DEBUG_PRINT_VEC(col_drive, n);
  cblas_daxpy(n, -theta, col_drive, 1, q_tilde, 1);
  q_tilde[block] = theta;

  /* Update the lexico_mat
   * XXX check if this is correct. The value of the pivot may be wrong --xhub */
  double pivot = col_drive[block];
  unsigned block_row_indx = block*n;

  cblas_dscal(n, 1./pivot, &lexico_mat[block_row_indx], 1);
  for (unsigned i = 0; i < n*n; i += n)
  {
    if (i == block_row_indx) continue;
    cblas_daxpy(n, -col_drive[i]/pivot, &lexico_mat[block_row_indx], 1, &lexico_mat[i], 1);
  }
}

void lcp_pivot_diagnose_info(int info)
{
  switch (info)
  {
    case 0:
      printf("lcp_pivot :: a solution to the problem has been found\n");
      break;
    case 1:
      printf("lcp_pivot :: no solution have been found\n");
      break;
    case LCP_PIVOT_NUL:
      printf("lcp_pivot :: the pivot is nul\n");
      break;
    case LCP_PATHSEARCH_LEAVING_T:
      printf("lcp_pivot :: the pathsearch was not successful since the t variable was leaving\n");
      break;
    case LCP_PATHSEARCH_NON_ENTERING_T:
      printf("lcp_pivot :: t could not become basic");
      break;
    default:
      printf("lcp_pivot :: the info code %d is not documented\n", info);
  }
}
