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
#include "SiconosLapack.h"

#include "sanitizer.h"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"
#include "pivot-utils.h"
#define TOL_LEXICO DBL_EPSILON*10000
#define MIN_INCREASE 10
#define SMALL_PIVOT 1e-10
#define LEXICO_TOL_DISPLAY 1e-10
#define SIZE_CANDIDATE_PIVOT 5

#define BASIS_OFFSET 1

char* basis_to_name(unsigned nb, unsigned n)
{
  if (nb < n + BASIS_OFFSET) return "w";
  else if (nb > n + BASIS_OFFSET) return "z"; 
  else return "e";
}

unsigned basis_to_number(unsigned nb, unsigned n)
{
  if (nb < n + BASIS_OFFSET) return nb + 1 - BASIS_OFFSET;
  else if (nb > n + BASIS_OFFSET) return nb - n - BASIS_OFFSET;
  else return 0;
}


DEBUG_GLOBAL_VAR_DECL(extern unsigned * basis_global;)
DEBUG_GLOBAL_VAR_DECL(static unsigned max_pivot_helped;)

int lexicosort_lumod(unsigned* vars, unsigned n, unsigned nb_candidates, SN_lumod_dense_data* restrict lumod_data, unsigned* restrict basis, double* restrict lexico_col, double* restrict driving_col, double lexico_tol);


int lexicosort_lumod(unsigned* vars, unsigned n, unsigned nb_candidates, SN_lumod_dense_data* restrict lumod_data, unsigned* restrict basis, double* restrict lexico_col, double* restrict driving_col, double lexico_tol)
{
  /* strategy used here:
   * - compare all candidate on each column of -B^{-1} (because we take it to
   *   be a submatrix of [-I,M,d] and not [I,-M,-d] like in CPS), the inverse of the
   *   matrix basis
   * - B^{-1} has the following structure: if the ith variable in the basis is
   *   w_i, then the column is the same as in the identity matrix, otherwise,
   *   we have to invert -e_i to get vector. Let's have fun
   */
  DEBUG_PRINTF("lexicosort_lumod :: starting lexicosort with %d candidates\n", nb_candidates);
  DEBUG_EXPR_WE( DEBUG_PRINT("lexicosort_lumod :: candidates: ");
      for (unsigned i = 0; i < nb_candidates; ++i)
      { DEBUG_PRINTF("%s%d ", basis_to_name(basis[vars[i]], n), basis_to_number(basis[vars[i]], n)); }
          DEBUG_PRINT("\n") );
  assert(nb_candidates > 1);
  unsigned nb_remaining_vars = nb_candidates;
  for (unsigned i = 0; i < n; ++i)
  {
    unsigned staying_var_indx = 0;
    unsigned var_indx = 0;
    if (basis[i] == i + BASIS_OFFSET) continue;
    else
    {
      while ((vars[var_indx] < i) && var_indx < nb_remaining_vars) ++var_indx;
      /* We are lucky the lexicomin search was trivial -> the last variable is
       * the arg lexicomin -> return it as the blocking variable*/
      if (var_indx > nb_remaining_vars-1)
      {
        return vars[nb_remaining_vars-1];
      }
      vars[0] = vars[var_indx];
      memset(lexico_col, 0, n*sizeof(double));
      /* Remember that our matrix is the opposite of the one usually considered, hence the -1 */
      lexico_col[i] = -1.;
      /* If the factorized basis was the initial one, we just have to correct
       * for the terms in Ck, TODO  */
      SN_lumod_dense_solve(lumod_data, lexico_col, NULL);
      double current_pivot = driving_col[vars[var_indx]];
      double current_lexico_col_element = lexico_col[vars[var_indx]];
      unsigned zero_current_lexico_col_element = current_lexico_col_element == 0. ? 1 : 0;
      ++var_indx;
      for (unsigned next_var = vars[var_indx]; var_indx < nb_remaining_vars; next_var = vars[++var_indx])
      {
        double next_lexico_col_element = lexico_col[next_var];
        if (zero_current_lexico_col_element)
        {
          if (next_lexico_col_element == 0.) /* equality  */
          {
            vars[++staying_var_indx] = next_var;
            DEBUG_PRINTF("lexicosort_lumod :: adding var %s%d because of equality (both row elements are 0)\n", basis_to_name(basis[next_var], n), basis_to_number(basis[next_var], n));
          }
          else if (next_lexico_col_element < -lexico_tol) /* delta_lexico > 0 (since pivot are always >0., => new lexicomin  */
          {
            zero_current_lexico_col_element = 0;
            goto update_current_blocking_variable;
          }
          /* Implicit else -> delta_lexico < 0, next_var gets discarded  */
        }
        /* delta_lexico > 0 (since pivot are always >0., => new lexicomin  */
        else if ((next_lexico_col_element == 0.) && (current_lexico_col_element > lexico_tol))
        {
          zero_current_lexico_col_element = 1;
          goto update_current_blocking_variable;
        }
        else /* we need to actually compute delta_block, both element of the lexico_col are non-zero */
        {
          double candidate_pivot = driving_col[next_var];
          double delta_lexico  = current_lexico_col_element * candidate_pivot - next_lexico_col_element * current_pivot;
          if (delta_lexico > lexico_tol) /* Difference is significant */
          {
            DEBUG_PRINTF("lexicosort_lumod :: lexicomin change var block changes to %s%d from %s%d, delta_lexico = %2.2e, new pivot = %e\n",
                basis_to_name(basis[next_var], n), basis_to_number(basis[next_var], n), basis_to_name(basis[vars[staying_var_indx]], n),
                basis_to_number(basis[vars[staying_var_indx]], n), delta_lexico, candidate_pivot);
            zero_current_lexico_col_element = 0;
            goto update_current_blocking_variable;
          }
          else if (delta_lexico < - lexico_tol)
          {
            continue; /* Delete variable  */
          }
          else if (delta_lexico != 0.) /* XXX This is a hack */
          {
            if ((current_pivot < SMALL_PIVOT) && (candidate_pivot > MIN_INCREASE*current_pivot))
            {
              DEBUG_PRINTF("lexicosort_lumod :: lexicomin small difference %2.2e, taking largest pivot %e > %e (var %s%d vs %s%d)\n",
                  delta_lexico, candidate_pivot, current_pivot, basis_to_name(basis[next_var], n), basis_to_number(basis[next_var], n),
                  basis_to_name(basis[vars[staying_var_indx]], n), basis_to_number(basis[vars[staying_var_indx]], n));
              DEBUG_EXPR_WE(max_pivot_helped = 1;);
              zero_current_lexico_col_element = 0;
              goto update_current_blocking_variable;
            }
            else
            {
              vars[++staying_var_indx] = next_var;
              DEBUG_PRINTF("lexicosort_lumod :: adding var %s%d because of inconclusive test: delta_lexico = %2.2e\n", basis_to_name(basis[next_var], n), basis_to_number(basis[next_var], n), delta_lexico);
            }
          }
        }
        continue;
update_current_blocking_variable:
        DEBUG_PRINTF("lexicosort_lumod :: lexicomin blocking variable changes to %s%d from %s%d\n", basis_to_name(basis[next_var], n), basis_to_number(basis[next_var], n),
            basis_to_name(basis[vars[0]], n), basis_to_number(basis[vars[0]], n));
        vars[0] = next_var;
        staying_var_indx = 0;
        current_pivot = driving_col[next_var];
        current_lexico_col_element = next_lexico_col_element;
      }
    }
    if (staying_var_indx == 0)
    {
//      assert(0 && "Do I stop here ?");
      return vars[0];
    }
//    if (staying_var_indx == 1) /* We are done here  */
//    {
//      return vars[0];
//    }
    else
    {
      nb_remaining_vars = staying_var_indx + 1;
    }
  }
  assert(0 && " We should not be here, never");
  return 0;
}

static inline int lexicosort_rowmajor(unsigned var, int block, unsigned n, double* lexico_mat, double candidate_pivot, double current_pivot, double lexico_tol)
{
  assert(lexico_tol >= 0.);
  DEBUG_EXPR_WE(if (!basis_global) {DEBUG_PRINT("basis_global not initialized!");});
  double* lexico_row_block = &lexico_mat[block*n];
  double* lexico_row_var = &lexico_mat[var*n];
  for (unsigned j = 0; j < n; ++j)
  {
    assert(block >= 0 && "pivot_selection_lemke2: block < 0");
    double dblock = lexico_row_block[j] * candidate_pivot - lexico_row_var[j] * current_pivot;
    DEBUG_EXPR_WE(if ((dblock != 0.) && (fabs(dblock) < LEXICO_TOL_DISPLAY)) { printf("lexicosort_rowmajor :: very small difference in lexicomin search: %2.2e\n", dblock);
        unsigned block_number = basis_to_number(basis_global[block], n); unsigned var_number =  basis_to_number(basis_global[var], n);
        char* block_name = basis_to_name(basis_global[block], n); char* var_name = basis_to_name(basis_global[var], n);
        printf("lexicomin: A[%s%d][j] / A[%s%d][drive] = %e / %e vs A[%s%d][j] / A[%s%d][drive] = %e / %e\n",
            block_name, block_number, block_name, block_number, lexico_row_block[j], current_pivot, var_name,
            var_number, var_name, var_number, lexico_row_var[j], candidate_pivot);});
    if (dblock < -lexico_tol) return block; /* Unchanged blocking variable */
    else if (dblock > lexico_tol) /* Difference is significant */
    {
      DEBUG_PRINTF("pivot_selection_lemke :: lexicomin change var block changes to %s%d from %s%d, dblock = %2.2e, new pivot = %e\n",
          basis_to_name(basis_global[var], n), basis_to_number(basis_global[var], n), basis_to_name(basis_global[block], n),
          basis_to_number(basis_global[block], n), dblock, candidate_pivot);
      return var;
    }
    else if (dblock != 0.)
    {
      if ((current_pivot < SMALL_PIVOT) && (candidate_pivot > MIN_INCREASE*current_pivot))
      {
        DEBUG_PRINTF("pivot_selection_lemke :: lexicomin small difference %2.2e, taking largest pivot %e > %e (var %s%d vs %s%d)\n",
            dblock, candidate_pivot, current_pivot, basis_to_name(basis_global[var], n), basis_to_number(basis_global[var], n),
            basis_to_name(basis_global[block], n), basis_to_number(basis_global[block], n));
        DEBUG_EXPR_WE(max_pivot_helped = 1;);
        return var;
      }
    }
  }
  DEBUG_PRINTF("lexicosort_rowmajor :: reaching end of function, which means that the lexicomin failed!, block = %d, var = %d\n",
      block, var);
  return block; /* This should almost never happen  */
}

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
        if (block == (int)aux_indx || i == aux_indx)
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
            dblock = mat[block + j*dim] * candidate_pivot - mat[i + j*dim] * current_pivot;
            DEBUG_EXPR_WE(if (fabs(dblock) < TOL_LEXICO) printf("pivot_selection_lemke :: very small difference in lexicomin: %2.2e\n", dblock););
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
int pivot_selection_lemke2(unsigned n, double* restrict col_drive, double* restrict q_tilde, double* restrict lexico_mat, unsigned aux_indx, double lexico_tol)
{
  int block = -1;
  unsigned candidate_indx[SIZE_CANDIDATE_PIVOT];
  unsigned nb_candidate = 0;
  double candidate_pivot, current_pivot =0.0, candidate_ratio;
  double ratio = INFINITY;
  DEBUG_EXPR_WE(max_pivot_helped = 0;);
  for (unsigned i = 0 ; i < n ; ++i)
  {
    candidate_pivot = col_drive[i];
    if (candidate_pivot > 0.)
    {
      candidate_ratio = q_tilde[i] / candidate_pivot;
      if (candidate_ratio > ratio) continue;
      else if (candidate_ratio < ratio)
      {
        ratio = candidate_ratio;
        block = i;
        current_pivot = candidate_pivot;
        nb_candidate = 0;
        DEBUG_EXPR_WE(max_pivot_helped = 0;);
      }
      else
      {
        if (block == (int)aux_indx || i == aux_indx)
        {
          /* We want the auxilliary variable to exit before any other.
           * see CPS p. 279 and example 4.4.16 */
          block = aux_indx;
          nb_candidate = 0;
          DEBUG_EXPR_WE(max_pivot_helped = 0;);
        }
        else
        {
          if (nb_candidate < SIZE_CANDIDATE_PIVOT-1) candidate_indx[nb_candidate++] = i;
          else
          {
            candidate_indx[nb_candidate] = i;
            DEBUG_PRINTF("pivot_selection_lemke :: lexicomin 5 candidates, ratio = %e\n", ratio);
            for (unsigned k = 0; k < 5; ++k)
            {
              unsigned var = candidate_indx[k];
              double test_pivot = col_drive[var];
              block = lexicosort_rowmajor(var, block, n, lexico_mat, test_pivot, current_pivot, lexico_tol);
              if (block == (int)var)
              {
                current_pivot = test_pivot;
              }
            }
            nb_candidate = 0;
          }
        }
      }
    }
  }
  if (nb_candidate > 0)
  {
    DEBUG_PRINTF("pivot_selection_lemke2 :: lexicomin %d candidate, ratio = %e\n", nb_candidate, ratio);
    for (unsigned k = 0; k < nb_candidate; ++k)
    {
      unsigned var = candidate_indx[k];
      double test_pivot = col_drive[var];
      block = lexicosort_rowmajor(var, block, n, lexico_mat, test_pivot, current_pivot, lexico_tol);
      if (block == (int)var)
      {
        current_pivot = test_pivot;
      }
    }
  }
  DEBUG_EXPR_WE(if (max_pivot_helped) {DEBUG_PRINT("pivot_selection_lemke2 :: lexicomin MAX PIVOT HELPED!\n");});
  return block;
}

int pivot_selection_lemke3(unsigned n, double* restrict col_drive, double* restrict q_tilde, double* restrict lexico_col, unsigned* restrict basis, unsigned* restrict candidate_indx, SN_lumod_dense_data* restrict lumod_data, unsigned aux_indx, double lexico_tol)
{
  int block = -1;
  unsigned nb_candidates = 0;
  double candidate_pivot , candidate_ratio;
  /* double current_pivot; */
  double ratio = INFINITY;
  for (unsigned i = 0 ; i < n ; ++i)
  {
    candidate_pivot = col_drive[i];
    if (candidate_pivot > 0.)
    {
      candidate_ratio = q_tilde[i] / candidate_pivot;
      if (candidate_ratio > ratio) continue;
      else if (candidate_ratio < ratio)
      {
        ratio = candidate_ratio;
        block = i;
        /* current_pivot = candidate_pivot; */
        nb_candidates = 0;
        candidate_indx[0] = i;
      }
      else
      {
        if (block == (int)aux_indx || i == aux_indx)
        {
          /* We want the auxilliary variable to exit before any other.
           * see CPS p. 279 and example 4.4.16 */
          block = aux_indx;
          nb_candidates = 0;
          DEBUG_EXPR_WE(max_pivot_helped = 0;);
        }
        else
        {
          candidate_indx[nb_candidates++] = i;
        }
      }
    }
  }
  if (nb_candidates > 1)
  {
    DEBUG_PRINTF("pivot_selection_lemke3 :: lexicomin %d candidates, ratio = %e\n", nb_candidates, ratio);
    return lexicosort_lumod(candidate_indx, n, nb_candidates, lumod_data, basis, lexico_col, col_drive, lexico_tol);
    DEBUG_EXPR_WE(if (max_pivot_helped) {DEBUG_PRINT("pivot_selection_lemke3 :: lexicomin MAX PIVOT HELPED!\n");});
  }
  else
  {
    return block;
  }
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
  cblas_dcopy_msan(n, q, 1, q_bar, 1);

  /* covering vector for the auxiliary variable */
  double* d = &mat[(n+1)*n];
  if (cov_vec) cblas_dcopy_msan(n, cov_vec, 1, d, 1);
  else for (int i = 0; i < n; ++i) d[i] = 1.0;

  /* take care of M */
  double* mat_basic = &mat[n];
  double* mat_nonbasic = &mat[(n+2)*n];
  for (int i = 0; i < n; ++i)
  {
    if (u[i] > DBL_EPSILON) // M_bas[:, i] = M[:, i]
    {
      basis[i] = i + 2 + n;
      cblas_dcopy_msan(n, &M[i*n], 1, &mat_basic[i*n], 1);
      memset(&mat_nonbasic[i*n], 0, sizeof(double) * n);
      mat_nonbasic[i*n + i] = -1.0;
    }
    else /* s[i] > 0.0 and if both s and u are nonbasic, we insert a column from the identity matrix
          * this choice could ne different
          * M_bas[:, i] = -I[:, i] */
    {
      basis[i] = i + 1;
      cblas_dcopy_msan(n, &M[i*n], 1, &mat_nonbasic[i*n], 1);
      memset(&mat_basic[i*n], 0, sizeof(double) * n);
      mat_basic[i*n + i] = -1.0;
    }
  }

  /* data for LAPACKE */
  lapack_int *ipiv = (lapack_int*)malloc((n+1)*sizeof(lapack_int));
  lapack_int info = 0;

  /* Compute LU factorisation of basis */
  DGETRF(n, n, mat_basic, n, ipiv, &info);

  assert(info <= 0 && "crash_pivot_basis :: info from DGETRF > 0, this should not append !\n");
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
  for (int i = 0; i < n; ++i) mat[i + n*(i + 1)] =  1.0;

  free(ipiv);
  return info;
}


void do_pivot_driftless(double* mat, unsigned int dim, unsigned int dim2, unsigned int block, unsigned int drive)
{
  if (fabs(mat[block + drive*dim]) < DBL_EPSILON)
  {
    printf("do_pivot_driftless :: pivot value too small %e; q[block] = %e; theta = %e\n", mat[block + drive*dim], mat[block], mat[block]/mat[block + drive*dim]);
  }
  double pivot_inv = 1.0/mat[block + drive*dim];
  unsigned ncols = dim*dim2;

#ifdef DEBUG_MESSAGES
  double* driving_col = (double*) malloc(dim*sizeof(double));
  cblas_dcopy(dim, &mat[drive*dim], 1, driving_col, 1);
#endif
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
  DEBUG_PRINT_MAT_SMALL_STR("lexico_mat", (&mat[dim]), dim, dim, driving_col);
#ifdef DEBUG_MESSAGES
  free(driving_col);
#endif
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
  unsigned leaving_var_indx = basis[block] - BASIS_OFFSET;
  unsigned entering_var_indx = drive - BASIS_OFFSET;
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
      assert(lumod_data->row_col_indx[leaving_var_indx] <= 0);
      assert(lumod_data->row_col_indx[entering_var_indx] == 0);
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
        SN_lumod_add_row_col(lumod_data, pos_leaving_var_in_factorized_basis - BASIS_OFFSET, col_tilde);
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
      /*  update the index, since a row and col were deleted */
      if (index_row < (int)lumod_data->k)
      {
        unsigned changed_row = SN_lumod_find_arg_var(lumod_data->row_col_indx, lumod_data->k, n);
        DEBUG_PRINTF("Changing row for variable %d to %d\n", changed_row, index_row);
        lumod_data->row_col_indx[changed_row] = index_row;
      }
      if (-index_col < (int)lumod_data->k)
      {
        unsigned changed_col = SN_lumod_find_arg_var(lumod_data->row_col_indx, -lumod_data->k, n);
        DEBUG_PRINTF("Changing col for variable %d to %d\n", changed_col, index_col);
        lumod_data->row_col_indx[changed_col] = index_col;
      }
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
      assert(lumod_data->row_col_indx[entering_var_indx]>=0);
      assert(lumod_data->row_col_indx[leaving_var_indx] == 0);
      SN_lumod_replace_row(lumod_data, lumod_data->row_col_indx[entering_var_indx], pos_leaving_var_in_factorized_basis-BASIS_OFFSET);
      lumod_data->row_col_indx[leaving_var_indx] = lumod_data->row_col_indx[entering_var_indx];
      lumod_data->row_col_indx[entering_var_indx] = 0;
      break;
    default:
      printf("do_pivot_lumod :: unkown update type occuring\n");
      exit(EXIT_FAILURE);
  }

  /* Update q_tilde */
  double theta = q_tilde[block]/col_drive[block];
  DEBUG_PRINTF("theta = %e; pivot = %e\n", theta, col_drive[block]);
  DEBUG_PRINT_VEC(col_drive, n);
  cblas_daxpy(n, -theta, col_drive, 1, q_tilde, 1);
  q_tilde[block] = theta;

  /* Update the lexico_mat
   * XXX check if this is correct. The value of the pivot may be wrong --xhub */
  double pivot = col_drive[block];
  unsigned block_row_indx = block*n;

  for (unsigned i = 0, j = 0; i < n; ++i, j += n)
  {
    if (j == block_row_indx) continue;
    cblas_daxpy(n, -col_drive[i]/pivot, &lexico_mat[block_row_indx], 1, &lexico_mat[j], 1);
  }
  cblas_dscal(n, 1./pivot, &lexico_mat[block_row_indx], 1);
  DEBUG_PRINT_MAT_ROW_MAJOR_NCOLS_SMALL2_STR("lexico_mat", lexico_mat, n, n, n, col_drive);
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
