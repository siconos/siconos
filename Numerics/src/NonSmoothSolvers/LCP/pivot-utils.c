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
#include "pivot-utils.h"
#include "lcp_pivot.h"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
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
        if (dblock < 0) break;
        else if (dblock > 0)
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

int pivot_selection_lemke(double* mat, unsigned int dim, unsigned int drive)
{
  int block = -1;
  double zb, z0, dblock;
  double pivot = 0.0;
  for (unsigned int i = 0 ; i < dim ; ++i)
  {
    zb = mat[i + drive*dim];
    if (zb > DBL_EPSILON)
    {
      z0 = mat[i] / zb;
      if ((block >= 0) && (z0 > pivot)) continue;
      else if ((block == -1) || (z0 < pivot))
      {
        pivot = z0;
        block = i;
      }
      else
      {
        for (unsigned int j = 1; j <= dim; ++j)
        {
          assert(block >= 0 && "pivot_selection_lemke: block < 0");
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
  return block;
}

int pivot_selection_pathsearch(double* mat, unsigned int dim, unsigned int drive, unsigned int t_indx)
{
  int block = pivot_selection_lemke(mat, dim, drive);
  double zb = mat[t_indx + drive*dim];
  if (zb <= -DBL_EPSILON)
  {
    /* try to set t to 1 */
    double z0 = (mat[t_indx] - 1.0)/zb;
    if (block >= 0)
    {
      double pivot_lemke = mat[block]/mat[block + drive*dim];
      if (z0 <= pivot_lemke)
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
  double tmp;
  double pivot_inv = 1.0/mat[block + drive*dim];

  /* Update column mat[block, :] */
  mat[block + drive*dim] = 1.0; /* nm_rs = 1 */
  /* nm_rj = m_rj/m_rs */
  for (unsigned int i = 0        ; i < drive; ++i) mat[block + i*dim] *= pivot_inv;
  for (unsigned int i = drive + 1; i < dim2 ; ++i) mat[block + i*dim] *= pivot_inv;

  /* Update other columns*/
  for (unsigned int i = 0; i < block; ++i)
  {
    tmp = mat[i + drive*dim];
    /* nm_ij = m_ij + (m_ir/m_rs)m_rj = m_ij - m_is*nm_rj */
    for (unsigned int j = 0; j < dim2; ++j) mat[i + j*dim] -= tmp*mat[block + j*dim];
  }
  for (unsigned int i = block + 1; i < dim; ++i)
  {
    tmp = mat[i + drive*dim];
    /* nm_ij = m_ij + (m_ir/m_rs)m_rj = m_ij - m_is*nm_rj */
    for (unsigned int j = 0; j < dim2; ++j) mat[i + j*dim] -= tmp*mat[block + j*dim];
  }
}

void do_pivot_driftless2(double* mat, unsigned int dim, unsigned int dim2, unsigned int block, unsigned int drive)
{
  double tmp;
  double pivot_inv = 1.0/mat[block + drive*dim];

  /* Update column mat[block, :] */
  mat[block + drive*dim] = 1.0; /* nm_rs = 1 */
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
