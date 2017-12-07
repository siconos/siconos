/* Siconos, Copyright INRIA 2005-2016.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <float.h>

#include "SparseMatrix_internal.h"
#include "SiconosCompat.h"

#if defined(__cplusplus)
#undef restrict
#define restrict __restrict
#endif
//#define VERBOSE_DEBUG

/* add an entry to triplet matrix only if value is not (nearly) null */
CS_INT cs_zentry(CSparseMatrix *T, CS_INT i, CS_INT j, double x)
{
  if(fabs(x) >= DBL_EPSILON)
  {
    return cs_entry(T, i, j, x);
  }
  else
  {
    return 1;
  }
}

int cs_lu_factorization(CS_INT order, const cs *A, double tol, cs_lu_factors * cs_lu_A )
{
  assert(A);
  cs_lu_A->n = A->n;
  css* S = cs_sqr (order, A, 0);
  cs_lu_A->S = S;
  cs_lu_A->N = cs_lu(A, S, tol);

  return (S && cs_lu_A->N);
}

void cs_sparse_free(cs_lu_factors* cs_lu_A)
{
  assert(cs_lu_A);
  if (cs_lu_A)
  {
    cs_lu_A->n = -1;

    cs_sfree(cs_lu_A->S);
    cs_lu_A->S = NULL;

    cs_nfree(cs_lu_A->N);
    cs_lu_A->N = NULL;

    free(cs_lu_A);
  }
}

/* Solve Ax = b with the factorization of A stored in the cs_lu_A
 * This is extracted from cs_lusol, you need to synchronize any changes! */
CS_INT cs_solve (cs_lu_factors* cs_lu_A, double* x, double *b)
{
  assert(cs_lu_A);

  CS_INT ok;
  CS_INT n = cs_lu_A->n;
  css* S = cs_lu_A->S;
  csn* N = cs_lu_A->N;
  ok = (S && N && x) ;
  if (ok)
  {
    cs_ipvec (N->pinv, b, x, n) ;       /* x = b(p) */
    cs_lsolve (N->L, x) ;               /* x = L\x */
    cs_usolve (N->U, x) ;               /* x = U\x */
    cs_ipvec (S->q, x, b, n) ;          /* b(q) = x */
  }
  return (ok);
}

int cs_check_triplet(CSparseMatrix *T)
{
  if (T->nz < 0)
  {
    fprintf(stderr, "cs_check_triplet :: given CSparseMatrix is not in a triplet form: nz = " CS_ID, T->nz);
    return 1;
  }
  CS_INT nb_row = T->m;
  CS_INT nb_col = T->n;
  CS_INT* Ti = T->i;
  CS_INT* Tp = T->p;
  int info = 0;
  CS_INT cc = 0;
  CS_INT max_row = -1;
  CS_INT max_col = -1;

  for (CS_INT indx = 0; indx < T->nz; ++indx)
  {
    if (Ti[indx] >= nb_row)
    {
      printf("cs_check_triplet :: matrix element " CS_ID " has a row number " CS_ID "  > " CS_ID " the number of rows\n", indx, Ti[indx], nb_row);
      info = 1;
    }
    if (Tp[indx] >= nb_col)
    {
      printf("cs_check_triplet :: matrix element " CS_ID " has a row number " CS_ID " > " CS_ID " the number of rows\n", indx, Tp[indx], nb_col);
      info = 1;
    }
    if (Tp[indx] < max_col)
    {
       printf("cs_check_csc :: " CS_ID " at index " CS_ID " > " CS_ID "\n", Tp[indx], indx, max_col);
    }
    if (Tp[indx] == cc)
    {
      if (Ti[indx] <= max_row)
      {
        printf("cs_check_triplet :: matrix element at column " CS_ID " has a row number " CS_ID "  > " CS_ID " the max on that column\n", cc, Ti[indx], max_row);
      }
      else
      {
        max_row = Ti[indx];
      }
    }
    else
    {
      cc = Tp[indx];
      max_row = -1;
    }
  }
  return info;
}

int cs_check_csc(CSparseMatrix *T)
{
  if (T->nz != -1)
  {
    fprintf(stderr, "cs_check_csc :: given CSparseMatrix is not in a csc form: nz = " SN_PTRDIFF_T_F "\n", T->nz);
    return 1;
  }

  CS_INT nb_row = T->m;
  CS_INT nb_col = T->n;
  CS_INT* Ti = T->i;
  CS_INT* Tp = T->p;
  int info = 0;

  for (size_t j = 0; j < nb_col; ++j)
  {
    CS_INT max_indx = -1;
    if (Tp[j] > Tp[j+1])
    {
       printf("cs_check_csc :: " SN_PTRDIFF_T_F " at index " SN_PTRDIFF_T_F " smaller than " SN_PTRDIFF_T_F "\n", Tp[j+1], j, Tp[j]);
    }
    for (size_t p = Tp[j]; p < Tp[j+1]; ++p)
    {
      if (Ti[p] <= max_indx)
      {
        printf("cs_check_csc :: matrix element ("SN_PTRDIFF_T_F","SN_PTRDIFF_T_F") at index " SN_PTRDIFF_T_F " has a row number < " SN_PTRDIFF_T_F " the previous max\n", Ti[p], j, p, max_indx);
        info = 1;
      }
      else if (Ti[p] >= nb_row)
      {
        printf("cs_check_csc :: matrix element ("SN_PTRDIFF_T_F","SN_PTRDIFF_T_F") at index " SN_PTRDIFF_T_F " has a row number > " SN_PTRDIFF_T_F " the max\n", Ti[p], j, p, nb_row);
        info = 1;
      }
      else
      {
        max_indx = Ti[p];
      }
    }
  }
  return info;
}

/* from sparse to dense */
double* cs_dense(CSparseMatrix *A)
{
  CS_INT m = A->m;
  CS_INT n = A->n;
  CS_INT *Ap = A->p;
  CS_INT *Ai = A->i;
  double *Ax = A->x;
  CS_INT nzmax = A->nzmax;
  CS_INT nz = A->nz;
  double *r = (double*) malloc(A->m * A->n * sizeof(double));
  for(int i = 0; i<m*n; ++i)
  {
    r[i] = 0.;
  }
  if(nz < 0)
  {
    for(int j = 0 ; j < n ; j++)
    {
      printf("    col %d : locations " CS_ID " to " CS_ID "\n",
             j, Ap [j], Ap [j+1]-1);
      for(CS_INT p = Ap [j] ; p < Ap [j+1] ; p++)
      {
        printf("      " CS_ID " : %g\n", Ai [p], Ax ? Ax [p] : 1) ;
        r[Ai[p] + j*m] = Ax ? Ax [p] : 1;
      }
    }
  }
  else
  {
    printf("triplet: " CS_ID "-by-" CS_ID ", nzmax: " \
        CS_ID " nnz: " CS_ID "\n", m, n, nzmax, nz) ;
    for(int p = 0 ; p < nz ; p++)
    {
      printf("    " CS_ID " " CS_ID " : %g\n", Ai [p], Ap [p], Ax ? Ax [p] : 1) ;
      r[Ai[p] + Ap[p] * m] = Ax ? Ax[p] : 1;
    }
  }
  return r;
}

/* y = alpha*A*x+beta*y */
int cs_aaxpy(const double alpha, const cs *A, const double * restrict x,
             const double beta, double * restrict y)
{
  CS_INT p, n, *Ap, *Ai ;
  int j;
  double *Ax ;
  if(!A || !x || !y) return (0) ;	     /* check inputs */
  n = A->n ;
  Ap = A->p ;
  Ai = A->i ;
  Ax = A->x ;
  assert(fabs(beta-1.) < 100*DBL_EPSILON && "cs_aaxpy is broken for beta != 1, if you don't see why, ask Olivier");
  for(j = 0 ; j < n ; j++)
  {
    for(p = Ap [j] ; p < Ap [j+1] ; p++)
    {
      y [Ai [p]] *= beta;
      y [Ai [p]] += alpha * Ax [p] * x [j] ;
    }
  }
  return 1;
}

CSparseMatrix* cs_spfree_on_stack(CSparseMatrix* A)
{
  if(!A) return NULL; /* do nothing if A already NULL */
  cs_free(A->p);
  A->p = NULL;
  cs_free(A->i);
  A->i = NULL;
  cs_free(A->x);
  A->x = NULL;
  return NULL;
}

int cs_printInFile(const cs *A, int brief, FILE* file)
{
  CS_INT m, n, nzmax, nz, p, j, *Ap, *Ai ;
  double *Ax ;
  if(!A)
  {
    fprintf(file,"(null)\n") ;
    return (0) ;
  }
  m = A->m ;
  n = A->n ;
  Ap = A->p ;
  Ai = A->i ;
  Ax = A->x ;
  nzmax = A->nzmax ;
  nz = A->nz ;
  fprintf(file,"CSparse Version %d.%d.%d, %s.  %s\n", CS_VER, CS_SUBVER,
         CS_SUBSUB, CS_DATE, CS_COPYRIGHT) ;
  if(nz < 0)
  {
    fprintf(file,"%lld-by-%lld, nzmax: %lld nnz: %lld, 1-norm: %g\n",
            (long long int)m,  (long long int)n,  (long long int)nzmax,
            (long long int)Ap [n],  cs_norm(A)) ;
    for(j = 0 ; j < n ; j++)
    {
      fprintf(file,"    col %lld : locations %lld to %lld\n",  (long long int)j,  (long long int)Ap [j],  (long long int)Ap [j+1]-1);
      for(p = Ap [j] ; p < Ap [j+1] ; p++)
      {
        fprintf(file,"      %lld : %g\n",  (long long int)Ai [p], Ax ? Ax [p] : 1) ;
        if(brief && p > 20)
        {
          fprintf(file,"  ...\n") ;
          return (1) ;
        }
      }
    }
  }
  else
  {
    fprintf(file,"triplet: %lld-by-%lld, nzmax: %lld nnz: %lld\n",  (long long int)m,  (long long int)n,
            (long long int)nzmax,  (long long int)nz) ;
    for(p = 0 ; p < nz ; p++)
    {
      fprintf(file,"    %lld %lld : %g\n",  (long long int)Ai [p],  (long long int)Ap [p], Ax ? Ax [p] : 1) ;
      if(brief && p > 20)
      {
        fprintf(file,"  ...\n") ;
        return (1) ;
      }
    }
  }
  return (1) ;
}


