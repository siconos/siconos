/* Copyright 2021 INRIA.
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
#define _POSIX_C_SOURCE 200112L

#include "CSparseMatrix_internal.h"
#include <assert.h>            // for assert
#include <cs.h>                // for CS_INT, CS_ID, cs_dl_spalloc, cs_dl_free
#include <float.h>             // for DBL_EPSILON
#include <math.h>              // for fabs
#include <stdio.h>             // for fprintf, sscanf, printf, NULL, fgets
#include <stdlib.h>            // for realloc, exit, free, malloc, EXIT_FAILURE
#include <string.h>            // for strtok_r, memcpy, strncmp
#include "SiconosCompat.h"     // for SN_PTRDIFF_T_F
#include "numerics_verbose.h"  // for CHECK_IO
#define LDL_LONG
#include "ldl.h"
#if defined(__cplusplus)
#undef restrict
#include <sys/cdefs.h> // for __restrict
#define restrict __restrict
#endif
/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "siconos_debug.h" // for DEBUG_PRINTF


#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"
#endif

double CSparseMatrix_get_value(const CSparseMatrix *A, CS_INT i, CS_INT j)
{
  CS_INT * Ai =   A->i;
  CS_INT * Ap =   A->p;
  double * Ax =   A->x;

  for(CS_INT row = Ap[j]; row < Ap[j+1] ; row++)
  {
    if(i == Ai[row])
      return  Ax[row];
  }
  return 0.0;
}

void CSparseMatrix_write_in_file_python(const CSparseMatrix* const m, FILE* file)
{

  fprintf(file, "m = %ld; \n", m->m);
  fprintf(file, "n = %ld; \n", m->n);
  fprintf(file, "data= [");
  for(int i = 0; i < m->m; i++)
  {
    fprintf(file, "[");
    for(int j = 0; j < m->n; j++)
    {
      fprintf(file, "%32.24e,\t ", CSparseMatrix_get_value((CSparseMatrix*) m,i,j));
    }
    fprintf(file, "],\n");
  }
  fprintf(file, "]");
}



/* y = alpha*A*x+beta*y */
int CSparseMatrix_aaxpby(const double alpha, const CSparseMatrix *A,
                         const double *restrict x,
                         const double beta, double *restrict y)
{

  CS_INT n, m, *Ap, *Ai ;
  double *Ax ;
  if(!CS_CSC(A) || !x || !y) return (0);	     /* check inputs */
  {
    n = A->n;
    m = A->m;
    Ap = A->p;
    Ai = A->i;
    Ax = A->x;


    for(int j=0; j<m; j++)
    {
      y[j] *= beta;
    }

    for(int j=0 ; j<n ; j++)
    {
      for(CS_INT p = Ap [j] ; p < Ap [j+1] ; p++)
      {
        y [Ai [p]] += alpha * Ax [p] * x [j];
      }
    }

  }
  return 1;

}
/* A <-- alpha*A */
int CSparseMatrix_scal(const double alpha, const CSparseMatrix *A)
{

  CS_INT n, *Ap;
  double *Ax ;
  if(!CS_CSC(A)) return (0);	     /* check inputs */
  {
    n = A->n;
    Ap = A->p;
    Ax = A->x;
    for(int j=0 ; j<n ; j++)
    {
      for(CS_INT p = Ap [j] ; p < Ap [j+1] ; p++)
      {
        Ax [p] = alpha * Ax [p];
      }
    }

  }
  return 1;
}

int CSparseMatrix_check_triplet(CSparseMatrix *T)
{
  if(T->nz < 0)
  {
    fprintf(stderr, "CSparseMatrix_check_triplet :: given CSparseMatrix is not in a triplet form: nz = " CS_ID, T->nz);
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

  for(CS_INT indx = 0; indx < T->nz; ++indx)
  {
    if(Ti[indx] >= nb_row)
    {
      printf("CSparseMatrix_check_triplet :: matrix element " CS_ID " has a row number " CS_ID "  > " CS_ID " the number of rows\n", indx, Ti[indx], nb_row);
      info = 1;
    }
    if(Tp[indx] >= nb_col)
    {
      printf("CSparseMatrix_check_triplet :: matrix element " CS_ID " has a row number " CS_ID " > " CS_ID " the number of rows\n", indx, Tp[indx], nb_col);
      info = 1;
    }
    if(Tp[indx] < max_col)
    {
      printf("CSparseMatrix_check_csc :: " CS_ID " at index " CS_ID " > " CS_ID "\n", Tp[indx], indx, max_col);
    }
    if(Tp[indx] == cc)
    {
      if(Ti[indx] <= max_row)
      {
        printf("CSparseMatrix_check_triplet :: matrix element at column " CS_ID " has a row number " CS_ID "  > " CS_ID " the max on that column\n", cc, Ti[indx], max_row);
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

int CSparseMatrix_check_csc(CSparseMatrix *T)
{
  if(T->nz != -1)
  {
    fprintf(stderr, "CSparseMatrix_check_csc :: given CSparseMatrix is not in a csc form: nz = " SN_PTRDIFF_T_F "\n", T->nz);
    return 1;
  }

  CS_INT nb_row = T->m;
  CS_INT nb_col = T->n;
  CS_INT* Ti = T->i;
  CS_INT* Tp = T->p;
  int info = 0;

  for(CS_INT j = 0; j < nb_col; ++j)
  {
    CS_INT max_indx = -1;
    if(Tp[j] > Tp[j+1])
    {
      printf("CSparseMatrix_check_csc :: " SN_PTRDIFF_T_F " at index " SN_PTRDIFF_T_F " smaller than " SN_PTRDIFF_T_F "\n", Tp[j+1], j, Tp[j]);
    }
    for(CS_INT p = Tp[j]; p < Tp[j+1]; ++p)
    {
      if(Ti[p] <= max_indx)
      {
        printf("CSparseMatrix_check_csc :: matrix element (" SN_PTRDIFF_T_F"," SN_PTRDIFF_T_F") at index " SN_PTRDIFF_T_F " has a row number < " SN_PTRDIFF_T_F " the previous max\n", Ti[p], j, p, max_indx);
        info = 1;
      }
      else if(Ti[p] >= nb_row)
      {
        printf("CSparseMatrix_check_csc :: matrix element (" SN_PTRDIFF_T_F"," SN_PTRDIFF_T_F") at index " SN_PTRDIFF_T_F " has a row number > " SN_PTRDIFF_T_F " the max\n", Ti[p], j, p, nb_row);
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

CSparseMatrix* CSparseMatrix_spfree_on_stack(CSparseMatrix* A)
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

int CSparseMatrix_lu_factorization(CS_INT order, const cs *A, double tol, CSparseMatrix_factors * cs_lu_A)
{
  assert(A);
  cs_lu_A->n = A->n;
  css* S = cs_sqr(order, A, 0);
  cs_lu_A->S = S;
  cs_lu_A->N = cs_lu(A, S, tol);

  return (S && cs_lu_A->N);
}
int CSparseMatrix_chol_factorization(CS_INT order, const cs *A,  CSparseMatrix_factors * cs_chol_A)
{
  assert(A);
  FILE * foutput = fopen("cs.py", "w");
  CSparseMatrix_print_in_file(A, 0, foutput);
  fclose(foutput);
  cs_chol_A->n = A->n;
  css* S = cs_schol(order, A);
  cs_chol_A->S = S;
  cs_chol_A->N = cs_chol(A, S);

  return (S && cs_chol_A->N);
}
int CSparseMatrix_ldlt_factorization(CS_INT order, const cs *A,  CSparseMatrix_factors * cs_ldlt_A)
{
  assert(A);

  CS_INT *Ap, *Ai, *Lp, *Li;
  CS_INT *Parent;
  CS_ENTRY *Ax, *Lx;
  css *S;
  csn *N;
  CS_INT n, lnz;
  DEBUG_EXPR(cs_print(A,1););
  Ap=A->p;
  Ai=A->i;
  Ax=A->x;

  cs_ldlt_A->n = n =  A->n;

  /* could be good to good cs_spalloc */

  cs_ldlt_A->N = N        =  cs_calloc (1, sizeof (csn)) ;       /* allocate result N */
  cs_ldlt_A->N->L         = cs_calloc (1, sizeof (cs)) ;         /* allocate the cs struct */
  cs_ldlt_A->N->L->m =n;
  cs_ldlt_A->N->L->n =n;
  cs_ldlt_A->N->L->nz =-1;
  cs_ldlt_A->N->L->p = Lp = cs_malloc (n+1, sizeof (CS_INT)) ;


  cs_ldlt_A->S = S = cs_calloc (1, sizeof (css)) ;              /* allocate result S */
  cs_ldlt_A->S->parent = Parent = cs_malloc (n+1, sizeof (CS_INT)) ;

  CS_INT* Lnz = cs_malloc (n, sizeof (CS_INT)) ;
  CS_INT* Flag =  cs_malloc (n, sizeof (CS_INT)) ;

  /* ordering with amd */
  CS_INT * Perm, *PermInv;
  cs_ldlt_A->N->pinv = Perm  = cs_amd (order, A) ;  /* We used pinv to store Perm !! */
  PermInv = cs_malloc (n, sizeof (CS_INT)) ;


  DEBUG_EXPR(for (int k =0; k< n+1; k++){printf("%li\t", Perm[k]);}printf("\n"););

  /* symbolic factorization to get Lp, Parent, Lnz, and Pinv */
  LDL_symbolic (n, Ap, Ai, Lp, Parent, Lnz, Flag, Perm, PermInv) ;
  DEBUG_EXPR(for (int k =0; k< n+1; k++){printf("%li\t", Lp[k]);}printf("\n"););
  DEBUG_EXPR(for (int k =0; k< n; k++){printf("%li\t", Flag[k]);}printf("\n"););

  /* factorization */
  lnz = Lp [n] ;
  DEBUG_PRINTF("Lp[n] = %ld\n", Lp[n]);
  cs_ldlt_A->N->L->i = Li = cs_malloc (lnz, sizeof (CS_INT)) ;
  cs_ldlt_A->N->L->x = Lx = cs_malloc (lnz, sizeof (CS_ENTRY)) ;

  CS_INT *Pattern =  cs_malloc (n, sizeof (CS_INT)) ;
  CS_ENTRY* D;
  cs_ldlt_A->N->B = D= cs_malloc (n, sizeof (CS_ENTRY)) ; /* We use cs_ldlt_A->N->B  for storing D !! */
  CS_ENTRY* Y = cs_malloc (n, sizeof (CS_ENTRY)) ;
  LDL_numeric (n, Ap, Ai, Ax, Lp, Parent, Lnz, Li, Lx, D,
               Y, Flag, Pattern, Perm, PermInv) ;

  DEBUG_EXPR(cs_print(cs_ldlt_A->N->L,1););
  DEBUG_EXPR(NV_display(D,n));


  cs_free(Lnz);
  cs_free(Flag);
  cs_free(PermInv);
  cs_free(Pattern);
  cs_free(Y);

  return (S && cs_ldlt_A->N);
}

void CSparseMatrix_free_lu_factors(CSparseMatrix_factors* cs_lu_A)
{
  assert(cs_lu_A);
  if(cs_lu_A)
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
CS_INT CSparseMatrix_solve(CSparseMatrix_factors* cs_lu_A, double* x, double *b)
{
  assert(cs_lu_A);

  CS_INT ok;
  CS_INT n = cs_lu_A->n;
  css* S = cs_lu_A->S;
  csn* N = cs_lu_A->N;
  ok = (S && N && x) ;
  if(ok)
  {
    cs_ipvec(N->pinv, b, x, n) ;        /* x = b(p) */
    cs_lsolve(N->L, x) ;                /* x = L\x */
    cs_usolve(N->U, x) ;                /* x = U\x */
    cs_ipvec(S->q, x, b, n) ;           /* b(q) = x */
  }
  return (ok);
}

/* Solve Ax = B with the factorization of A stored in the cs_lu_A
 * B is a sparse matrix (CSparseMatrix_factors)
 * This is extracted from cs_lusol, you need to synchronize any changes! */
CS_INT CSparseMatrix_spsolve(CSparseMatrix_factors* cs_lu_A,  CSparseMatrix* X, CSparseMatrix* B)
{
  assert(cs_lu_A);
  DEBUG_BEGIN("CSparseMatrix_spsolve(...)\n");

  if(!CS_CSC(X)) return 1 ;                  /* check inputs */
  if(!CS_CSC(B)) return 1 ;                  /* check inputs */

  CS_INT ok;
  CS_INT n = cs_lu_A->n;
  csn* N = cs_lu_A->N;
  if(!N) return 1;
  css* S = cs_lu_A->S;
  if(!S) return 1;

  CS_ENTRY *x, *b, *Xx, *Bx ;
  CS_INT *xi, *pinv, *q, top, k, i, p, *Bp, *Bi, *Xp, *Xi;

  x = cs_malloc(n, sizeof(CS_ENTRY)) ;              /* get CS_ENTRY workspace */
  b = cs_malloc(n, sizeof(CS_ENTRY)) ;              /* get CS_ENTRY workspace */
  xi = cs_malloc(2*n, sizeof(CS_INT)) ;             /* get CS_INT workspace */

  CS_INT xnz=0;
  Xp= X->p;
  Xi= X->i;
  Xx= X->x;
  Xp[0]=0;
  pinv = N->pinv;
  /* --- 1. First step X = L\B ---------------------------------------------- */
  for(k = 0 ; k < B->n ; k++)
  {
    /* permutation of the rows  of B(:,k) */
    for(p = B->p [k] ; p < B->p [k+1] ; p++) x [B->i[p]] = B->x [p] ; /* scatter B  */
    for(p = B->p [k] ; p < B->p [k+1] ; p++)
    {
      CS_INT i_old= B->i[p];
      B->i[p] = pinv[i_old]; /* permute row indices with N->pinv */
      B->x[p] = x[i_old];
    }
    /* call spsolve */
    top = cs_spsolve(N->L, B, k, xi, x, NULL, 1) ;    /* x = L\B(:,col) */
    /* store the result in X */
    if(Xp[k]+ n-top > X->nzmax && !cs_sprealloc(X, 2*(X->nzmax)+ n-top))    /* realloc X if need */
    {
      return 1;  /* (cs_done(X, w, x, 0)) ;   */            /* out of memory */
    }
    Xp= X->p;
    Xi= X->i;
    Xx= X->x;
    for(p = top ; p < n ; p++)
    {
      i = xi [p] ;/* x(i) is nonzero */
      Xi[xnz]=i;          /* store the result in X */
      Xx[xnz++] = x[i];
    }
    Xp[k+1] =Xp[k] + n-top;
  }
  DEBUG_EXPR(cs_print(X,0););

  /* --- 2. First step B = U\B ---------------------------------------------- */
  DEBUG_PRINT("2- Second step B = U\\X\n");
  CS_INT bnz=0;
  Bp= B->p;
  Bi= B->i;
  Bx= B->x;
  Bp[0]=0;
  q = S->q;
  for(k = 0 ; k < X->n ; k++)
  {
    top = cs_spsolve(N->U, X, k, xi, x, NULL, 0) ;   /* x = U\X(:,col) */

    /* store the result in B */
    if(Bp[k]+ n-top > B->nzmax && !cs_sprealloc(B, 2*(B->nzmax)+ n-top))
    {
      return 1;  /* (cs_done(X, w, x, 0)) ;   */            /* out of memory */
    }
    Bp= B->p;
    Bi= B->i;
    Bx= B->x;
    for(p = top ; p < n ; p++)
    {
      i = xi[p] ;           /* x(i) is nonzero */
      Bi[bnz] = q[i];       /* permute the indices with S->q */
      Bx[bnz++] = x[i];
    }
    Bp[k+1] = Bp[k] + n-top;
  }
  DEBUG_EXPR(cs_print(B,0));
  ok =1;
  free(x);
  free(b);
  free(xi);
  DEBUG_END("CSparseMatrix_spsolve(...)\n");
  return (ok);
}





CS_INT CSparseMatrix_chol_solve(CSparseMatrix_factors* cs_chol_A, double* x, double *b)
{
  assert(cs_chol_A);
  
  CS_INT ok;
  CS_INT n = cs_chol_A->n;
  css* S = cs_chol_A->S;
  csn* N = cs_chol_A->N;
  ok = (S && N && x) ;
  if(ok)
  {
    cs_ipvec(S->pinv, b, x, n) ;    /* x = P*b */
    cs_lsolve(N->L, x) ;            /* x = L\x */
    cs_ltsolve(N->L, x) ;           /* x = L'\x */
    cs_pvec(S->pinv, x, b, n) ;     /* b = P'*x */
  }
  return (ok);
}



/* Solve Ax = B with the factorization of A stored in the cs_chol_A
 * B is a sparse matrix (CSparseMatrix_factors)
 * This is extracted from cs_lusol, you need to synchronize any changes! */
CS_INT CSparseMatrix_chol_spsolve(CSparseMatrix_factors* cs_chol_A,  CSparseMatrix* X, CSparseMatrix* B)
{
  DEBUG_BEGIN("CSparseMatrix_chol_spsolve(...)\n");

  if(!CS_CSC(X)) return 1 ;                  /* check inputs */
  if(!CS_CSC(B)) return 1 ;                  /* check inputs */

  CS_INT ok;
  CS_INT n = cs_chol_A->n;
  csn* N = cs_chol_A->N;
  if(!N) return 1;
  css* S = cs_chol_A->S;
  if(!S) return 1;

  CS_ENTRY *x, *b, *Xx, *Bx ;
  CS_INT *xi, *pinv, top, k, i, p, *Bp, *Bi, *Xp, *Xi;

  x = cs_malloc(n, sizeof(CS_ENTRY)) ;              /* get CS_ENTRY workspace */
  b = cs_malloc(n, sizeof(CS_ENTRY)) ;              /* get CS_ENTRY workspace */
  xi = cs_malloc(2*n, sizeof(CS_INT)) ;             /* get CS_INT workspace */

  CS_INT xnz=0;
  Xp= X->p;
  Xi= X->i;
  Xx= X->x;
  Xp[0]=0;
  pinv = S->pinv;

  /* for (i =0; i <n; i++) printf("pinv[%li] = %li\t", i, pinv[i]); */
  /* printf("\n"); */

  /* --- 1. First step X = L\B ---------------------------------------------- */
  for(k = 0 ; k < B->n ; k++)
  {
    /* permutation of the rows  of B(:,k) */
    for(p = B->p [k] ; p < B->p [k+1] ; p++) x [B->i[p]] = B->x [p] ; /* scatter B  */
    for(p = B->p [k] ; p < B->p [k+1] ; p++)
    {
      CS_INT i_old= B->i[p];
      B->i[p] = pinv[i_old]; /* permute row indices with N->pinv */
      B->x[p] = x[i_old];
    }
    /* call spsolve */
    top = cs_spsolve(N->L, B, k, xi, x, NULL, 1) ;    /* x = L\B(:,col) */
    /* store the result in X */
    if(Xp[k]+ n-top > X->nzmax && !cs_sprealloc(X, 2*(X->nzmax)+ n-top))    /* realloc X if need */
    {
      return 1;  /* (cs_done(X, w, x, 0)) ;   */            /* out of memory */
    }
    Xp= X->p;
    Xi= X->i;
    Xx= X->x;
    for(p = top ; p < n ; p++)
    {
      i = xi [p] ;/* x(i) is nonzero */
      Xi[xnz]=i;          /* store the result in X */
      Xx[xnz++] = x[i];
    }
    Xp[k+1] =Xp[k] + n-top;
  }
  DEBUG_EXPR(cs_print(X,0););

  /* --- 2. First step B = L'\B ---------------------------------------------- */
  DEBUG_PRINT("2- Second step B = L'\\X\n");

  CSparseMatrix* LT = cs_transpose(N->L,1);
  CS_INT bnz=0;
  Bp= B->p;
  Bi= B->i;
  Bx= B->x;
  Bp[0]=0;
  for(k = 0 ; k < X->n ; k++)
  {
    top = cs_spsolve(LT, X, k, xi, x, NULL, 0) ;   /* x = U\X(:,col) */

    /* store the result in B */
    if(Bp[k]+ n-top > B->nzmax && !cs_sprealloc(B, 2*(B->nzmax)+ n-top))
    {
      return 1;  /* (cs_done(X, w, x, 0)) ;   */            /* out of memory */
    }
    Bp= B->p;
    Bi= B->i;
    Bx= B->x;
    /* permutation with S->pinv */
    cs_pvec(pinv, x, b, n) ;     /* b = P'*x */
    //for(p = top ; p < n ; p++) b [q ? q [p] : p] = x [p] ;  /* b(q) = x */
    for(p = top ; p < n ; p++)
    {
      i = xi[p] ;               /* x(i) is nonzero */
      i = pinv[i];              /* permute the indices with S->pinv */ // to be checked carefully
      Bi[bnz] = i;
      Bx[bnz++] = b[i];
    }
    Bp[k+1] = Bp[k] + n-top;
  }
  DEBUG_EXPR(cs_print(B,0));
  ok =1;
  free(x);
  free(b);
  free(xi);
  DEBUG_END("CSparseMatrix_chol_spsolve(...)\n");
  return (ok);
}


CS_INT CSparseMatrix_ldlt_solve(CSparseMatrix_factors* cs_ldlt_A, double* x, double *b)
{
  assert(cs_ldlt_A);

  CS_INT ok;

  CS_INT n = cs_ldlt_A->n;
  css* S = cs_ldlt_A->S;
  csn* N = cs_ldlt_A->N;
  ok = (S && N && x) ;


  cs* L;
  CS_INT  *Lp, *Li;
  CS_ENTRY *Lx;

  if(ok)
  {
    CS_INT * P = cs_ldlt_A->N->pinv; /* We used pinv to store Perm !! */
    CS_ENTRY* D = cs_ldlt_A->N->B;   /* We use cs_ldlt_A->N->B  for storing D !! */

    L = cs_ldlt_A->N->L;
    Lp=L->p;
    Li=L->i;
    Lx=L->x;

    /* solve Ax=b */
		if (P)
		{
      /* the factorization is LDL' = PAP' */
      LDL_perm (n, x, b, P) ;			/* y = Pb */
      LDL_lsolve (n, x, Lp, Li, Lx) ;		/* y = L\y */
      LDL_dsolve (n, x, D) ;			/* y = D\y */
      LDL_ltsolve (n, x, Lp, Li, Lx) ;		/* y = L'\y */
      LDL_permt (n, b, x, P) ;			/* x = P'y */
		}
		else
		{
      /* the factorization is LDL' = A */
      DEBUG_EXPR(NV_display(b,n));
      DEBUG_EXPR(NV_display(D,n));
      LDL_lsolve (n, b, Lp, Li, Lx) ;		/* x = L\x */
      LDL_dsolve (n, b, D) ;			/* x = D\x */
      LDL_ltsolve (n, b, Lp, Li, Lx) ;		/* x = L'\x */
      DEBUG_EXPR(NV_display(b,n));
		}


  }
  return (ok);
}

CSparseMatrix * CSparseMatrix_new_from_file(FILE* file)
{
  CS_INT m=0, n=0, nzmax=0, nz, p, j, *Ap, *Ai ;
  long long foo;
  double *Ax ;
  int info = 0;
  char line[2048];

  /* CHECK_IO(fscanf(file, "%20[^\n]", line ), &info); */
  /* fscanf(file, "%2047[^\n]", line ); */

  if(fgets(line, 2047, file)!=NULL)
  {
    DEBUG_PRINTF("line = %s\n",line);
  }
  if(fgets(line, 2047, file)!=NULL)
  {
    DEBUG_PRINTF("line = %s\n",line);
  }
  if(fgets(line, 2047, file)!=NULL)
  {
    DEBUG_PRINTF("line = %s\n",line);
  }

  char *str1, *str2, *token, *subtoken;
  char *saveptr1, *saveptr2;
  const char s_1[2] = " ", s_2[2] = "-";
  int k;
  int is_triplet = 0;
  for(k = 1, str1 = line; ; k++, str1 = NULL)
  {
    token = strtok_r(str1, s_1, &saveptr1);
    if(token == NULL)
      break;
    DEBUG_PRINTF("%d: %s\n", k, token);
    if(strncmp(token, "triplet:",8) == 0)
    {
      DEBUG_PRINT(" triplet matrix\n");
      is_triplet =1;
    }
    if(k==1+is_triplet)
    {
      int kk =0;
      for(str2 = token; ; str2 = NULL)
      {
        subtoken = strtok_r(str2, s_2, &saveptr2);
        if(kk==0)
        {
          if(1 == sscanf(subtoken, "%lld", &foo))
            m = (CS_INT)foo;
        }
        if(kk==2)
        {
          if(1 == sscanf(subtoken, "%lld", &foo) && kk)
            n = (CS_INT)foo;
        }
        kk=kk+1;
        if(subtoken == NULL)
          break;
        DEBUG_PRINTF(" --> %s\n", subtoken);
      }
      DEBUG_PRINTF("m = %li, n = %li \n", m, n);
    }

    if(k==3+is_triplet)
    {
      if(1 == sscanf(token, "%lld", &foo))
        nzmax = (CS_INT)foo;
    }
    if(k==6  && is_triplet)
    {
      if(1 == sscanf(token, "%lld", &foo))
        nz = (CS_INT)foo;
    }
  }

  CSparseMatrix * out = cs_spalloc(m, n, nzmax, 1, is_triplet);

  if(is_triplet)
  {
    out->nz=nz;
  }
  else
  {
    out->nz=-1;
  }

  Ai = out->i;
  Ap = out->p;
  Ax = out->x;

  long long int val1, val2;
  double val3;
  if(out->nz < 0)
  {
    for(j = 0 ; j < n ; j++)
    {
      /* fprintf(file,"    col %lld : locations %lld to %lld\n",  (long long int)j,  (long long int)Ap [j],  (long long int)Ap [j+1]-1); */

      if(fgets(line, 2047, file)!=NULL)
      {
        DEBUG_PRINTF("line 1 = %s\n",line);
      }
      for(k = 1, str1 = line; ; k++, str1 = NULL)
      {
        token = strtok_r(str1, s_1, &saveptr1);
        if(token == NULL)
          break;
        DEBUG_PRINTF("%d: %s\n", k, token);
        if(k==2)
        {
          if(1 == sscanf(token, "%lld", &val1))
          {
            DEBUG_PRINTF(" j- col = %lli\n", j-val1);
          }
          assert(j-val1 == 0);
        }
        if(k==5)
        {
          if(1 == sscanf(token, "%lld", &val1))
            Ap [j] = (CS_INT)val1;
        }
        if(k==7)
        {
          if(1 == sscanf(token, "%lld", &val2))
            Ap [j+1] = (CS_INT)val2+1;
        }
      }
      DEBUG_PRINTF("    col %lld : locations %lld to %lld\n", (long long int)j, (long long int)Ap [j], (long long int)Ap [j+1]-1);

      for(p = Ap [j] ; p < Ap [j+1] ; p++)
      {
        if(fgets(line, 2047, file)!=NULL)
        {
          DEBUG_PRINTF("line 2 = %s\n",line);
        }

        for(k = 1, str1 = line; ; k++, str1 = NULL)
        {
          token = strtok_r(str1, s_1, &saveptr1);
          if(token == NULL)
            break;
          DEBUG_PRINTF("%d: %s\n", k, token);
          if(k==1)
          {
            if(1 == sscanf(token, "%lld", &val1))
              Ai[p] = (CS_INT)val1;
          }
          if(k==3)
          {
            if(1 == sscanf(token, "%32le", &val3))
              Ax[p] = val3;
          }
        }
      }
    }
  }
  else
  {
    /* fprintf(file,"triplet: %lld-by-%lld, nzmax: %lld nnz: %lld\n",  (long long int)m,  (long long int)n, */
    /*         (long long int)nzmax,  (long long int)nz) ; */


    for(p = 0 ; p < nz ; p++)
    {
      CHECK_IO(fscanf(file,"    %lld %lld : %32le\n", &val1,  &val2, &val3), &info);
      DEBUG_PRINTF("    %lld %lld : %32le\n", val1,  val2, val3);
      Ai [p] = (CS_INT) val1;
      Ap [p] = (CS_INT) val2;
      Ax [p] = val3;
      DEBUG_PRINTF("    %lld %lld : %g\n", (long long int)Ai [p], (long long int)Ap [p], (double)Ax[p]);
    }
  }
  return out ;
}


/* add an entry to triplet matrix only if value is not (nearly) null */
CS_INT CSparseMatrix_zentry(CSparseMatrix *T, CS_INT i, CS_INT j, double x, double threshold)
{
  if(fabs(x) >= threshold)
  {
    return cs_entry(T, i, j, x);
  }
  else
  {
    return 1;
  }
}

/* add an entry to a symmetric triplet matrix only if value is not (nearly) null */
CS_INT CSparseMatrix_symmetric_zentry(CSparseMatrix *T, CS_INT i, CS_INT j, double x, double threshold)
{
  if(fabs(x) >= threshold)
  {
    if(j<=i)
    {
      return cs_entry(T, i, j, x);
    }
  }
  return 1;
}


int CSparseMatrix_print(const CSparseMatrix *A, int brief)
{
  cs_print(A, brief);
  return 0;
}

/* add an entry to triplet matrix */
CS_INT CSparseMatrix_entry(CSparseMatrix *T, CS_INT i, CS_INT j, double x)
{
    return cs_entry(T, i, j, x);
}

/* add an entry to a symmetric triplet matrix */
CS_INT CSparseMatrix_symmetric_entry(CSparseMatrix *T, CS_INT i, CS_INT j, double x)
{
  if(j<=i)
  {
    return cs_entry(T, i, j, x);
  }
  return 1;
}


int CSparseMatrix_print_in_file(const CSparseMatrix *A, int brief, FILE* file)
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
            (long long int)m, (long long int)n, (long long int)nzmax,
            (long long int)Ap [n],  cs_norm(A)) ;
    for(j = 0 ; j < n ; j++)
    {
      fprintf(file,"    col %lld : locations %lld to %lld\n", (long long int)j, (long long int)Ap [j], (long long int)Ap [j+1]-1);
      for(p = Ap [j] ; p < Ap [j+1] ; p++)
      {
        fprintf(file,"      %lld : %g\n", (long long int)Ai [p], Ax ? Ax [p] : 1) ;
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
    fprintf(file,"triplet: %lld-by-%lld, nzmax: %lld nnz: %lld\n", (long long int)m, (long long int)n,
            (long long int)nzmax, (long long int)nz) ;
    for(p = 0 ; p < nz ; p++)
    {
      fprintf(file,"    %lld %lld : %g\n", (long long int)Ai [p], (long long int)Ap [p], Ax ? Ax [p] : 1) ;
      if(brief && p > 20)
      {
        fprintf(file,"  ...\n") ;
        return (1) ;
      }
    }
  }
  return (1) ;
}

CS_INT CSparseMatrix_to_dense(const CSparseMatrix* const A, double * B)
{

  CS_INT p, j, m, n, nz, *Ap, *Ai ;
  CS_ENTRY *Ax ;

  if(!A)
  {
    printf("CSparseMatrix_to_dense :: A = null\n") ;
    return (1) ;
  }

  m = A->m ;
  n = A->n ;
  nz = A->nz ;
  Ap = A->p ;
  Ai = A->i ;
  Ax = A->x ;
  Ax = A->x;

  if(nz >= 0)
  {
    for(p = 0 ; p < nz ; p++)
    {
      DEBUG_PRINTF("    %g %g : ", (double)(Ai [p]), (double)(Ap [p])) ;
      DEBUG_PRINTF("%g\n", Ax ? Ax [p] : 1) ;
      DEBUG_PRINTF("B %g\n", B[Ai[p] + Ap[p]*m]) ;
      B[Ai[p] + Ap[p]*m] = Ax [p];
    }
  }
  else
  {
    for(j = 0 ; j < n ; j++)
    {
      DEBUG_PRINTF("    col %g : locations %g to %g\n", (double) j,
                   (double)(Ap [j]), (double)(Ap [j+1]-1)) ;
      for(p = Ap [j] ; p < Ap [j+1] ; p++)
      {
        DEBUG_PRINTF("      %g %g : ", (double)(Ai [p]), Ax ? Ax [p] : 1) ;
        B[Ai[p]+ j * m] = Ax [p];
      }
    }
  }
  return (0) ;
}

CSparseMatrix* CSparseMatrix_alloc_for_copy(const CSparseMatrix* const m)
{
  assert(m);
  CSparseMatrix* out = NULL;
  if(m->nz >= 0)  /* triplet  */
  {
    out = cs_spalloc(m->m, m->n, m->nzmax, 1, 1);
  }
  else if(m->nz == -1)  /* csc */
  {
    out = cs_spalloc(m->m, m->n, m->nzmax, 1, 0);
  }
  else if(m->nz == -2)  /* csr  */
  {
    out = cs_spalloc(m->n, m->m, m->nzmax, 1, 0);
    out->nz = -2;
    out->m = m->m;
    out->n = m->n;
  }
  else
  {
    fprintf(stderr, "NM_copy :: error unknown type " CS_ID
            " for CSparse matrix\n", m->nz);
    exit(EXIT_FAILURE);
  }

  return out;
}

void CSparseMatrix_copy(const CSparseMatrix* const A, CSparseMatrix* B)
{
  assert(A);
  assert(B);

  if(B->nzmax < A->nzmax)
  {
    B->x = (double *) realloc(B->x, A->nzmax * sizeof(double));
    B->i = (CS_INT *) realloc(B->i, A->nzmax * sizeof(CS_INT));
  }
  else if(!(B->x))
  {
    B->x = (double *) malloc(A->nzmax * sizeof(double));
  }

  if(A->nz >= 0)
  {
    /* triplet */
    B->p = (CS_INT *) realloc(B->p, A->nzmax * sizeof(CS_INT));
  }
  else if((A->nz == -1) && (B->n < A->n))
  {
    /* csc */
    B->p = (CS_INT *) realloc(B->p, (A->n + 1) * sizeof(CS_INT));
  }
  else if((A->nz == -2) && (B->m < A->m))
  {
    /* csr */
    B->p = (CS_INT *) realloc(B->p, (A->m + 1) * sizeof(CS_INT));
  }


  B->nzmax = A->nzmax;
  B->nz = A->nz;
  B->m = A->m;
  B->n = A->n;

  memcpy(B->x, A->x, A->nzmax * sizeof(double));
  memcpy(B->i, A->i, A->nzmax * sizeof(CS_INT));

  size_t size_cpy = -1;
  if(A->nz >= 0)
  {
    size_cpy = A->nzmax;
  }
  else if(A->nz == -1)
  {
    size_cpy = A->n + 1;
  }
  else if(A->nz == -2)
  {
    size_cpy = A->m + 1;
  }

  memcpy(B->p, A->p, size_cpy * sizeof(CS_INT));
}

int CSparseMatrix_max_by_columns(const CSparseMatrix *A, double * max)
{
  CS_INT p, j, n, *Ap ;
  CS_ENTRY *Ax ;
  double s ;

  if(!CS_CSC(A) || !A->x) return (-1) ;               /* check inputs */
  n = A->n ;
  Ap = A->p ;
  Ax = A->x ;
  /* loop over the column */
  for(j = 0 ; j < n ; j++)
  {
    /* loop over the element of the row */
    p = Ap [j] ;
    s = Ax [p] ;
    for(p = Ap [j] ; p < Ap [j+1] ; p++)
      s = CS_MAX(Ax[p], s) ;
    max[j] = s;
  }
  return 0 ;
}
int CSparseMatrix_max_abs_by_columns(const CSparseMatrix *A, double * max)
{
  CS_INT p, j, n, *Ap ;
  CS_ENTRY *Ax ;
  double s ;

  if(!CS_CSC(A) || !A->x) return (-1) ;               /* check inputs */
  n = A->n ;
  Ap = A->p ;
  Ax = A->x ;
  /* loop over the column */
  for(j = 0 ; j < n ; j++)
  {
    /* loop over the element of the row */
    p = Ap [j] ;
    s = fabs(Ax [p]) ;
    for(p = Ap [j] ; p < Ap [j+1] ; p++)
      s = CS_MAX(fabs(Ax[p]), s) ;
    max[j] = s;
  }
  return 0 ;
}
