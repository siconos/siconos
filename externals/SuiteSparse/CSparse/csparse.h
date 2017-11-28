#ifndef _CS_H
#define _CS_H
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stddef.h>
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif
#define CS_VER 3                    /* CSparse Version */
#define CS_SUBVER 1
#define CS_SUBSUB 4
#define CS_DATE "Oct 10, 2014"    /* CSparse release date */
#define CS_COPYRIGHT "Copyright (c) Timothy A. Davis, 2006-2014"

#ifdef MATLAB_MEX_FILE
#undef CS_INT
#define CS_INT mwSignedIndex
#endif
#ifndef CS_INT
#include "SiconosConfig.h"
#include <stdint.h>

#ifdef SICONOS_INT64
#define CS_INT int64_t
#define CS_ID "%ld"
#else
#define CS_INT int32_t
#define CS_ID "%d"
#endif
//#define CS_INT ptrdiff_t
#endif

/* --- primary CSparse routines and data structures ------------------------- */
typedef struct cs_sparse    /* matrix in compressed-column or triplet form */
{
    CS_INT nzmax ;     /* maximum number of entries */
    CS_INT m ;         /* number of rows */
    CS_INT n ;         /* number of columns */
    CS_INT *p ;        /* column pointers (size n+1) or col indices (size nzmax) */
    CS_INT *i ;        /* row indices, size nzmax */
    double *x ;     /* numerical values, size nzmax */
    CS_INT nz ;        /* # of entries in triplet matrix, -1 for compressed-col */
} cs ;

cs *cs_add (const cs *A, const cs *B, double alpha, double beta) ;
CS_INT cs_cholsol (CS_INT order, const cs *A, double *b) ;
cs *cs_compress (const cs *T) ;
CS_INT cs_dupl (cs *A) ;
CS_INT cs_entry (cs *T, CS_INT i, CS_INT j, double x) ;
CS_INT cs_gaxpy (const cs *A, const double *x, double *y) ;
cs *cs_load (FILE *f) ;
CS_INT cs_lusol (CS_INT order, const cs *A, double *b, double tol) ;
cs *cs_multiply (const cs *A, const cs *B) ;
double cs_norm (const cs *A) ;
CS_INT cs_print (const cs *A, CS_INT brief) ;
CS_INT cs_qrsol (CS_INT order, const cs *A, double *b) ;
cs *cs_transpose (const cs *A, CS_INT values) ;
/* utilities */
void *cs_calloc (CS_INT n, size_t size) ;
void *cs_free (void *p) ;
void *cs_realloc (void *p, CS_INT n, size_t size, CS_INT *ok) ;
cs *cs_spalloc (CS_INT m, CS_INT n, CS_INT nzmax, CS_INT values, CS_INT triplet) ;
cs *cs_spfree (cs *A) ;
CS_INT cs_sprealloc (cs *A, CS_INT nzmax) ;
void *cs_malloc (CS_INT n, size_t size) ;

/* --- secondary CSparse routines and data structures ----------------------- */
typedef struct cs_symbolic  /* symbolic Cholesky, LU, or QR analysis */
{
    CS_INT *pinv ;     /* inverse row perm. for QR, fill red. perm for Chol */
    CS_INT *q ;        /* fill-reducing column permutation for LU and QR */
    CS_INT *parent ;   /* elimination tree for Cholesky and QR */
    CS_INT *cp ;       /* column pointers for Cholesky, row counts for QR */
    CS_INT *leftmost ; /* leftmost[i] = min(find(A(i,:))), for QR */
    CS_INT m2 ;        /* # of rows for QR, after adding fictitious rows */
    CS_INT lnz ;    /* # entries in L for LU or Cholesky; in V for QR */
    CS_INT unz ;    /* # entries in U for LU; in R for QR */
} css ;

typedef struct cs_numeric   /* numeric Cholesky, LU, or QR factorization */
{
    cs *L ;         /* L for LU and Cholesky, V for QR */
    cs *U ;         /* U for LU, R for QR, not used for Cholesky */
    CS_INT *pinv ;     /* partial pivoting for LU */
    double *B ;     /* beta [0..n-1] for QR */
} csn ;

typedef struct cs_dmperm_results    /* cs_dmperm or cs_scc output */
{
    CS_INT *p ;        /* size m, row permutation */
    CS_INT *q ;        /* size n, column permutation */
    CS_INT *r ;        /* size nb+1, block k is rows r[k] to r[k+1]-1 in A(p,q) */
    CS_INT *s ;        /* size nb+1, block k is cols s[k] to s[k+1]-1 in A(p,q) */
    CS_INT nb ;        /* # of blocks in fine dmperm decomposition */
    CS_INT rr [5] ;    /* coarse row decomposition */
    CS_INT cc [5] ;    /* coarse column decomposition */
} csd ;

CS_INT *cs_amd (CS_INT order, const cs *A) ;
csn *cs_chol (const cs *A, const css *S) ;
csd *cs_dmperm (const cs *A, CS_INT seed) ;
CS_INT cs_droptol (cs *A, double tol) ;
CS_INT cs_dropzeros (cs *A) ;
CS_INT cs_happly (const cs *V, CS_INT i, double beta, double *x) ;
CS_INT cs_ipvec (const CS_INT *p, const double *b, double *x, CS_INT n) ;
CS_INT cs_lsolve (const cs *L, double *x) ;
CS_INT cs_ltsolve (const cs *L, double *x) ;
csn *cs_lu (const cs *A, const css *S, double tol) ;
cs *cs_permute (const cs *A, const CS_INT *pinv, const CS_INT *q, CS_INT values) ;
CS_INT *cs_pinv (const CS_INT *p, CS_INT n) ;
CS_INT cs_pvec (const CS_INT *p, const double *b, double *x, CS_INT n) ;
csn *cs_qr (const cs *A, const css *S) ;
css *cs_schol (CS_INT order, const cs *A) ;
css *cs_sqr (CS_INT order, const cs *A, CS_INT qr) ;
cs *cs_symperm (const cs *A, const CS_INT *pinv, CS_INT values) ;
CS_INT cs_updown (cs *L, CS_INT sigma, const cs *C, const CS_INT *parent) ;
CS_INT cs_usolve (const cs *U, double *x) ;
CS_INT cs_utsolve (const cs *U, double *x) ;
/* utilities */
css *cs_sfree (css *S) ;
csn *cs_nfree (csn *N) ;
csd *cs_dfree (csd *D) ;

/* --- tertiary CSparse routines -------------------------------------------- */
CS_INT *cs_counts (const cs *A, const CS_INT *parent, const CS_INT *post, CS_INT ata) ;
double cs_cumsum (CS_INT *p, CS_INT *c, CS_INT n) ;
CS_INT cs_dfs (CS_INT j, cs *G, CS_INT top, CS_INT *xi, CS_INT *pstack, const CS_INT *pinv) ;
CS_INT cs_ereach (const cs *A, CS_INT k, const CS_INT *parent, CS_INT *s, CS_INT *w) ;
CS_INT *cs_etree (const cs *A, CS_INT ata) ;
CS_INT cs_fkeep (cs *A, CS_INT (*fkeep) (CS_INT, CS_INT, double, void *), void *other) ;
double cs_house (double *x, double *beta, CS_INT n) ;
CS_INT cs_leaf (CS_INT i, CS_INT j, const CS_INT *first, CS_INT *maxfirst, CS_INT *prevleaf,
    CS_INT *ancestor, CS_INT *jleaf) ;
CS_INT *cs_maxtrans (const cs *A, CS_INT seed) ;
CS_INT *cs_post (const CS_INT *parent, CS_INT n) ;
CS_INT *cs_randperm (CS_INT n, CS_INT seed) ;
CS_INT cs_reach (cs *G, const cs *B, CS_INT k, CS_INT *xi, const CS_INT *pinv) ;
CS_INT cs_scatter (const cs *A, CS_INT j, double beta, CS_INT *w, double *x, CS_INT mark,
    cs *C, CS_INT nz) ;
csd *cs_scc (cs *A) ;
CS_INT cs_spsolve (cs *G, const cs *B, CS_INT k, CS_INT *xi, double *x,
    const CS_INT *pinv, CS_INT lo) ;
CS_INT cs_tdfs (CS_INT j, CS_INT k, CS_INT *head, const CS_INT *next, CS_INT *post,
    CS_INT *stack) ;
/* utilities */
csd *cs_dalloc (CS_INT m, CS_INT n) ;
csd *cs_ddone (csd *D, cs *C, void *w, CS_INT ok) ;
cs *cs_done (cs *C, void *w, void *x, CS_INT ok) ;
CS_INT *cs_idone (CS_INT *p, cs *C, void *w, CS_INT ok) ;
csn *cs_ndone (csn *N, cs *C, void *w, void *x, CS_INT ok) ;

#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define CS_FLIP(i) (-(i)-2)
#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define CS_MARKED(w,j) (w [j] < 0)
#define CS_MARK(w,j) { w [j] = CS_FLIP (w [j]) ; }
#define CS_CSC(A) (A && (A->nz == -1))
#define CS_TRIPLET(A) (A && (A->nz >= 0))
#endif
