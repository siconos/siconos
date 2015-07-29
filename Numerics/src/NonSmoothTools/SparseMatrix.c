#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "SparseMatrix.h"
#include <math.h>
#include <float.h>

//#define VERBOSE_DEBUG

/* add an entry to triplet matrix only if value is not (nearly) null */
int cs_zentry(CSparseMatrix *T, int i, int j, double x)
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

/* from sparse to dense */
double* cs_dense(CSparseMatrix *A)
{
  int m = A->m;
  int n = A->n;
  int *Ap = A->p;
  int *Ai = A->i;
  double *Ax = A->x;
  int nzmax = A->nzmax;
  int nz = A->nz;
  double *r = (double*) malloc(A->m * A->n * sizeof(double));
  for(unsigned int i = 0; i<m*n; ++i)
  {
    r[i] = 0.;
  }
  if(nz < 0)
  {
    for(unsigned int j = 0 ; j < n ; j++)
    {
      printf("    col %d : locations %d to %d\n", j, Ap [j], Ap [j+1]-1);
      for(unsigned int p = Ap [j] ; p < Ap [j+1] ; p++)
      {
        printf("      %d : %g\n", Ai [p], Ax ? Ax [p] : 1) ;
        r[Ai[p] + j*m] = Ax ? Ax [p] : 1;
      }
    }
  }
  else
  {
    printf("triplet: %d-by-%d, nzmax: %d nnz: %d\n", m, n, nzmax, nz) ;
    for(unsigned int p = 0 ; p < nz ; p++)
    {
      printf("    %d %d : %g\n", Ai [p], Ap [p], Ax ? Ax [p] : 1) ;
      r[Ai[p] + Ap[p] * m] = Ax ? Ax[p] : 1;
    }
  }
  return r;
}

/* y = alpha*A*x+beta*y */
int cs_aaxpy(const double alpha, const cs *A, const double *x,
             const double beta, double *y)
{
  int p, j, n, *Ap, *Ai ;
  double *Ax ;
  if(!A || !x || !y) return (0) ;	     /* check inputs */
  n = A->n ;
  Ap = A->p ;
  Ai = A->i ;
  Ax = A->x ;
  for(j = 0 ; j < n ; j++)
  {
    for(p = Ap [j] ; p < Ap [j+1] ; p++)
    {
      y [Ai [p]] *= beta;
      y [Ai [p]] += alpha * Ax [p] * x [j] ;
    }
  }
  return (1) ;
}

void printSparse(const CSparseMatrix* const m)
{
  if (! m)
  {
    fprintf(stderr, "Numerics, CSparseMatrix display failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  assert(m->p);
  assert(m->i);
  assert(m->x);

  printf("Sparse matrix of size %dX%d\n", m->m, m->n);
  printf("maximum number of entries %i\n", m->nzmax);
  printf("storage (nz) %li\n", (long int)(m->nz));

  printf("pointer (p) of size %li= {", (long int)(m->m + 1));
  for (int i = 0 ; i < (m->m + 1); i++)
  {
    printf("%li,  ", (long int)(m->p[i]));
  }
  printf("}\n");
  printf("index (i) of size %li= {", (long int)m->nzmax);
  for (int i = 0 ; i < m->nzmax; i++) printf("%li,  ", (long int)(m->i[i]));
  printf("}\n");

  printf("data (x) of size %li = [\n", (long int)m->nzmax);
  for (int i = 0 ; i < m->m; i++)
  {
    for (int j = m->p[i]; j < m->p[i + 1] ; j++)
    {
      assert(m->nzmax > j);
      printf("%12.8e,  ", m->x[j]);
    }
    printf("\n");
  }
  printf("]\n");


}

void freeSparse(CSparseMatrix* const M)
{

  assert(M);
  if (M->p)
  {
    free(M->p);
    M->p = NULL;
  }
  if (M->i)
  {
    free(M->i);
    M->i = NULL;
  }
  if (M->x)
  {
    free(M->x);
    M->x = NULL;
  }
}
