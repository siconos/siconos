#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "SparseMatrix.h"
#include <math.h>
#include <float.h>

#include "SiconosCompat.h"

//#define VERBOSE_DEBUG

/* add an entry to triplet matrix only if value is not (nearly) null */
csi cs_zentry(CSparseMatrix *T, csi i, csi j, double x)
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
  csi m = A->m;
  csi n = A->n;
  csi *Ap = A->p;
  csi *Ai = A->i;
  double *Ax = A->x;
  csi nzmax = A->nzmax;
  csi nz = A->nz;
  double *r = (double*) malloc(A->m * A->n * sizeof(double));
  for(int i = 0; i<m*n; ++i)
  {
    r[i] = 0.;
  }
  if(nz < 0)
  {
    for(int j = 0 ; j < n ; j++)
    {
      printf("    col %d : locations " SN_PTRDIFF_T_F " to " SN_PTRDIFF_T_F "\n", j, Ap [j], Ap [j+1]-1);
      for(csi p = Ap [j] ; p < Ap [j+1] ; p++)
      {
        printf("      " SN_PTRDIFF_T_F " : %g\n", Ai [p], Ax ? Ax [p] : 1) ;
        r[Ai[p] + j*m] = Ax ? Ax [p] : 1;
      }
    }
  }
  else
  {
    printf("triplet: " SN_PTRDIFF_T_F "-by-" SN_PTRDIFF_T_F ", nzmax: " \
        SN_PTRDIFF_T_F " nnz: " SN_PTRDIFF_T_F "\n", m, n, nzmax, nz) ;
    for(int p = 0 ; p < nz ; p++)
    {
      printf("    " SN_PTRDIFF_T_F " " SN_PTRDIFF_T_F " : %g\n", Ai [p], Ap [p], Ax ? Ax [p] : 1) ;
      r[Ai[p] + Ap[p] * m] = Ax ? Ax[p] : 1;
    }
  }
  return r;
}

/* y = alpha*A*x+beta*y */
int cs_aaxpy(const double alpha, const cs *A, const double *x,
             const double beta, double *y)
{
  csi p, n, *Ap, *Ai ;
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
  return (1) ;
}

CSparseMatrix* cs_spfree_on_stack(CSparseMatrix* A)
{
  if(!A) return (NULL) ;	/* do nothing if A already NULL */
  cs_free(A->p) ;
  cs_free(A->i) ;
  cs_free(A->x) ;
  return NULL;
}

NumericsSparseLinearSolverParams* newNumericsSparseLinearSolverParams(void)
{
  NumericsSparseLinearSolverParams* p = (NumericsSparseLinearSolverParams*)
    malloc(sizeof(NumericsSparseLinearSolverParams));

  p->solver = NS_CS_LUSOL;

  p->iparam = NULL;
  p->dparam = NULL;
  p->iWork = NULL;
  p->dWork = NULL;

#ifdef HAVE_MPI
  p->mpi_com_init = 0;
  p->mpi_com = MPI_COMM_NULL;
#endif

  p->solver_data = NULL;

  p->iSize = 0;
  p->dSize = 0;
  p->iWorkSize = 0;
  p->dWorkSize = 0;

  return p;
}

NumericsSparseLinearSolverParams* freeNumericsSparseLinearSolverParams(NumericsSparseLinearSolverParams* p)
{
  if (p->iparam)
  {
    assert(p->iSize>0);
    free(p->iparam);
    p->iparam = NULL;
  }
  if (p->dparam)
  {
    assert(p->dSize>0);
    free(p->dparam);
    p->dparam = NULL;
  }
  if (p->iWork)
  {
    assert(p->iWorkSize>0);
    free(p->iWork);
    p->iWork = NULL;
  }
  if (p->dWork)
  {
    assert(p->dWorkSize>0);
    free(p->dWork);
    p->dWork = NULL;
  }
  if (p->solver_data)
  {
    free(p->solver_data);
    p->solver_data = NULL;
  }

#ifdef HAVE_MPI
  if (p->mpi_com)
  {
    /* MPI_Finalize called only if initialization has been done for
     * this matrix */
    if (p->mpi_com_init)
    {
      MPI_Finalize();
    }
  }
#endif

  free(p);
  return NULL;
}

NumericsSparseMatrix* newNumericsSparseMatrix(void)
{
  NumericsSparseMatrix* p = (NumericsSparseMatrix*)
    malloc(sizeof(NumericsSparseMatrix));

  p->linearSolverParams = newNumericsSparseLinearSolverParams();

  p->triplet = NULL;
  p->csc = NULL;
  p->trans_csc = NULL;

  return p;
}

NumericsSparseMatrix* freeNumericsSparseMatrix(NumericsSparseMatrix* A)
{
  if (A->linearSolverParams)
  {
    freeNumericsSparseLinearSolverParams(A->linearSolverParams);
    A->linearSolverParams = NULL;
  }
  if (A->triplet)
  {
    cs_spfree(A->triplet);
    A->triplet = NULL;
  }
  if (A->csc)
  {
    cs_spfree(A->csc);
    A->csc = NULL;
  }
  if (A->trans_csc)
  {
    cs_spfree(A->trans_csc);
    A->trans_csc = NULL;
  }
  return NULL;
}

