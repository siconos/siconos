#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "SparseMatrix.h"
#include <math.h>
#include <float.h>

#include "SiconosCompat.h"
#include "NumericsSparseMatrix.h"

void NM_sparse_null(NumericsSparseMatrix* A)
{
  A->linearSolverParams = NULL;
  A->triplet = NULL;
  A->csc = NULL;
  A->trans_csc = NULL;
  A->csr = NULL;
  A->origin = NS_UNKNOWN;
}

NumericsSparseMatrix* newNumericsSparseMatrix(void)
{
  NumericsSparseMatrix* p = (NumericsSparseMatrix*)
    malloc(sizeof(NumericsSparseMatrix));

  NM_sparse_null(p);
  p->linearSolverParams = newNumericsSparseLinearSolverParams();

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
  if (A->csr)
  {
    cs_spfree(A->csr);
    A->csr = NULL;
  }
  return NULL;
}


NumericsSparseLinearSolverParams* newNumericsSparseLinearSolverParams(void)
{
  NumericsSparseLinearSolverParams* p = (NumericsSparseLinearSolverParams*)
    malloc(sizeof(NumericsSparseLinearSolverParams));

#if defined(WITH_MKL) && 0
  p->solver = NS_PARDISO;
#elif defined(WITH_MUMPS)
  p->solver = NS_MUMPS;
#elif defined(WITH_UMFPACK)
  p->solver = NS_UMFPACK;
#else
  p->solver = NS_CS_LUSOL;
#endif

  p->iparam = NULL;
  p->dparam = NULL;
  p->iWork = NULL;
  p->dWork = NULL;

  p->solver_data = NULL;
  p->solver_free_hook = NULL;

  p->iSize = 0;
  p->dSize = 0;
  p->iWorkSize = 0;
  p->dWorkSize = 0;

  return p;
}

NumericsSparseLinearSolverParams* freeNumericsSparseLinearSolverParams(NumericsSparseLinearSolverParams* p)
{
  /* First free solver_data if some additional information has been given  */
  if (p->solver_free_hook)
    (*p->solver_free_hook)(p);
  p->solver_free_hook = NULL;

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

  free(p);
  return NULL;
}

void NM_sparse_free(void *p)
{
  assert(p);
  NumericsSparseLinearSolverParams* ptr = (NumericsSparseLinearSolverParams*) p;
  cs_lu_factors* cs_lu_A = (cs_lu_factors*)NM_sparse_solver_data(ptr);

  cs_sparse_free(cs_lu_A);

  ptr->solver_data = NULL;
  }
