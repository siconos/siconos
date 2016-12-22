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
  A->diag_indx = NULL;
  A->origin = NS_UNKNOWN;
}

double* NM_sparse_data(NumericsSparseMatrix* A)
{
  switch (A->origin)
  {
  case NS_CSC:
  {
    assert(A->csc);
    return A->csc->x;
    break;
  }
  case NS_CSR:
  {
    assert(A->csr);
    return A->csr->x;
    break;
  }
  case NS_TRIPLET:
  {
    assert(A->triplet);
    return A->triplet->x;
    break;
  }
  default:
    printf("NM_sparse_data :: unknown matrix origin %d", A->origin);
    exit(EXIT_FAILURE);
  }
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
  if (A->diag_indx)
  {
    free(A->diag_indx);
    A->diag_indx = NULL;
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
#elif defined(WITH_SUPERLU)
  p->solver = NS_SUPERLU;
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

CSparseMatrix* NM_csparse_alloc_for_copy(const CSparseMatrix* const m)
{
  assert(m);
  CSparseMatrix* out = NULL;
  if (m->nz >= 0) /* triplet  */
  {
    out = cs_spalloc(m->m, m->n, m->nzmax, 1, 1);
  }
  else if (m->nz == -1) /* csc */
  {
    out = cs_spalloc(m->m, m->n, m->nzmax, 1, 0);
  }
  else if (m->nz == -2) /* csr  */
  {
    out = cs_spalloc(m->n, m->m, m->nzmax, 1, 0);
    out->nz = -2;
    out->m = m->m;
    out->n = m->n;
  }
  else
  {
    fprintf(stderr, "NM_copy :: error unknown type %lld for CSparse matrix\n",  (long long int)m->nz);
    exit(EXIT_FAILURE);
  }

  return out;
}

void NM_sparse_free(void *p)
{
  assert(p);
  NumericsSparseLinearSolverParams* ptr = (NumericsSparseLinearSolverParams*) p;
  cs_lu_factors* cs_lu_A = (cs_lu_factors*)NM_sparse_solver_data(ptr);

  cs_sparse_free(cs_lu_A);

  ptr->solver_data = NULL;
}

size_t NM_sparse_nnz(const CSparseMatrix* const A)
{
  if (A->nz >= 0)
  {
    return (size_t)A->nz;
  }
  else if (A->nz == NS_CS_CSC)
  {
    return (size_t)A->p[A->n];
  }
  else if (A->nz == NS_CS_CSR)
  {
    return (size_t)A->p[A->m];
  }
  else
  {
    fprintf(stderr, "NM_sparse_nnz :: unsupported nz number %d", A->nz);
    exit(EXIT_FAILURE);
  }
}

