
/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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
#include <assert.h>                   // for assert
#include <float.h>                    // for DBL_EPSILON
#include <math.h>                     // for fabs, fmax, NAN
#include <stdint.h>                   // for SIZE_MAX
#include <stdio.h>                    // for printf, fprintf, size_t, fscanf
#include <stdlib.h>                   // for exit, malloc, free, EXIT_FAILURE
#include <string.h>                   // for memcpy, memset
#include "CSparseMatrix_internal.h"            // for CSparseMatrix, CS_INT, cs_dl_sp...
#include "NM_MPI.h"                   // for NM_MPI_copy
#include "NM_MUMPS.h"                 // for NM_MUMPS_copy
#include "NM_conversions.h"           // for NM_csc_to_csr, NM_csc_to_triplet
#include "NumericsFwd.h"              // for NumericsMatrix, NumericsSparseM...
#include "NumericsMatrix.h"           // for NumericsMatrix, NumericsMatrixI...
#include "NumericsMatrix_internal.h"  // for NM_internalData_free
#include "NumericsSparseMatrix.h"     // for NumericsSparseMatrix, NSM_new
#include "SiconosCompat.h"            // for SN_SIZE_T_F
#include "SiconosBlas.h"              // for cblas_ddot, cblas_dgemv, CblasN...
#include "SiconosLapack.h"            // for lapack_int, DGESV, DGETRF, DGETRS, LA_NOTRANS
#include "SparseBlockMatrix.h"        // for SparseBlockStructuredMatrix
/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "siconos_debug.h"                    // for DEBUG_EXPR, DEBUG_BEGIN, DEBUG_...
#include "numerics_verbose.h"         // for numerics_error, numerics_printf...
#include "sanitizer.h"                // for cblas_dcopy_msan
#include "NumericsVector.h"           // for NV_max

#ifdef WITH_OPENSSL
#include <openssl/sha.h>
#endif

#ifdef WITH_MKL_SPBLAS
#include "MKL_common.h"
#include "NM_MKL_spblas.h"
#endif

#ifdef WITH_MA57
#include "lbl.h"
#include "NM_MA57.h"
#endif


#ifdef __cplusplus
#undef restrict
#include <sys/cdefs.h>                // for __restrict
#define restrict __restrict
#endif

void NM_null(NumericsMatrix* A)
{

  A->storageType= NM_UNKNOWN;
  A->size0 =-1;
  A->size1 =-1;
  A->matrix0 = NULL;
  A->matrix1 = NULL;
  A->matrix2 = NULL;
  A->internalData = NULL;
  NDV_reset(&(A->version));
  A->destructible = A; /* by default, the destructible matrix is itself */
}

void NM_internalData_new(NumericsMatrix* M)
{
  M->internalData = (NumericsMatrixInternalData *)malloc(sizeof(NumericsMatrixInternalData));
  M->internalData->iWork = NULL;
  M->internalData->iWorkSize = 0;
  M->internalData->dWork = NULL;
  M->internalData->dWorkSize = 0;
  M->internalData->isInversed = false ;
  M->internalData->isLUfactorized = false ;
  M->internalData->isCholeskyfactorized = false ;
  M->internalData->isLDLTfactorized = false ;
#ifdef SICONOS_HAS_MPI
  M->internalData->mpi_comm = MPI_COMM_NULL;
#endif
#ifdef WITH_OPENSSL
  M->internalData->values_sha1_count = 0;
#endif
}


static version_t NM_version(const NumericsMatrix* M, NM_types id)
{
  switch (id)
  {
  case NM_DENSE:
  {
    return NDV_value(&(M->version));
  }
  case NM_SPARSE_BLOCK:
  {
    if (M->matrix1)
    {
      return NDV_value(&(M->matrix1->version));
    }
    else
    {
      return 0;
    }
  }
  case NM_SPARSE:
  {
    if (M->matrix2)
    {
      return NSM_max_version(M->matrix2);
    }
    else
    {
      return 0;
    }
  }
  default:
    numerics_error("NM_version", "unknown id");
    return 0;
  }
  assert (false);
}

void NM_reset_version(NumericsMatrix* M, NM_types id)
{
  switch (id)
  {
  case NM_DENSE:
  {
    NDV_reset(&(M->version));
    break;
  }
  case NM_SPARSE_BLOCK:
  {
  if (M->matrix1)
      NDV_reset(&(M->matrix1->version));
    break;
  }
  case NM_SPARSE:
  {
    if (M->matrix2)
    {
      NSM_reset_versions(M->matrix2);
    }
    break;
  }
  default: numerics_error("NM_reset_version", "unknown id");
  }
}

void NM_reset_versions(NumericsMatrix* M)
{
  NM_reset_version(M, NM_DENSE);
  NM_reset_version(M, NM_SPARSE_BLOCK);
  NM_reset_version(M, NM_SPARSE);
}

static void NM_set_version(NumericsMatrix* M, NM_types id, version_t value)
{
  switch (id)
  {
  case NM_DENSE:
  {
    NDV_set_value(&(M->version), value);
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    NDV_set_value(&(M->matrix1->version), value);
    break;
  }
  case NM_SPARSE:
  {
    numerics_error("NM_set_version", "cannot set version of sparse matrix, use NSM_set_version");
    break;
  }
  default: numerics_error("NM_set_version", "unknown id");
  }
}

/* internal compare function */
static NM_types nm_max(const NumericsMatrix* M,
             NM_types type1,
             NM_types type2)
{
  return NM_version(M, type1) > NM_version(M, type2) ?
    type1 : type2;
}

static NM_types NM_latest_id(const NumericsMatrix* M)
{
  return nm_max(M, nm_max(M, NM_DENSE, NM_SPARSE_BLOCK), NM_SPARSE);
}


static version_t NM_max_version(const NumericsMatrix* M)
{
  return NM_version(M, NM_latest_id(M));
}

static void NM_inc_version(NumericsMatrix* M, NM_types id)
{
  version_t new_version = NM_max_version(M) + 1;

  switch(id)
  {
  case NM_DENSE:
  case NM_SPARSE_BLOCK:
  {
    NM_set_version(M, id, new_version);
    break;
  }
  case NM_SPARSE:
  {
    numerics_error("NM_inc_version",
                   "cannot increment version of sparse matrix, use NSM_inc_version instead");
    break;
  }
  default: numerics_error("NM_inc_version", "unknown storage");
  }
}


void NM_prod_mv_3x3(int sizeX, int sizeY, NumericsMatrix* A,
                    double* const x, double* y)
{

  assert(A);
  assert(x);
  assert(y);
  assert(A->size0 == sizeY);
  assert(A->size1 == sizeX);
  double alpha=1;
  double beta=1;

  NM_types storage = A->storageType;

  /* double* storage */

  assert (NM_version(A, A->storageType) ==
          NM_version(A, NM_latest_id(A)));
  switch(storage)
  {
  case NM_DENSE:
    cblas_dgemv(CblasColMajor, CblasNoTrans, sizeY, sizeX, alpha, A->matrix0, sizeY, x, 1, beta, y, 1);
    break;
  /* SparseBlock storage */
  case NM_SPARSE_BLOCK:
    SBM_gemv_3x3(sizeX, sizeY, A->matrix1, x, y);
    break;
  /* coordinate */
  case NM_SPARSE:
    CSparseMatrix_aaxpby(alpha, NM_csc(A), x, beta, y);
    break;

  default:
    fprintf(stderr, "Numerics, NumericsMatrix, product matrix - vector prod(A,x,y) failed, unknown storage type for A.\n");
    exit(EXIT_FAILURE);
  }
}

void NM_row_prod(int sizeX, int sizeY, int currentRowNumber, const NumericsMatrix* A, const double* const x, double* y, int init)
{

  assert(A);
  assert(x);
  assert(y);
  assert(A->size0 >= sizeY);
  assert(A->size1 == sizeX);

  NM_types storage = A->storageType;

  /* double* storage */
  if(storage == NM_DENSE)
  {
    int incx = A->size0, incy = 1;
    double* mat = A->matrix0;
    if(init == 0)  /* y += subAx */
    {
      for(int row = 0; row < sizeY; row++)
        y[row] += cblas_ddot(sizeX, &mat[currentRowNumber + row], incx, x, incy);
    }
    else
    {
      for(int row = 0; row < sizeY; row++)
        y[row] = cblas_ddot(sizeX, &mat[currentRowNumber + row], incx, x, incy);
    }

  }
  /* SparseBlock storage */
  else if(storage == NM_SPARSE_BLOCK)
    SBM_row_prod(sizeX, sizeY, currentRowNumber, A->matrix1, x, y, init);
  else
  {
    fprintf(stderr, "Numerics, NumericsMatrix, product matrix - vector NM_row_prod(A,x,y) failed, unknown storage type for A.\n");
    exit(EXIT_FAILURE);
  }

}

void NM_row_prod_no_diag(size_t sizeX, size_t sizeY, int block_start, size_t row_start, NumericsMatrix* A, double* restrict x, double* restrict y, double* restrict xsave, bool init)
{
  assert(A);
  assert(x);
  assert(y);
  assert((size_t) A->size0 >= sizeY);
  assert((size_t)A->size1 == sizeX);

  switch(A->storageType)
  {
  case NM_DENSE:
  {
    double * xSave;
    if(xsave) xSave = xsave;
    else xSave = (double*) malloc(sizeY * sizeof(double));

    memcpy(xSave, &x[row_start], sizeY * sizeof(double));
    memset(&x[row_start], 0, sizeY * sizeof(double));

    double * M = A->matrix0 + row_start;
    int incx = A->size0;
    int incy = 1;
    if(init)
    {
      memset(y, 0, sizeY*sizeof(double));
    }
    for(size_t i = 0; i < sizeY; ++i)
    {
      y[i] += cblas_ddot(A->size0, M, incx, x, incy);
      ++M;
    }

    memcpy(&x[row_start], xSave, sizeY*sizeof(double));

    if(!xsave)
    {
      free(xSave);
    }

    break;
  }
  case NM_SPARSE_BLOCK:
  {
    SBM_row_prod_no_diag(sizeX, sizeY, block_start, A->matrix1, x, y, init);
    break;
  }
  case NM_SPARSE:
  {
    double * xSave;
    if(xsave) xSave = xsave;
    else xSave = (double*) malloc(sizeY * sizeof(double));

    memcpy(xSave, &x[row_start], sizeY * sizeof(double));
    memset(&x[row_start], 0, sizeY * sizeof(double));

    if(init)
    {
      memset(y, 0, sizeY*sizeof(double));
    }

    CSparseMatrix* M;
    if(A->matrix2->origin == NSM_CSR)
    {
      M = NM_csr(A);
    }
    else
    {
      M = NM_csc_trans(A);
    }

    CS_INT* Mp = M->p;
    CS_INT* Mi = M->i;
    double* Mx = M->x;

    for(size_t i = 0, j = row_start; i < sizeY; ++i, ++j)
    {
      for(CS_INT p = Mp[j]; p < Mp[j+1]; ++p)
      {
        y[i] += Mx[p] * x[Mi[p]];
      }
    }

    memcpy(&x[row_start], xSave, sizeY*sizeof(double));

    if(!xsave)
    {
      free(xSave);
    }
    break;
  }
  default:
  {
    fprintf(stderr, "Numerics, NumericsMatrix, product matrix - vector NM_row_prod_no_diag(A,x,y) failed, unknown storage type for A.\n");
    exit(EXIT_FAILURE);
  }
  }
}


void NM_row_prod_no_diag3(size_t sizeX, int block_start, size_t row_start, NumericsMatrix* A, double* x, double* y, bool init)
{
  assert(A);
  assert(x);
  assert(y);
  assert((size_t)A->size0 >= 3);
  assert((size_t)A->size1 == sizeX);

  switch(A->storageType)
  {
  case NM_DENSE:
  {
    if(init)
    {
      y[0] = 0.;
      y[1] = 0.;
      y[2] = 0.;
    }
    double* M = A->matrix0;
    assert(M);
    int incx = sizeX, incy = 1;
    size_t in = row_start, it = row_start + 1, is = row_start + 2;
    double rin = x[in] ;
    double rit = x[it] ;
    double ris = x[is] ;
    x[in] = 0.;
    x[it] = 0.;
    x[is] = 0.;
    y[0] += cblas_ddot(sizeX, &M[in], incx, x, incy);
    y[1] += cblas_ddot(sizeX, &M[it], incx, x, incy);
    y[2] += cblas_ddot(sizeX, &M[is], incx, x, incy);
    x[in] = rin;
    x[it] = rit;
    x[is] = ris;
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    /* qLocal += rowMB * x
     * with rowMB the row of blocks of MGlobal which corresponds
     * to the current contact
     */
    SBM_row_prod_no_diag_3x3(sizeX, 3, block_start, A->matrix1, x, y);
    break;
  }
  case NM_SPARSE:
  {
    if(init)
    {
      y[0] = 0.;
      y[1] = 0.;
      y[2] = 0.;
    }

    size_t in = row_start, it = row_start + 1, is = row_start + 2;
    double rin = x[in] ;
    double rit = x[it] ;
    double ris = x[is] ;
    x[in] = 0.;
    x[it] = 0.;
    x[is] = 0.;

    CSparseMatrix* M;
    if(A->matrix2->origin == NSM_CSR)
    {
      M = NM_csr(A);
    }
    else
    {
      M = NM_csc_trans(A);
    }

    CS_INT* Mp = M->p;
    CS_INT* Mi = M->i;
    double* Mx = M->x;

    for(size_t i = 0, j = row_start; i < 3; ++i, ++j)
    {
      for(CS_INT p = Mp[j]; p < Mp[j+1]; ++p)
      {
        y[i] += Mx[p] * x[Mi[p]];
      }
    }

    x[in] = rin;
    x[it] = rit;
    x[is] = ris;

    break;
  }
  default:
  {
    fprintf(stderr, "NM_row_prod_no_diag3 :: unknown matrix storage %d", A->storageType);
    exit(EXIT_FAILURE);
  }
  }
}
void NM_row_prod_no_diag2(size_t sizeX, int block_start, size_t row_start, NumericsMatrix* A, double* x, double* y, bool init)
{
  assert(A);
  assert(x);
  assert(y);
  assert((size_t)A->size0 >= 2);
  assert((size_t)A->size1 == sizeX);

  switch(A->storageType)
  {
  case NM_DENSE:
  {
    if(init)
    {
      y[0] = 0.;
      y[1] = 0.;
    }
    double* M = A->matrix0;
    assert(M);
    int incx = sizeX, incy = 1;
    size_t in = row_start, it = row_start + 1;
    double rin = x[in] ;
    double rit = x[it] ;
    x[in] = 0.;
    x[it] = 0.;
    y[0] += cblas_ddot(sizeX, &M[in], incx, x, incy);
    y[1] += cblas_ddot(sizeX, &M[it], incx, x, incy);
    x[in] = rin;
    x[it] = rit;
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    /* qLocal += rowMB * x
     * with rowMB the row of blocks of MGlobal which corresponds
     * to the current contact
     */
    SBM_row_prod_no_diag_2x2(sizeX, 2, block_start, A->matrix1, x, y);
    break;
  }
  case NM_SPARSE:
  {
    if(init)
    {
      y[0] = 0.;
      y[1] = 0.;
    }

    size_t in = row_start, it = row_start + 1;
    double rin = x[in] ;
    double rit = x[it] ;
    x[in] = 0.;
    x[it] = 0.;

    CSparseMatrix* M;
    if(A->matrix2->origin == NSM_CSR)
    {
      M = NM_csr(A);
    }
    else
    {
      M = NM_csc_trans(A);
    }

    CS_INT* Mp = M->p;
    CS_INT* Mi = M->i;
    double* Mx = M->x;

    for(size_t i = 0, j = row_start; i < 2; ++i, ++j)
    {
      for(CS_INT p = Mp[j]; p < Mp[j+1]; ++p)
      {
        y[i] += Mx[p] * x[Mi[p]];
      }
    }

    x[in] = rin;
    x[it] = rit;

    break;
  }
  default:
  {
    fprintf(stderr, "NM_row_prod_no_diag2 :: unknown matrix storage %d", A->storageType);
    exit(EXIT_FAILURE);
  }
  }
}

void NM_row_prod_no_diag1x1(size_t sizeX, int block_start, size_t row_start, NumericsMatrix* A, double* x, double* y, bool init)
{
  assert(A);
  assert(x);
  assert(y);
  assert((size_t)A->size0 >= 1);
  assert((size_t)A->size1 == sizeX);

  switch(A->storageType)
  {
  case NM_DENSE:
  {
    if(init)
    {
      y[0] = 0.;
    }
    double* M = A->matrix0;
    assert(M);
    int incx = sizeX, incy = 1;
    double xin = x[row_start] ;

    x[row_start] = 0.;
    y[0] += cblas_ddot(sizeX, &M[row_start], incx, x, incy);
    x[row_start] = xin;
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    /* qLocal += rowMB * x
     * with rowMB the row of blocks of MGlobal which corresponds
     * to the current contact
     */
    SBM_row_prod_no_diag_1x1(sizeX, 1, block_start, A->matrix1, x, y);
    break;
  }
  case NM_SPARSE:
  {
    if(init)
    {
      y[0] = 0.;
    }

    double rin = x[row_start] ;
    x[row_start] = 0.;

    CSparseMatrix* M;
    if(A->matrix2->origin == NSM_CSR)
    {
      M = NM_csr(A);
    }
    else
    {
      M = NM_csc_trans(A);
    }

    CS_INT* Mp = M->p;
    CS_INT* Mi = M->i;
    double* Mx = M->x;

    for(CS_INT p = Mp[row_start]; p < Mp[row_start+1]; ++p)
    {
      y[0] += Mx[p] * x[Mi[p]];
    }
    x[row_start] = rin;
    break;
  }
  default:
  {
    fprintf(stderr, "NM_row_prod_no_diag3 :: unknown matrix storage %d", A->storageType);
    exit(EXIT_FAILURE);
  }
  }
}
void NM_internalData_free(NumericsMatrix* m)
{
  assert(m && "NM_internalData_free, m == NULL");
  if(m->internalData)
  {
    if(m->internalData->iWork)
    {
      assert(m->internalData->iWorkSize > 0);
      free(m->internalData->iWork);
    }
    m->internalData->iWork = NULL;
    if(m->internalData->dWork)
    {
      assert(m->internalData->dWorkSize > 0);
      free(m->internalData->dWork);
    }
    m->internalData->dWork = NULL;
    free(m->internalData);
    m->internalData = NULL;
  }
}
void NM_internalData_copy(const NumericsMatrix* const A, NumericsMatrix* B)
{

  if(!A->internalData) /* no internal data in A */
  {
    if(B->internalData)
    {
      NM_internalData_free(B);
    }
  }
  else
  {
    if(!B->internalData)
    {
      NM_internalData_new(B);
    }


    size_t sizeof_elt = A->internalData->sizeof_elt;

    if(A->internalData->iWork)
    {

      size_t size = A->internalData->iWorkSize / sizeof_elt;


      if(! B->internalData->iWork)
      {
        B->internalData->iWork = malloc(size*sizeof_elt);
        B->internalData->iWorkSize = A->internalData->iWorkSize;
      }
      else
      {
        if(A->internalData->iWorkSize != B->internalData->iWorkSize)
        {
          B->internalData->iWorkSize = A->internalData->iWorkSize;
          B->internalData->iWork = realloc(B->internalData->iWork,size*sizeof_elt);
        }
      }

      memcpy(B->internalData->iWork, A->internalData->iWork, A->internalData->iWorkSize);
    }
    else
    {
      if(B->internalData->iWork)
        free(B->internalData->iWork);
      B->internalData->iWorkSize=0;
    }



    if(A->internalData->dWork)
    {
      if(! B->internalData->dWork)
      {
        B->internalData->dWork = (double*)malloc(A->internalData->dWorkSize*sizeof(double));
        B->internalData->dWorkSize = A->internalData->dWorkSize;
      }
      else
      {
        if(A->internalData->dWorkSize != B->internalData->dWorkSize)
        {
          B->internalData->dWorkSize = A->internalData->dWorkSize;
          B->internalData->dWork = realloc(B->internalData->dWork,A->internalData->dWorkSize*sizeof(double));
        }
      }

      for(unsigned i =0; i < A->internalData->dWorkSize; i++)
        B->internalData->dWork[i] = A->internalData->dWork[i];
    }
    else
    {
      if(B->internalData->dWork)
        free(B->internalData->dWork);
      B->internalData->dWorkSize=0;
    }
    B->internalData->isLUfactorized = A->internalData->isLUfactorized;
    B->internalData->isCholeskyfactorized = A->internalData->isCholeskyfactorized;
    B->internalData->isLDLTfactorized = A->internalData->isLDLTfactorized;
    B->internalData->isInversed = A->internalData->isInversed;
  }

}
void NM_clear(NumericsMatrix* m)
{
  assert(m && "NM_clear, m == NULL");

  NM_clearDense(m);
  NM_clearSparseBlock(m);
  NM_clearSparse(m);

  NM_internalData_free(m);

  /* restore the destructible pointer */
  if (!NM_destructible(m))
  {
    NM_clear(m->destructible);
    m->destructible = m;
  }

}

NumericsMatrix*  NM_free(NumericsMatrix* m)
{
  assert(m && "NM_free, m == NULL");

  NM_clear(m);
  free(m);
  return NULL;

}

void  NM_clear_not_dense(NumericsMatrix* m)
{
  DEBUG_BEGIN("NM_clear_not_dense(NumericsMatrix* m)\n");
  assert(m && "NM_clear_not_dense, m == NULL");

  //NM_clearDense(m);
  NM_clearSparseBlock(m);
  NM_clearSparse(m);

  NM_internalData_free(m);
  /* restore the destructible pointer */
  if (!NM_destructible(m))
  {
    NM_clear(m->destructible);
    m->destructible = m;
  }
  DEBUG_END("NM_clear_not_dense(NumericsMatrix* m)\n");
}
NumericsMatrix*  NM_free_not_dense(NumericsMatrix* m)
{
  assert(m && "NM_free_not_dense, m == NULL");
  NM_clear_not_dense(m);
  free(m);
  return NULL;
}

void  NM_clear_not_SBM(NumericsMatrix* m)
{
  assert(m && "NM_clear_not_SBM, m == NULL");

  NM_clearDense(m);
/*  NM_clearSparseBlock(m); */
  NM_clearSparse(m);
  NM_internalData_free(m);
  /* restore the destructible pointer */
  if (!NM_destructible(m))
  {
    NM_clear(m->destructible);
    m->destructible = m;
  }
}

NumericsMatrix*  NM_free_not_SBM(NumericsMatrix* m)
{
  assert(m && "NM_free_not_SBM, m == NULL");
  NM_clear_not_SBM(m);
  free(m);
  return NULL;
}









void NM_clear_other_storages(NumericsMatrix* M, NM_types storageType)
{
  assert(M && "NM_clear, M == NULL");

  switch(storageType)
  {
  case NM_DENSE:
  {
    NM_clearSparseBlock(M);
    NM_clearSparse(M);
    NM_internalData_free(M);
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    NM_clearDense(M);
    NM_clearSparse(M);
    NM_internalData_free(M);
    break;
  }
  case NM_SPARSE:
  {
    NM_clearDense(M);
    NM_clearSparseBlock(M);
    NM_internalData_free(M);
    break;
  }
  default:
    numerics_error("NM_clear_other_storages ","unknown storageType %d for matrix\n", M->storageType);
  }
}


void NM_dense_display_matlab(double * m, int nRow, int nCol, int lDim)
{
  if(m)
  {
    int lin, col;
    if(lDim == 0)
      lDim = nRow;
    printf("Matrix of size\t%d\t x \t%d =\n[", nRow, nCol);
    if(nRow == 0)
    {
      printf("]\n");
    }
    if(nCol == 0)
    {
      printf("]\n");
    }

    for(lin = 0; lin < nRow; lin++)
    {
      for(col = 0; col < nCol; col++)
      {
        printf(" %.15e", m[lin + col * lDim]);
        if(col != nCol - 1)
          printf(",");
      }
      if(lin != nRow - 1)
        printf(";\n");
      else
        printf("]\n");
    }
  }
  else
    printf("Matrix : NULL\n");

}

void NM_dense_display(double * m, int nRow, int nCol, int lDim)
{
  if(m)
  {
    int lin, col;
    if(lDim == 0)
      lDim = nRow;
    printf("Matrix of size\t%d\t x \t%d =\n[", nRow, nCol);
    if(nRow == 0)
    {
      printf("]\n");
    }
    if(nCol == 0)
    {
      printf("]\n");
    }

    for(lin = 0; lin < nRow; lin++)
    {
      printf("[");
      for(col = 0; col < nCol; col++)
      {
        printf(" %.15e", m[lin + col * lDim]);
        if(col != nCol - 1)
          printf(",");
      }
      if(lin != nRow - 1)
        printf("],\n");
      else
        printf("]\t ]\n");
    }
  }
  else
    printf("Matrix : NULL\n");

}

void NM_zentry(NumericsMatrix* M, int i, int j, double val, double threshold)
{
  int insertion=0;
  switch(M->storageType)
  {
  case NM_DENSE:
  {
    // column major
    if (fabs(val) >= threshold)
    {
      M->matrix0[i+j*M->size0] = val;
      insertion=1;
      NM_inc_version(M, NM_DENSE);
    }
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    if (fabs(val) >= threshold)
    {
      /* version is incremented in SBM_zentry */
      CHECK_RETURN(SBM_entry(M->matrix1, i, j, val));
      insertion=1;
    }
    break;
  }
  case NM_SPARSE:
  {
    assert(M->matrix2);
    switch(M->matrix2->origin)
    {
    case NSM_TRIPLET:
    {
      assert(M->matrix2->triplet);
      insertion = 1 - CSparseMatrix_zentry(M->matrix2->triplet, i, j, val, threshold);
      NSM_inc_version(M->matrix2, NSM_TRIPLET);
      break;
    }
    case NSM_HALF_TRIPLET:
    {
      assert(M->matrix2->half_triplet);
      insertion = 1 - CSparseMatrix_symmetric_zentry(M->matrix2->triplet, i, j, val, threshold);
      NSM_inc_version(M->matrix2, NSM_HALF_TRIPLET);
      break;
    }
    case NSM_CSC:
    {
      assert(M->matrix2->csc);
      insertion = 1 -CSparseMatrix_zentry(NM_triplet(M), i, j, val, threshold);
      NSM_inc_version(M->matrix2, NSM_TRIPLET);

      M->matrix2->origin= NSM_TRIPLET;
      NM_clearCSC(M);
      NM_csc(M);
      M->matrix2->origin= NSM_CSC;
      NM_clearTriplet(M);
      NM_clearHalfTriplet(M);
      break;
    }
    default:
    {
      numerics_error("NM_zentry","unknown origin %d for sparse matrix\n", M->matrix2->origin);
      break;
    }
    }
    break;
  }
  default:
    numerics_error("NM_zentry  ","unknown storageType %d for matrix\n", M->storageType);
  }
  if(insertion)
  {
    if(i>M->size0-1)
    {
      M->size0 = i+1;
    }
    if(j>M->size1-1)
    {
      M->size1 = j+1;
    }
  }
}

void NM_entry(NumericsMatrix* M, int i, int j, double val)
{
  switch(M->storageType)
  {
  case NM_DENSE:
  {
    // column major
    M->matrix0[i+j*M->size0] = val;
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    CHECK_RETURN(SBM_entry(M->matrix1, i, j, val));
    break;
  }
  case NM_SPARSE:
  {
    assert(M->matrix2);
    switch(M->matrix2->origin)
    {
    case NSM_TRIPLET:
    {
      assert(M->matrix2->triplet);
      CHECK_RETURN(CSparseMatrix_entry(M->matrix2->triplet, i, j, val));
      break;
    }
    case NSM_HALF_TRIPLET:
    {
      assert(M->matrix2->half_triplet);
      CHECK_RETURN(CSparseMatrix_symmetric_entry(M->matrix2->triplet, i, j, val));
      break;
    }
    case NSM_CSC:
    {
      assert(M->matrix2->csc);
      CHECK_RETURN(CSparseMatrix_entry(NM_triplet(M), i, j, val));
      M->matrix2->origin= NSM_TRIPLET;
      NM_clearCSC(M);
      NM_csc(M);
      M->matrix2->origin= NSM_CSC;
      NM_clearTriplet(M);
      NM_clearHalfTriplet(M);
      break;
    }
    default:
    {
      numerics_error("NM_entry","unknown origin %d for sparse matrix\n", M->matrix2->origin);
      break;
    }
    }
    break;
  }
  default:
    numerics_error("NM_entry  ","unknown storageType %d for matrix\n", M->storageType);
  }

  if(i>M->size0-1)
  {
    M->size0 = i+1;
  }
  if(j>M->size1-1)
  {
    M->size1 = j+1;
  }
}



double NM_get_value(const NumericsMatrix* const M, int i, int j)
{
  assert(M);

  if((i + 1 > M->size0) || (j + 1 > M->size1))
  {
    fprintf(stderr, "NM_get_value :: out of range \n");
    exit(EXIT_FAILURE);
  }
  switch(M->storageType)
  {
  case NM_DENSE:
  {
    assert(M->matrix0);
    return M->matrix0[i+j*M->size0];
  }
  case NM_SPARSE_BLOCK:
  {
    assert(M->matrix1);
    return SBM_get_value(M->matrix1,i,j);
  }
  case NM_SPARSE:
  {
    assert(M->matrix2);
    switch(M->matrix2->origin)
    {
    case NSM_TRIPLET:
    {
      assert(M->matrix2->triplet);
      assert(i <M->matrix2->triplet->m);
      assert(j <M->matrix2->triplet->n);
      CS_INT * Mi =   M->matrix2->triplet->i;
      CS_INT * Mp =   M->matrix2->triplet->p;
      double * Mx =   M->matrix2->triplet->x;

      for(int idx = 0 ; idx < M->matrix2->triplet->nz  ; idx++)
      {
        if(Mi[idx] ==i && Mp[idx] == j)
        {
          return Mx[idx];
        }
      }
      return 0.0;
      break;
    }
    case NSM_HALF_TRIPLET:
    {
      assert(M->matrix2->triplet);
      assert(i <M->matrix2->triplet->m);
      assert(j <M->matrix2->triplet->n);
      CS_INT * Mi =   M->matrix2->triplet->i;
      CS_INT * Mp =   M->matrix2->triplet->p;
      double * Mx =   M->matrix2->triplet->x;

      for(int idx = 0 ; idx < M->matrix2->triplet->nz  ; idx++)
      {
        if((Mi[idx] ==i && Mp[idx] == j)||(Mi[idx] ==j && Mp[idx] == i))
        {
          return Mx[idx];
        }
      }
      return 0.0;
      break;
    }
    case NSM_CSC:
    {
      assert(M->matrix2->csc);
      CS_INT * Mi =   M->matrix2->csc->i;
      CS_INT * Mp =   M->matrix2->csc->p;
      double * Mx =   M->matrix2->csc->x;

      for(CS_INT row = Mp[j]; row < Mp[j+1] ; row++)
      {
        if(i == Mi[row])
          return  Mx[row];
      }
      return 0.0;
      break;
    }
    default:
    {
      fprintf(stderr, "NM_get_value ::  unknown origin %d for sparse matrix\n", M->matrix2->origin);
    }
    }
    break;
  }
  default:
    fprintf(stderr, "NM_get_value ::  unknown matrix storage = %d\n", M->storageType);
  }

  return 0.0;

}

bool NM_equal(NumericsMatrix* A, NumericsMatrix* B)
{
  return NM_compare(A, B, DBL_EPSILON*2);
}

bool NM_compare(NumericsMatrix* A, NumericsMatrix* B, double tol)
{
  DEBUG_BEGIN("NM_compare(NumericsMatrix* A, NumericsMatrix* B, double tol)\n");
  assert(A);
  assert(B);
  if(A->size0 != B->size0)
  {
    return false;
  }
  if(A->size1 != B->size1)
  {
    return false;
  }
  for(int i =0; i< A->size0 ; i++)
  {
    for(int j =0; j< A->size1 ; j++)
    {
      /* DEBUG_PRINTF("error %i %i = %e\n",i,j, fabs(NM_get_value(A, i, j) - NM_get_value(B, i, j))); */
      if(fabs(NM_get_value(A, i, j) - NM_get_value(B, i, j)) >= tol)
      {

        DEBUG_PRINTF("A(%i,%i) = %e\t, B(%i,%i) = %e\t,  error = %e\n",
                     i,j, NM_get_value(A, i, j),
                     i,j, NM_get_value(B, i, j),
                     fabs(NM_get_value(A, i, j) - NM_get_value(B, i, j)));

        printf("A(%i,%i) = %e\t, B(%i,%i) = %e\t,  error = %e\n",
                     i,j, NM_get_value(A, i, j),
                     i,j, NM_get_value(B, i, j),
                     fabs(NM_get_value(A, i, j) - NM_get_value(B, i, j)));
        return false;
      }
    }
  }
  DEBUG_END("NM_compare(NumericsMatrix* A, NumericsMatrix* B, double tol)\n");
  return true;
}


void NM_vector_display(double * m, int nRow)
{
  int lin;
  printf("vector of size\t%d\t =\n[", nRow);
  if(nRow == 0)
  {
    printf("]\n");
  }
  for(lin = 0; lin < nRow; lin++)
  {
    printf(" %.15e", m[lin]);
    if(lin != nRow - 1)
      printf(", ");
    else
      printf("]\n");
  }

}

void NM_display_storageType(const NumericsMatrix* const m)
{
  if(! m)
  {
    fprintf(stderr, "Numerics, NumericsMatrix display failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  printf("========== Numerics Matrix\n");

  printf("========== size0 = %i, size1 = %i\n", m->size0, m->size1);

  switch(m->storageType)
  {
  case NM_DENSE:
  {
    printf("========== storageType = NM_DENSE\n");
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    assert(m->matrix1);
    printf("========== storageType =  NM_SPARSE_BLOCK\n");
    break;
  }
  case NM_SPARSE:
  {
    assert(m->matrix2);
    printf("========== storageType = NM_SPARSE\n");
    switch(m->matrix2->origin)
    {
    case NSM_TRIPLET:
    {
      printf("========== origin =  NSM_TRIPLET\n");
      break;
    }
    case NSM_CSC:
    {
      printf("========== origin =  NSM_CSC\n");
      break;
    }
    case NSM_CSR:
    {
      printf("========== origin =  NSM_CSR\n");
      break;
    }
    default:
    {
      fprintf(stderr, "NM_display ::  unknown origin %d for sparse matrix\n", m->matrix2->origin);
    }
    }

    printf("========== size0 = %i, size1 = %i\n", m->size0, m->size1);
    if(m->matrix2->triplet)
    {
      printf("========== a matrix in format triplet is stored\n");
    }
    if(m->matrix2->csc)
    {
      printf("========== a matrix in format csc is stored\n");
    }
    if(m->matrix2->trans_csc)
    {
      printf("========== a matrix in format trans_csc is stored\n");
    }

    break;
  }
  default:
  {
    fprintf(stderr, "display_storageType for NumericsMatrix: matrix type %d not supported!\n", m->storageType);
  }
  }
}
void NM_display(const NumericsMatrix* const m)
{
  if(! m)
  {
    fprintf(stderr, "Numerics, NumericsMatrix display failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  printf("========== Numerics Matrix\n");

  printf("========== size0 = %i, size1 = %i\n", m->size0, m->size1);

  switch(m->storageType)
  {
  case NM_DENSE:
  {
    printf("========== storageType = NM_DENSE\n");
    NM_dense_display(m->matrix0, m->size0, m->size1, m->size0);
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    assert(m->matrix1);
    printf("========== storageType =  NM_SPARSE_BLOCK\n");
    SBM_print(m->matrix1);
    break;
  }
  case NM_SPARSE:
  {
    assert(m->matrix2);
    printf("========== storageType = NM_SPARSE\n");
    switch(m->matrix2->origin)
    {
    case NSM_TRIPLET:
    {
      printf("========== origin =  NSM_TRIPLET\n");
      break;
    }
    case NSM_HALF_TRIPLET:
    {
      printf("========== origin =  NSM_HALF_TRIPLET\n");
      break;
    }
    case NSM_CSC:
    {
      printf("========== origin =  NSM_CSC\n");
      break;
    }
    case NSM_CSR:
    {
      printf("========== origin =  NSM_CSR\n");
      break;
    }
    default:
    {
      fprintf(stderr, "NM_display ::  unknown origin %d for sparse matrix\n", m->matrix2->origin);
    }
    }

    printf("========== size0 = %i, size1 = %i\n", m->size0, m->size1);
    if(m->matrix2->triplet)
    {
      printf("========== a matrix in format triplet is stored\n");
      cs_print(m->matrix2->triplet, 0);
    }
    if(m->matrix2->csc)
    {
      printf("========== a matrix in format csc is stored\n");
      cs_print(m->matrix2->csc, 0);
    }
    if(m->matrix2->trans_csc)
    {
      printf("========== a matrix in format trans_csc is stored\n");
      cs_print(m->matrix2->trans_csc, 0);
    }
    /* else */
    /* { */
    /*   fprintf(stderr, "display for sparse matrix: no matrix found!\n"); */
    /* } */
    if(m->matrix2->diag_indx)
    {
      printf("========== m->matrix2->diag_indx = %p\n", m->matrix2->diag_indx);
      for(int  i = 0; i < m->size0; ++i) printf("diag_indices[%i] = %li\t ", i, m->matrix2->diag_indx[i]);
    }
    else
    {
      printf("========== m->matrix2->diag_indx --> NULL\n");
    }
    if(m->matrix2->linearSolverParams)
    {
      printf("========== m->matrix2->linearSolverParams = %p \n", m->matrix2->linearSolverParams);
    }
    else
    {
      printf("========== m->matrix2->linearSolverParams --> NULL\n");
    }


    break;
  }
  default:
  {
    fprintf(stderr, "display for NumericsMatrix: matrix type %d not supported!\n", m->storageType);
  }
  }

  if(m->internalData)
  {
    printf("========== internalData->iWorkSize = %lu\n", m->internalData->iWorkSize);
    printf("========== internalData->iWork = %p\n", m->internalData->iWork);

    printf("========== internalData->dWorkSize = %lu\n", (unsigned long) m->internalData->dWorkSize);
    printf("========== internalData->dWork = %p\n", m->internalData->dWork);
  }
  else
  {
    printf("========== internalData = NULL\n");
  }
  if (NM_destructible(m))
    printf("========== is destructible \n");
  else
    printf("========== is not destructible \n");
  if (NM_LU_factorized(m))
    printf("========== is LU factorized \n");
  else
    printf("========== is not LU factorized \n");
  if (NM_Cholesky_factorized(m))
    printf("========== is Cholesky factorized \n");
  else
    printf("========== is not Cholesky factorized \n");
  if (NM_LDLT_factorized(m))
    printf("========== is LDLT factorized \n");
  else
    printf("========== is not LDLT factorized \n");
}

void NM_display_row_by_row(const NumericsMatrix* const m)
{
  if(! m)
  {
    fprintf(stderr, "Numerics, NumericsMatrix display failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  NM_types storageType = m->storageType;
  if(storageType == NM_DENSE)
  {
    printf("\n ========== Numerics Matrix of dim %dX%d\n", m->size0, m->size1);
    for(int lin = 0; lin < m->size0; lin++)
    {
      for(int col = 0; col < m->size1; col++)
        printf("%lf ", m->matrix0[lin + col * m->size1]);
      printf("\n");
    }
  }
  else if(storageType == NM_SPARSE_BLOCK)
    SBM_print(m->matrix1);
  else
  {
    fprintf(stderr, "NM_display_row_by_row :: unknown matrix storage");
    exit(EXIT_FAILURE);
  }
}

void NM_write_in_file(const NumericsMatrix* const m, FILE* file)
{
  DEBUG_PRINT("\n  ========== NM_write_in_file(const NumericsMatrix* const m, FILE* file) start\n");

  if(! m)
  {
    fprintf(stderr, "Numerics, NM_write_in_file failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }

  fprintf(file, "%d\n", m->storageType);
  fprintf(file, "%d\n", m->size0);
  fprintf(file, "%d\n", m->size1);
  DEBUG_PRINTF("\n ========== storageType = %i\n", m->storageType);

  switch(m->storageType)
  {
  case NM_DENSE:
  {
    fprintf(file, "%i\t%i\n", m->size0, m->size1);
    for(int i = 0; i < m->size1 * m->size0; i++)
    {
      fprintf(file, "%32.24e ", m->matrix0[i]);
      if((i + 1) % m->size1 == 0)
        fprintf(file, "\n");
    }
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    assert(m->matrix1);
    SBM_write_in_file(m->matrix1, file);
    break;
  }
  case NM_SPARSE:
  {
    NSM_write_in_file(m->matrix2, file);
    break;
  }
  default:
  {
    fprintf(stderr, "Numerics, NM_write_in_file failed, unknown storage type .\n");
    exit(EXIT_FAILURE);
  }
  }
  DEBUG_PRINT("\n  ========== NM_write_in_file(const NumericsMatrix* const m, FILE* file) end\n");

}
void NM_write_in_file_python(const NumericsMatrix* const m, FILE* file)
{
  if(! m)
  {
    fprintf(stderr, "Numerics, NumericsMatrix_write_in_file_python  failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  fprintf(file, "storageType = %d ; \n", m->storageType);
  fprintf(file, "size0 = %d; \n", m->size0);
  fprintf(file, "size1 = %d; \n", m->size1);
  fprintf(file, "data= [");
  for(int i = 0; i < m->size0; i++)
  {
    fprintf(file, "[");
    for(int j = 0; j < m->size1; j++)
    {
      fprintf(file, "%32.24e,\t ", NM_get_value((NumericsMatrix*) m,i,j));
    }
    fprintf(file, "],\n");
  }
  fprintf(file, "]");
}

void NM_write_in_file_scilab(const NumericsMatrix* const m, FILE* file)
{
  if(! m)
  {
    fprintf(stderr, "Numerics, NumericsMatrix printInFile failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  fprintf(file, "storageType = %d ; \n", m->storageType);
  fprintf(file, "size0 = %d; \n", m->size0);
  fprintf(file, "size1 = %d; \n", m->size1);
  fprintf(file, "data= [");
  for(int i = 0; i < m->size0; i++)
  {
    fprintf(file, "[");
    for(int j = 0; j < m->size1; j++)
    {
      fprintf(file, "%32.24e,\t ", NM_get_value((NumericsMatrix*) m,i,j));
    }
    fprintf(file, "];\n");
  }
  fprintf(file, "]");
}

void NM_write_in_filename(const NumericsMatrix* const m, const char *filename)
{
  FILE* foutput = fopen(filename, "w");
  NM_write_in_file(m, foutput);
  fclose(foutput);
}

void NM_read_in_filename(NumericsMatrix* const m, const char *filename)
{
  FILE* finput = fopen(filename, "r");
  if(finput == NULL)
  {
    puts("Error while opening file");
  }
  else
  {
    NM_read_in_file(m, finput);
    fclose(finput);
  }
}

void NM_read_in_file(NumericsMatrix* const m, FILE *file)
{
  if(!m)
  {
    fprintf(stderr, "Numerics, NumericsMatrix NM_read_in_file failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  CHECK_IO(fscanf(file, "%d", &(m->storageType)));
  CHECK_IO(fscanf(file, "%d", &(m->size0)));
  CHECK_IO(fscanf(file, "%d", &(m->size1)));
  NM_types storageType = m->storageType;

  if(storageType == NM_DENSE)
  {
    CHECK_IO(fscanf(file, "%d\t%d\n", &(m->size0), &(m->size1)));

    for(int i = 0; i < m->size1 * m->size0; i++)
    {
      CHECK_IO(fscanf(file, "%le ", &(m->matrix0[i])));
    }
  }
  else if(storageType == NM_SPARSE_BLOCK)
  {
    SBM_read_in_file(m->matrix1, file);
  }
  else if(storageType == NM_SPARSE)
  {
    NumericsMatrix * tmp = NM_new_from_file(file);
    NM_copy(tmp,m);
    NM_clear(tmp);
  }
  else
  {
    printf("NM_read_in_file :: unknown matrix storage");
    exit(EXIT_FAILURE);
  }
}


NumericsMatrix* NM_new_from_file(FILE *file)
{
  NumericsMatrix* m = NM_new();

  NM_types storageType;
  size_t size0;
  size_t size1;
  int info = 0;
  void* data = NULL;

  CHECK_IO(fscanf(file, "%d", &storageType), &info);
  CHECK_IO(fscanf(file, SN_SIZE_T_F, &size0), &info);
  CHECK_IO(fscanf(file, SN_SIZE_T_F, &size1), &info);

  if(storageType == NM_DENSE)
  {
    CHECK_IO(fscanf(file, SN_SIZE_T_F "\t" SN_SIZE_T_F "\n", &size0, &size1), &info);

    data =  malloc(size1 * size0 * sizeof(double));
    double* data_d = (double *) data;

    for(size_t i = 0; i < size1 * size0; ++i)
    {
      CHECK_IO(fscanf(file, "%le ", &(data_d[i])), &info);
    }
  }
  else if(storageType == NM_SPARSE_BLOCK)
  {
    data = SBM_new_from_file(file);
  }
  else
  {
    data = NSM_new_from_file(file);
  }

  NM_fill(m, storageType, (int)size0, (int)size1, data);

  return m;
}

NumericsMatrix* NM_new_from_filename(const char * filename)
{
  FILE* finput = fopen(filename, "r");
  if(finput == NULL)
  {
    puts("Error while opening file");
  }
  else
  {
    NumericsMatrix * A= NM_new_from_file(finput);
    fclose(finput);
    return A;

  }
  return NULL;
}

NumericsMatrix* NM_create_from_file(FILE *file)
{
  return NM_new_from_file(file);
}

NumericsMatrix* NM_create_from_filename(const char * filename)
{
  FILE* finput = fopen(filename, "r");
  if(finput == NULL)
  {
    puts("Error while opening file");
  }
  else
  {
    NumericsMatrix * A= NM_create_from_file(finput);
    fclose(finput);
    return A;

  }
  return NULL;
}



void NM_read_in_file_scilab(NumericsMatrix* const M, FILE *file)
{
  fprintf(stderr, "Numerics, NumericsMatrix,NM_read_in_file_scilab");
  exit(EXIT_FAILURE);
}

void NM_extract_diag_block(NumericsMatrix* M, int block_row_nb, size_t start_row, int size, double ** Block)
{
  NM_types storageType = M->storageType;
  switch(storageType)
  {
  case NM_DENSE:
  {
    double* Mptr = M->matrix0 + (M->size0 + 1)*start_row;
    double* Bmat = *Block;
    /* The part of MM which corresponds to the current block is copied into MLocal */
    for(size_t i = 0; i < (size_t) size; ++i)
    {
      memcpy(Bmat, Mptr, (size_t)size*sizeof(double));
      Mptr += M->size0;
      Bmat += size;
    }
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    int diagPos = SBM_diagonal_block_index(M->matrix1, block_row_nb);
    (*Block) = M->matrix1->block[diagPos];
    break;
  }
  case NM_SPARSE:
  {
    NSM_extract_block(M, *Block, start_row, start_row, size, size);
    break;
  }
  default:
  {
    printf("NM_extract_diag_block :: unknown matrix storage");
    exit(EXIT_FAILURE);
  }
  }
}
void NM_extract_diag_block3(NumericsMatrix* M, int block_row_nb, double ** Block)
{
  NM_types storageType = M->storageType;
  switch(storageType)
  {
  case NM_DENSE:
  {
    double* Mptr = M->matrix0 + (M->size0 + 1)*(block_row_nb + block_row_nb + block_row_nb);
    double* Bmat = *Block;
    /* The part of MM which corresponds to the current block is copied into MLocal */
    Bmat[0] = Mptr[0];
    Bmat[1] = Mptr[1];
    Bmat[2] = Mptr[2];
    Mptr += M->size0;
    Bmat[3] = Mptr[0];
    Bmat[4] = Mptr[1];
    Bmat[5] = Mptr[2];
    Mptr += M->size0;
    Bmat[6] = Mptr[0];
    Bmat[7] = Mptr[1];
    Bmat[8] = Mptr[2];
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    int diagPos = SBM_diagonal_block_index(M->matrix1, block_row_nb);
    (*Block) = M->matrix1->block[diagPos];
    break;
  }
  case NM_SPARSE:
  {
    size_t start_row = (size_t)block_row_nb + block_row_nb + block_row_nb;
    NSM_extract_block(M, *Block, start_row, start_row, 3, 3);
    break;
  }
  default:
  {
    printf("NM_extract_diag_block :: unknown matrix storage");
    exit(EXIT_FAILURE);
  }
  }
}
void NM_extract_diag_block2(NumericsMatrix* M, int block_row_nb, double ** Block)
{
  NM_types storageType = M->storageType;
  switch(storageType)
  {
  case NM_DENSE:
  {
    double* Mptr = M->matrix0 + (M->size0 + 1)*(block_row_nb + block_row_nb);
    double* Bmat = *Block;
    /* The part of MM which corresponds to the current block is copied into MLocal */
    Bmat[0] = Mptr[0];
    Bmat[1] = Mptr[1];
    Mptr += M->size0;
    Bmat[2] = Mptr[0];
    Bmat[3] = Mptr[1];
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    int diagPos = SBM_diagonal_block_index(M->matrix1, block_row_nb);
    (*Block) = M->matrix1->block[diagPos];
    break;
  }
  case NM_SPARSE:
  {
    size_t start_row = (size_t)block_row_nb + block_row_nb;
    NSM_extract_block(M, *Block, start_row, start_row, 2, 2);
    break;
  }
  default:
  {
    printf("NM_extract_diag_block :: unknown matrix storage");
    exit(EXIT_FAILURE);
  }
  }
}





void NM_extract_diag_block5(NumericsMatrix* M, int block_row_nb, double ** Block)
{
  NM_types storageType = M->storageType;
  switch(storageType)
  {
  case NM_DENSE:
  {
    double* Mptr = M->matrix0 + (M->size0 + 1)*(block_row_nb
                   + block_row_nb
                   + block_row_nb
                   + block_row_nb
                   + block_row_nb);
    double* Bmat = *Block;
    /* The part of MM which corresponds to the current block is copied into MLocal */
    Bmat[0] = Mptr[0];
    Bmat[1] = Mptr[1];
    Bmat[2] = Mptr[2];
    Bmat[3] = Mptr[3];
    Bmat[4] = Mptr[4];
    Mptr += M->size0;
    Bmat[5] = Mptr[0];
    Bmat[6] = Mptr[1];
    Bmat[7] = Mptr[2];
    Bmat[8] = Mptr[3];
    Bmat[9] = Mptr[4];
    Mptr += M->size0;
    Bmat[10] = Mptr[0];
    Bmat[11] = Mptr[1];
    Bmat[12] = Mptr[2];
    Bmat[13] = Mptr[3];
    Bmat[14] = Mptr[4];
    Mptr += M->size0;
    Bmat[15] = Mptr[0];
    Bmat[16] = Mptr[1];
    Bmat[17] = Mptr[2];
    Bmat[18] = Mptr[3];
    Bmat[19] = Mptr[4];
    Mptr += M->size0;
    Bmat[20] = Mptr[0];
    Bmat[21] = Mptr[1];
    Bmat[22] = Mptr[2];
    Bmat[23] = Mptr[3];
    Bmat[24] = Mptr[4];
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    assert(0); /* this has to be checked carefully */
    int diagPos = SBM_diagonal_block_index(M->matrix1, block_row_nb);
    (*Block) = M->matrix1->block[diagPos];
    break;
  }
  case NM_SPARSE:
  {
    size_t start_row = (size_t)5*block_row_nb;
    NSM_extract_block(M, *Block, start_row, start_row, 5, 5);
    break;
  }
  default:
  {
    printf("NM_extract_diag_block :: unknown matrix storage");
    exit(EXIT_FAILURE);
  }
  }
}

void NM_copy_diag_block3(NumericsMatrix* M, int block_row_nb, double ** Block)
{
  NM_types storageType = M->storageType;
  switch(storageType)
  {
  case NM_DENSE:
  {
    double* Mptr = M->matrix0 + (M->size0 + 1)*(block_row_nb + block_row_nb + block_row_nb);
    double* Bmat = *Block;
    /* The part of MM which corresponds to the current block is copied into MLocal */
    Bmat[0] = Mptr[0];
    Bmat[1] = Mptr[1];
    Bmat[2] = Mptr[2];
    Mptr += M->size0;
    Bmat[3] = Mptr[0];
    Bmat[4] = Mptr[1];
    Bmat[5] = Mptr[2];
    Mptr += M->size0;
    Bmat[6] = Mptr[0];
    Bmat[7] = Mptr[1];
    Bmat[8] = Mptr[2];
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    int diagPos = SBM_diagonal_block_index(M->matrix1, block_row_nb);
    double* Mptr = M->matrix1->block[diagPos];
    double* Bmat = *Block;
    /* The part of MM which corresponds to the current block is copied into MLocal */
    Bmat[0] = Mptr[0];
    Bmat[1] = Mptr[1];
    Bmat[2] = Mptr[2];
    Bmat[3] = Mptr[3];
    Bmat[4] = Mptr[4];
    Bmat[5] = Mptr[5];
    Bmat[6] = Mptr[6];
    Bmat[7] = Mptr[7];
    Bmat[8] = Mptr[8];
    break;
  }
  case NM_SPARSE:
  {
    size_t start_row = (size_t)block_row_nb + block_row_nb + block_row_nb;
    NSM_extract_block(M, *Block, start_row, start_row, 3, 3);
    break;
  }
  default:
  {
    printf("NM_extract_diag_block :: unknown matrix storage");
    exit(EXIT_FAILURE);
  }
  }
}

void NM_add_to_diag3(NumericsMatrix* M, double alpha)
{
  size_t n = M->size0;
  switch(M->storageType)
  {
  case NM_DENSE:
  {
    for(size_t indx = 0; indx < n*n; indx += n+1) M->matrix0[indx] += alpha;
    NM_inc_version(M, NM_DENSE);
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    for(size_t ic = 0; ic < n/3; ++ic)
    {
      int diagPos = SBM_diagonal_block_index(M->matrix1, ic);
      M->matrix1->block[diagPos][0] += alpha;
      M->matrix1->block[diagPos][4] += alpha;
      M->matrix1->block[diagPos][8] += alpha;
    }
    NM_inc_version(M, NM_SPARSE_BLOCK);
    break;
  }
  case NM_SPARSE:
  {
    /* NSM_diag_indices modifies M->matrix2->origin */
    CS_INT* diag_indices = NSM_diag_indices(M);

    DEBUG_EXPR(
      printf("diag_indices:\n");
      for(size_t i = 0; i < n; ++i) printf("diag_indices[%zu] = %li\t ", i, diag_indices[i]);
    );


    /* assert (NM_version(M, M->matrix2->origin) == */
    /*         NSM_version(M->matrix2, NSM_latest_id(M->matrix2))); */
    assert (NSM_version(M->matrix2, M->matrix2->origin) ==
            NSM_version(M->matrix2, NSM_latest_id(M->matrix2)));

    double* Mx = NSM_data(M->matrix2);
    for(size_t i = 0; i < n; ++i) Mx[diag_indices[i]] += alpha;

    NSM_inc_version(M->matrix2, M->matrix2->origin);
    break;
  }
  default:
    printf("NM_add_to_diag3 :: unsupported matrix storage %d", M->storageType);
    exit(EXIT_FAILURE);
  }
}
void NM_add_to_diag5(NumericsMatrix* M, double alpha)
{
  size_t n = M->size0;
  switch(M->storageType)
  {
  case NM_DENSE:
  {
    for(size_t indx = 0; indx < n*n; indx += n+1) M->matrix0[indx] += alpha;
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    for(size_t ic = 0; ic < n/5; ++ic)
    {
      int diagPos = SBM_diagonal_block_index(M->matrix1, ic);
      M->matrix1->block[diagPos][0] += alpha;
      M->matrix1->block[diagPos][6] += alpha;
      M->matrix1->block[diagPos][12] += alpha;
      M->matrix1->block[diagPos][18] += alpha;
      M->matrix1->block[diagPos][24] += alpha;
    }
    break;
  }
  case NM_SPARSE:
  {
    CS_INT* diag_indices = NSM_diag_indices(M);

    DEBUG_EXPR(
      printf("diag_indices:\n");
      for(size_t i = 0; i < n; ++i) printf("diag_indices[%zu] = %li\t ", i, diag_indices[i]);
    );

    double* Mx = NSM_data(M->matrix2);
    for(size_t i = 0; i < n; ++i) Mx[diag_indices[i]] += alpha;

    break;
  }
  default:
    printf("NM_add_to_diag5 :: unsupported matrix storage %d", M->storageType);
    exit(EXIT_FAILURE);
  }
}

NumericsMatrix *  NM_add(double alpha, NumericsMatrix* A, double beta, NumericsMatrix* B)
{


  assert(A->size0 == B->size0 && "NM_add :: A->size0 != B->size0 ");
  assert(A->size1 == B->size1 && "NM_add :: A->size1 != B->size1 ");


  /* The storageType  for C inherits from A except for NM_SPARSE_BLOCK */
  NumericsMatrix *C = NM_create(A->storageType, A->size0, A->size1);

  /* should we copy the whole internal data ? */
  /*NM_internalData_copy(A, C);*/
  NM_MPI_copy(A, C);
  NM_MUMPS_copy(A, C);

  switch(A->storageType)
  {
  case NM_DENSE:
  {
    int nm= A->size0*A->size1;
    cblas_dcopy(nm, A->matrix0, 1, C->matrix0, 1);
    cblas_dscal(nm, alpha, C->matrix0,1);
    switch(B->storageType)
    {
    case NM_DENSE:
    {
      cblas_daxpy(nm, beta, B->matrix0, 1, C->matrix0,1);
      NM_inc_version(C, NM_DENSE);
      break;
    }
    case NM_SPARSE_BLOCK:
    case NM_SPARSE:
    {
      NumericsMatrix* B_dense = NM_create(NM_DENSE, A->size0, A->size1);
      NM_to_dense(B, B_dense);
      /* MB: where is cleaned B_dense ? */
      cblas_daxpy(nm, beta, B_dense->matrix0, 1, C->matrix0,1);
      NM_inc_version(C, NM_DENSE);
      break;
    }
    default:
    {
      numerics_error("NM_add","Unsupported storage type %d, exiting!\n", B->storageType);
      exit(EXIT_FAILURE);
      break;
    }
    }
    break;
  }
  case NM_SPARSE_BLOCK:
  case NM_SPARSE:
  {

    CSparseMatrix* result = cs_add(NM_csc(A), NM_csc(B), alpha, beta);
    assert(result && "NM_add :: cs_add failed");
    NSM_fix_csc(result);
    NumericsSparseMatrix* C_nsm  = numericsSparseMatrix(C);

    C_nsm->csc = result;
    C_nsm->origin = NSM_CSC;
    C->storageType=NM_SPARSE;

    NSM_set_version(C->matrix2, NSM_CSC, NM_max_version(C));
    NSM_inc_version(C->matrix2, NSM_CSC);
    break;
  }
  default:
  {
    numerics_error("NM_add:","unsupported matrix storage %d", A->storageType);
  }
  }
  return C;

}

void  NM_scal(double alpha, NumericsMatrix* A)
{

  switch(A->storageType)
  {
  case NM_DENSE:
  {
    int nm= A->size0*A->size1;
    cblas_dscal(nm, alpha, A->matrix0,1);
    NM_inc_version(A, NM_DENSE);
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    /* version incremented in SBM_scal */
    SBM_scal(alpha, A->matrix1);
    break;
  }
  case NM_SPARSE:
  {
    CSparseMatrix_scal(alpha, NM_csc(A));
    NSM_inc_version(A->matrix2, NSM_CSC);
    A->matrix2->origin = NSM_CSC;
    /* Invalidations */
    NM_clearTriplet(A);
    NM_clearHalfTriplet(A);
    NM_clearCSCTranspose(A);
    NM_clearCSR(A);
    break;
  }
  default:
  {
    numerics_error("NM_scal:","unsupported matrix storage %d", A->storageType);
  }
  }
  return;

}


NumericsMatrix* NM_create_from_data(int storageType, int size0, int size1, void* data)
{
  NumericsMatrix* M = NM_new();

  NM_fill(M, storageType, size0, size1, data);

  return M;
}

NumericsMatrix* NM_duplicate(NumericsMatrix* mat)
{
  NumericsMatrix* M = NM_new();

  void* data;
  int size0 = mat->size0;
  int size1 = mat->size1;

  M->storageType = mat->storageType;
  switch(mat->storageType)
  {
  case NM_DENSE:
    data = malloc((size_t)size0*size1*sizeof(double));
    break;
  case NM_SPARSE_BLOCK:
    data = SBM_new();
    break;
  case NM_SPARSE:
    data = NSM_new();
    break;
  default:
    fprintf(stderr, "NM_duplicate :: storageType value %d not implemented yet !", mat->storageType);
    exit(EXIT_FAILURE);
  }

  NM_fill(M, mat->storageType, size0, size1, data);

  return M;
}

NumericsMatrix* NM_new(void)
{
  DEBUG_BEGIN("NumericsMatrix* NM_new(void)\n");
  NumericsMatrix* M = (NumericsMatrix*) malloc(sizeof(NumericsMatrix));
  M->storageType = NM_UNKNOWN;
  M->size0 = 0;
  M->size1 = 0;
  NM_null(M);
  DEBUG_END("NumericsMatrix* NM_new(void)\n");
  return M;
}

NumericsMatrix* NM_eye(int size)
{
  NumericsMatrix* M = NM_create(NM_SPARSE, size, size);
  /* version incremented in NSM_triplet_eye */
  NSM_clear(M->matrix2);
  free(M->matrix2);
  M->matrix2 = NSM_triplet_eye(size);
  return M;
}

NumericsMatrix* NM_scalar(int size, double s)
{
  NumericsMatrix* M = NM_create(NM_SPARSE, size, size);
  M->matrix2 = NSM_triplet_scalar(size, s);
  return M;
}

NumericsMatrix* NM_create(NM_types storageType, int size0, int size1)
{
  NumericsMatrix* M = NM_new();

  void* data;

  switch(storageType)
  {
  case NM_DENSE:
    data = calloc(size0*size1,sizeof(double));
    break;
  case NM_SPARSE_BLOCK:
    data = SBM_new();
    break;
  case NM_SPARSE:
    data = NSM_new();
    break;
  default:
    data=NULL;
    numerics_error("NM_create", "storageType value %d not implemented yet !", storageType);
  }

  NM_fill(M, storageType, size0, size1, data);

  return M;
}


void NM_fill(NumericsMatrix* M, NM_types storageType, int size0, int size1, void* data)
{

  assert(M);
  NM_null(M);
  M->storageType = storageType;
  M->size0 = size0;
  M->size1 = size1;



  if(data)
  {
    switch(storageType)
    {
    case NM_DENSE:
      M->matrix0 = (double*) data;
      NM_inc_version(M, NM_DENSE);
      break;
    case NM_SPARSE_BLOCK:
      M->matrix1 = (SparseBlockStructuredMatrix*) data;
      NM_inc_version(M, NM_SPARSE_BLOCK);
      break;
    case NM_SPARSE:
      M->matrix2 = (NumericsSparseMatrix*) data;
      if(data)
      {
        if(M->matrix2->origin == NSM_UNKNOWN)
        {
          if(M->matrix2->triplet)
          {
            M->matrix2->origin = NSM_TRIPLET;
            NSM_inc_version(M->matrix2, NSM_TRIPLET);
          }
          else if(M->matrix2->half_triplet)
          {
            M->matrix2->origin = NSM_HALF_TRIPLET;
            NSM_inc_version(M->matrix2, NSM_HALF_TRIPLET);
          }
          else if(M->matrix2->csc)
          {
            M->matrix2->origin = NSM_CSC;
            NSM_inc_version(M->matrix2, NSM_CSC);
          }
          else if(M->matrix2->csr)
          {
            M->matrix2->origin = NSM_CSR;
            NSM_inc_version(M->matrix2, NSM_CSR);
          }
        }
      }
      break;

    default:
      printf("NM_fill :: storageType value %d not implemented yet !", storageType);
      exit(EXIT_FAILURE);
    }
  }
}

NumericsMatrix* NM_new_SBM(int size0, int size1, SparseBlockStructuredMatrix* m1)
{
  return NM_create_from_data(NM_SPARSE_BLOCK, size0, size1, (void*)m1);
}
NumericsMatrix* NM_transpose(NumericsMatrix * A)
{
  NumericsMatrix* Atrans;
  switch(A->storageType)
  {
  case NM_DENSE:
  {
    Atrans = NM_create(NM_DENSE, A->size1, A->size0);
    for(int i = 0; i < Atrans->size0; i++)
    {
      for(int j = 0; j < Atrans->size1; j++)
      {
        Atrans->matrix0[i+j*Atrans->size0] = A->matrix0[j+i*A->size0];
      }
    }
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    Atrans = NM_create(NM_SPARSE_BLOCK, A->size1, A->size0);
    SBM_transpose(A->matrix1, Atrans->matrix1);
    break;
  }
  case NM_SPARSE:
  {
    assert(A->matrix2);
    Atrans = NM_create(NM_SPARSE,A->size1,A->size0);
    NM_csc_alloc(Atrans, 0);
    Atrans->matrix2->origin = NSM_CSC;
    // \todo should be a copy */
    CSparseMatrix_copy(NM_csc_trans(A), NM_csc(Atrans));
    //DEBUG_EXPR(NM_display(Atrans););
    break;
  }
  default:
  {
    numerics_error("NM_to_dense","Unsupported storage type %d, exiting!\n", A->storageType);
    exit(EXIT_FAILURE);
  }
  }
  NM_MPI_copy(A, Atrans);
  NM_MUMPS_copy(A, Atrans);

  return Atrans;
}

bool NM_destructible(const NumericsMatrix* A)
{
  return A->destructible == A;
}

RawNumericsMatrix* NM_preserve(NumericsMatrix* A)
{
  if (NM_destructible(A))
  {
    if (NM_LU_factorized(A) || NM_Cholesky_factorized(A) ||  NM_LDLT_factorized(A) )
    {
      numerics_warning("NM_preserve", "preservation is done on a factorized matrix!");
    }
    NumericsMatrix* B = NM_new();
    NM_copy(A, B);

    A->destructible = B;
    assert(A->destructible->destructible == A->destructible);
  };
  return A;
}

RawNumericsMatrix* NM_unpreserve(NumericsMatrix* A)
{
  if (A->destructible != A)
  {
    NM_clear(A->destructible);
    free(A->destructible);
    A->destructible = A;
  }
  return A;
}

bool NM_LU_factorized(const NumericsMatrix* const A)
{
  return NM_internalData(A->destructible)->isLUfactorized;
}
bool NM_Cholesky_factorized(const NumericsMatrix* const A)
{
  return NM_internalData(A->destructible)->isCholeskyfactorized;
}
bool NM_LDLT_factorized(const NumericsMatrix* const A)
{
  return NM_internalData(A->destructible)->isLDLTfactorized;
}

void NM_set_LU_factorized(NumericsMatrix* A, bool flag)
{
  NM_internalData(A->destructible)->isLUfactorized = flag;
}

void NM_set_Cholesky_factorized(NumericsMatrix* A, bool flag)
{
  NM_internalData(A->destructible)->isCholeskyfactorized = flag;
}

void NM_set_LDLT_factorized(NumericsMatrix* A, bool flag)
{
  NM_internalData(A->destructible)->isLDLTfactorized = flag;
}



void NM_clearDense(NumericsMatrix* A)
{
  if(A->matrix0)
  {
    free(A->matrix0);
  }
  A->matrix0 = NULL;
  NM_reset_version(A, NM_DENSE);
}

void NM_clearSparseBlock(NumericsMatrix* A)
{
  if(A->matrix1)
  {
    SBM_clear(A->matrix1);
    free(A->matrix1);
  }
  A->matrix1 = NULL;
  /* no need to reset version! */
}

void NM_clearSparse(NumericsMatrix* A)
{
  if(A->matrix2)
  {
    NSM_clear(A->matrix2);
    free(A->matrix2);
  }
  A->matrix2 = NULL;
  /* no need to reset version! */
}

void NM_clearTriplet(NumericsMatrix* A)
{
  if(A->matrix2)
  {
    if(A->matrix2->triplet)
    {
      cs_spfree(A->matrix2->triplet);
    }
    A->matrix2->triplet = NULL;
    NSM_reset_version(A->matrix2, NSM_TRIPLET);
  }
}

void NM_clearHalfTriplet(NumericsMatrix* A)
{
  if(A->matrix2)
  {
    if(A->matrix2->half_triplet)
    {
      cs_spfree(A->matrix2->half_triplet);
      A->matrix2->half_triplet = NULL;
    }
    NSM_reset_version(A->matrix2, NSM_HALF_TRIPLET);
  }
}

void NM_clearCSC(NumericsMatrix* A)
{
  if(A->matrix2)
  {
    if(A->matrix2->csc)
    {
      cs_spfree(A->matrix2->csc);
    }
    A->matrix2->csc = NULL;
    NSM_reset_version(A->matrix2, NSM_CSC);
  }
}

void NM_clearCSCTranspose(NumericsMatrix* A)
{
  if(A->matrix2)
  {
    if(A->matrix2->trans_csc)
    {
      cs_spfree(A->matrix2->trans_csc);
    }
    A->matrix2->trans_csc = NULL;
  }
  /* no version for csc transpose as it is a terminal format */
}

void NM_clearCSR(NumericsMatrix* A)
{
  if(A->matrix2)
  {
    if(A->matrix2->csr)
    {
      cs_spfree(A->matrix2->csr);
    }
    A->matrix2->csr = NULL;
    NSM_reset_version(A->matrix2, NSM_CSR);
  }
}

void NM_clearSparseStorage(NumericsMatrix *A)
{
  if(A->matrix2)
  {
    A->matrix2->origin = NSM_UNKNOWN;
    if(A->matrix2->linearSolverParams)
      A->matrix2->linearSolverParams = NSM_linearSolverParams_free(A->matrix2->linearSolverParams);
  }
  NM_clearTriplet(A);
  NM_clearHalfTriplet(A);
  NM_clearCSC(A);
  NM_clearCSCTranspose(A);
  NM_clearCSR(A);

  /* reset version done in NM_clear* */
}


void NM_dense_to_sparse(const NumericsMatrix* const A, NumericsMatrix* B, double threshold)
{
  assert(A->matrix0);
  assert(B->matrix2->triplet);
  for(int i = 0; i < A->size0; ++i)
  {
    for(int j = 0; j < A->size1; ++j)
    {
      CHECK_RETURN(CSparseMatrix_zentry(B->matrix2->triplet, i, j, A->matrix0[i + A->size0*j], threshold));
    }
  }
  if (A == B)
  {
    /* on the same matrix, the versions are the same */
    NSM_set_version(B->matrix2, NSM_TRIPLET, NM_version(A, NM_DENSE));
  }
  else
  {
    /* increment the version to the max */
    NSM_inc_version(B->matrix2, NSM_TRIPLET);
  }
}
int NM_to_dense(const NumericsMatrix* const A, NumericsMatrix* B)
{
  int info = 1;
  if(!B->matrix0)
  {
    B->matrix0 = (double *)calloc(A->size0*A->size1, sizeof(double));
  }
  else if(B->size0 != A->size0 || B->size0 != A->size0)
  {
    free(B->matrix0);
    B->matrix0 = (double *)calloc(A->size0*A->size1, sizeof(double));
  }

  assert(B->matrix0);

  B->size0 = A->size0;
  B->size1 = A->size1;
  B->storageType=NM_DENSE;

  unsigned long src_version;
  switch(A->storageType)
  {
  case NM_DENSE:
  {
    NM_copy(A, B);
    info=0;
    src_version = NM_version(A, NM_DENSE);
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    SBM_to_dense(A->matrix1, B->matrix0);
    info=0;
    src_version = NM_version(A, NM_SPARSE_BLOCK);
    break;
  }
  case NM_SPARSE:
  {
    assert(A->matrix2);
    info  = NSM_to_dense(A->matrix2, B->matrix0);
    assert (NSM_version(A->matrix2, NSM_latest_id(A->matrix2)) ==
            NSM_version(A->matrix2, A->matrix2->origin));
    src_version = NSM_version(A->matrix2, A->matrix2->origin);
    break;
  }
  default:
  {
    numerics_error("NM_to_dense","Unsupported storage type %d, exiting!\n", A->storageType);
    exit(EXIT_FAILURE);
  }
  }
  /* invalidations */
  NM_clearSparse(B);
  NM_clearSparseBlock(B);

  if (A == B)
  {
    /* on the same matrix, the versions are the same */
    NM_set_version(B, NM_DENSE, src_version);
  }
  else
  {
    /* increment the version to the max */
    NM_inc_version(B , NM_DENSE);
  }

  return info;


}


void NM_copy_to_sparse(const NumericsMatrix* const A, NumericsMatrix* B, double threshold)
{
  DEBUG_BEGIN("NM_copy_to_sparse(...)\n")
  assert(A);
  assert(B);
  B->size0 = A->size0;
  B->size1 = A->size1;

  assert(B->storageType == NM_SPARSE);
  if(!B->matrix2)
  {
    B->matrix2 = NSM_new();
  }

  switch(A->storageType)
  {
  case NM_DENSE:
  {
    B->matrix2->triplet = cs_spalloc(0,0,1,1,1);
    B->matrix2->origin = NSM_TRIPLET;
    /* version set in NM_dense_to_sparse */
    NM_dense_to_sparse(A, B, threshold );
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    // XXX this is suboptimal since the matrix A might have already been converted
    // to csc or triplet --xhub

    B->matrix1 = A->matrix1;
    B->storageType = NM_SPARSE_BLOCK;
    NM_triplet(B);
    B->matrix1 = NULL;
    B->storageType = NM_SPARSE;

    if (A == B)
    {
      NSM_set_version(B->matrix2, NSM_TRIPLET, NM_version(A, NM_SPARSE_BLOCK));
    }
    else
    {
      NSM_inc_version(B->matrix2, NSM_TRIPLET);
    }
    break;
  }
  case NM_SPARSE:
  {
    /* version set in NM_copy */
    NM_copy(A, B);
    break;
  }
  default:
  {
    printf("NM_copy_to_sparse :: Unsupported storage type %d, exiting!\n", A->storageType);
    exit(EXIT_FAILURE);
  }
  }
  DEBUG_END("NM_copy_to_sparse(...)\n")
}

void NM_version_copy(const NumericsMatrix* const A, NumericsMatrix* B)
{
  assert(A);
  assert(B);
  switch(A->storageType)
  {
  case NM_DENSE:
  {
    NM_set_version(B, NM_DENSE, NM_version(A, NM_DENSE));
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    NM_set_version(B, NM_SPARSE_BLOCK, NM_version(A, NM_SPARSE_BLOCK));
    break;
  }
  case NM_SPARSE:
  {
    assert(A->matrix2);
    assert(B->matrix2);
    NSM_version_copy(A->matrix2, B->matrix2);
    break;
  }
  default:
  {
    numerics_error("NM_version_copy", "unknown id");
  }
  assert (false);
  }
}

void NM_copy(const NumericsMatrix* const A, NumericsMatrix* B)
{
  assert(A);
  assert(B);
  int sizeA = A->size0 * A->size1;
  int sizeB = B->size0 * B->size1;
  B->size0 = A->size0;
  B->size1 = A->size1;

  NM_internalData_free(B);

  B->storageType = A->storageType;
  switch(A->storageType)
  {
  case NM_DENSE:
  {
    if(B->matrix0)
    {
      if(sizeB < sizeA)
      {
        B->matrix0 = (double*) realloc(B->matrix0, sizeA * sizeof(double));
      }
    }
    else
    {
      B->matrix0 = (double*) malloc(sizeA * sizeof(double));
    }
    cblas_dcopy(sizeA, A->matrix0, 1, B->matrix0, 1);

    /* invalidations */
    NM_clearSparseBlock(B);
    NM_clearSparseStorage(B);

    NM_set_version(B, NM_DENSE, NM_version(A, NM_DENSE));

    break;
  }
  case NM_SPARSE_BLOCK:
  {
    if(!B->matrix1)
    {
      B->matrix1 = SBM_new();
    }

    SparseBlockStructuredMatrix* A_ = A->matrix1;
    SparseBlockStructuredMatrix* B_ = B->matrix1;

    SBM_copy(A_,B_,1);

    /* invalidations */
    NM_clearDense(B);
    NM_clearSparseStorage(B);

    NM_set_version(B, NM_SPARSE_BLOCK, NM_version(A, NM_SPARSE_BLOCK));
    break;
  }
  case NM_SPARSE:
  {
    NumericsSparseMatrix * A_ = A->matrix2;
    NumericsSparseMatrix * B_ = numericsSparseMatrix(B);

    /* version management done in NSM_copy*/
    NSM_copy(A_,B_);

    /* invalidations */
    NM_clearDense(B);
    NM_clearSparseBlock(B);

    if(NSM_get_origin(B_)->nz >= 0)
    {
      NM_clearCSC(B);
      NM_clearCSCTranspose(B);
      NM_clearCSR(B);
    }
    else
    {
      NM_clearTriplet(B);
      NM_clearHalfTriplet(B);
      if(A->matrix2->origin == NSM_CSC)
      {
        NM_clearCSR(B);
      }
      else
      {
        NM_clearCSC(B);
      }
    }
    break;
  }
  default:
    numerics_error("NM_copy","The type of the source matrix is unknown.");

  }
  NM_internalData_copy(A, B);
  NM_MPI_copy(A, B);
  NM_MUMPS_copy(A, B);

  if (NM_destructible(A))
  {
    /* A is destructible, so B must be destructible */
    NM_unpreserve(B);
    B->destructible = B;
  }
  else
  {
    /* A is preserved, so B must be preserved */
    /* assert(!NM_destructible(B));  VA. 22-09-2020 I do not understand B must be destructible. It the case by default if B is created with NM_new */
    NM_preserve(B);
  }

  assert(NM_destructible(B) == NM_destructible(A));
  assert(NM_max_version(B) == NM_max_version(A));

}

NumericsSparseMatrix* numericsSparseMatrix(NumericsMatrix* A)
{
  if(!A->matrix2)
  {
    A->matrix2 = NSM_new();
  }
  return A->matrix2;
}



CSparseMatrix* NM_triplet(NumericsMatrix* A)
{
  if(!numericsSparseMatrix(A)->triplet)
  {
    switch(A->storageType)
    {
    case NM_DENSE:
    case NM_SPARSE_BLOCK:
    {

      /* Invalidation of previously constructed csc storage. */
      /* If we want to avoid this -> rewrite cs_compress with reallocation. */
      NM_clearCSC(A);
      NM_clearCSR(A);
      NM_clearCSCTranspose(A);

      A->matrix2->origin = NSM_TRIPLET;

      A->matrix2->triplet = cs_spalloc(0,0,1,1,1);

      if(A->matrix1)
      {

        /* iteration on row, cr : current row */
        for(unsigned int cr = 0; cr < A->matrix1->filled1-1; ++cr)
        {
          for(size_t bn = A->matrix1->index1_data[cr];
              bn < A->matrix1->index1_data[cr + 1]; ++bn)
          {
            /* cc : current column */
            size_t cc = A->matrix1->index2_data[bn];
            unsigned int inbr = A->matrix1->blocksize0[cr];
            unsigned int roffset = 0;
            unsigned int coffset = 0;
            if(cr != 0)
            {
              roffset = A->matrix1->blocksize0[cr - 1];
              inbr -= roffset;
            }
            unsigned int inbc = A->matrix1->blocksize1[cc];
            if(cc != 0)
            {
              coffset = A->matrix1->blocksize1[cc - 1];
              inbc -= coffset;
            }
            for(unsigned j = 0; j < inbc; ++j)
            {
              for(unsigned i = 0; i < inbr; ++i)
              {
                CHECK_RETURN(CSparseMatrix_entry(A->matrix2->triplet, i + roffset, j + coffset,
                                                  A->matrix1->block[bn][i + j*inbr]));
              }
            }
          }
        }
        NSM_set_version(A->matrix2, NSM_TRIPLET, NM_version(A,
                                                            NM_SPARSE_BLOCK));
      }
      else if(A->matrix0)
      {
        /* version set in NM_dense_to_sparse */
        NM_dense_to_sparse(A, A, DBL_EPSILON);
      }
      else if(A->size0 > 0 || A->size1 > 0)
      {
        fprintf(stderr, "NM_triplet: sparse matrix cannot be constructed.\n");
        exit(EXIT_FAILURE);
      }
      break;
    }
    case NM_SPARSE:
    {
      switch(A->matrix2->origin)
      {
      case NSM_CSC:
      {
        assert(A->matrix2->csc);
        A->matrix2->triplet = NM_csc_to_triplet(A->matrix2->csc);
        NSM_set_version(A->matrix2, NSM_TRIPLET, NSM_version(A->matrix2,
                                                             NSM_CSC));
        break;
      }
      case NSM_CSR:
      {
        assert(A->matrix2->csr);
        A->matrix2->triplet = NM_csr_to_triplet(A->matrix2->csr);
        NSM_set_version(A->matrix2, NSM_TRIPLET, NSM_version(A->matrix2,
                                                             NSM_CSR));
        break;
      }
      default:
      case NSM_UNKNOWN:
      {
        NSM_UNKNOWN_ERR("NM_triplet", A->matrix2->origin);
        exit(EXIT_FAILURE);
      }
      }
      break;
      default:
      {
        fprintf(stderr, "NM_triplet: unknown matrix type\n");
        exit(EXIT_FAILURE);
      }
    }

    }
  }
  assert(A->matrix2->triplet);

  assert(NM_max_version(A) == NSM_version(A->matrix2, NSM_TRIPLET));

  return A->matrix2->triplet;
}


CSparseMatrix* NM_half_triplet(NumericsMatrix* A)
{
  if(!numericsSparseMatrix(A)->half_triplet)
  {
    switch(A->storageType)
    {
    case NM_DENSE:
    case NM_SPARSE_BLOCK:
    {

      /* Invalidation of previously constructed csc storage. */
      /* If we want to avoid this -> rewrite cs_compress with reallocation. */
      NM_clearCSC(A);
      NM_clearCSR(A);
      NM_clearCSCTranspose(A);

      A->matrix2->origin = NSM_HALF_TRIPLET;

      A->matrix2->half_triplet = cs_spalloc(0,0,1,1,1);

      if(A->matrix1)
      {

        /* iteration on row, cr : current row */
        for(unsigned int cr = 0; cr < A->matrix1->filled1-1; ++cr)
        {
          for(size_t bn = A->matrix1->index1_data[cr];
              bn < A->matrix1->index1_data[cr + 1]; ++bn)
          {
            /* cc : current column */
            size_t cc = A->matrix1->index2_data[bn];
            unsigned int inbr = A->matrix1->blocksize0[cr];
            unsigned int roffset = 0;
            unsigned int coffset = 0;
            if(cr != 0)
            {
              roffset = A->matrix1->blocksize0[cr - 1];
              inbr -= roffset;
            }
            unsigned int inbc = A->matrix1->blocksize1[cc];
            if(cc != 0)
            {
              coffset = A->matrix1->blocksize1[cc - 1];
              inbc -= coffset;
            }
            for(unsigned j = 0; j < inbc; ++j)
            {
              for(unsigned i = 0; i < inbr; ++i)
              {
                CHECK_RETURN(CSparseMatrix_symmetric_entry(A->matrix2->half_triplet, i + roffset, j + coffset,
                             A->matrix1->block[bn][i + j*inbr]));
              }
            }
          }
        }
        NSM_set_version(A->matrix2, NSM_HALF_TRIPLET,
                        NM_version(A, NM_SPARSE_BLOCK));
      }
      else if(A->matrix0)
      {
        fprintf(stderr, "NM_half_triplet: conversion is not implemented");
        exit(EXIT_FAILURE);
      }
      break;
    }
    case NM_SPARSE:
    {
      switch(A->matrix2->origin)
      {
      case NSM_TRIPLET:
      {
        A->matrix2->half_triplet = NM_csc_to_half_triplet(NM_csc(A));
        NSM_set_version(A->matrix2, NSM_HALF_TRIPLET,
                        NSM_version(A->matrix2, NSM_TRIPLET));

        break;
      }
      case NSM_CSC:
      {
        assert(A->matrix2->csc);
        A->matrix2->half_triplet = NM_csc_to_half_triplet(A->matrix2->csc);
        NSM_set_version(A->matrix2, NSM_HALF_TRIPLET, NSM_version(A->matrix2,
                                                                  NSM_CSC));

        break;
      }
      case NSM_CSR:
      {
        fprintf(stderr, "NM_half_triplet: conversion is not implemented");
        exit(EXIT_FAILURE);
        break;
      }
      default:
      case NSM_UNKNOWN:
      {
        NSM_UNKNOWN_ERR("NM_half_triplet", (int) A->matrix2->origin);
        exit(EXIT_FAILURE);
      }
      }
      break;
      default:
      {
        fprintf(stderr, "NM_half_triplet: unknown matrix type\n");
        exit(EXIT_FAILURE);
      }
    }

    }
  }
  assert(A->matrix2->half_triplet);

  return A->matrix2->half_triplet;
}

CSparseMatrix* NM_csc(NumericsMatrix *A)
{
  DEBUG_BEGIN("NM_csc(NumericsMatrix *A)\n");
  assert(A);

  if(!numericsSparseMatrix(A)->csc)
  {
    assert(A->matrix2);
    switch(A->matrix2->origin)
    {
    case NSM_TRIPLET:
    case NSM_UNKNOWN:
    {
      /*  triplet -> csc with allocation */
      A->matrix2->csc = cs_compress(NM_triplet(A));
      NSM_set_version(A->matrix2, NSM_CSC,
                      NSM_version(A->matrix2, NSM_TRIPLET));
      break;
    }
    case NSM_CSR:
    {
      A->matrix2->csc = NM_csr_to_csc(NM_csr(A));
      NSM_set_version(A->matrix2, NSM_CSC,
                      NSM_version(A->matrix2, NSM_CSR));
      break;
    }
    case NSM_HALF_TRIPLET:
    {
      numerics_error("NM_csc", "cannot get csc from half triplet");
      break;
    }
    default:
    {
      NSM_UNKNOWN_ERR("NM_csc", A->matrix2->origin);
      exit(EXIT_FAILURE);
    }
    }

    assert(A->matrix2->csc);
    NM_clearCSCTranspose(A);
  }

  assert(A->matrix2->csc);
  assert(A->matrix2->csc->m == A->size0 && "inconsistent size of csc storage");
  assert(A->matrix2->csc->n == A->size1 && "inconsistent size of csc storage");



  DEBUG_END("NM_csc(NumericsMatrix *A)\n");

  assert(NSM_version(A->matrix2, NSM_TRIPLET) <=
         NSM_version(A->matrix2, NSM_CSC));

  assert(NSM_version(A->matrix2, NSM_CSR) <=
         NSM_version(A->matrix2, NSM_CSC));

  return A->matrix2->csc;
}

CSparseMatrix* NM_csc_trans(NumericsMatrix* A)
{
  if(!numericsSparseMatrix(A)->trans_csc)
  {
    assert(A->matrix2);
    A->matrix2->trans_csc = cs_transpose(NM_csc(A), 1); /* value = 1
                                                         * ->
                                                         * allocation */
  }

  return A->matrix2->trans_csc;
}

CSparseMatrix* NM_csr(NumericsMatrix *A)
{
  assert(A);

  if(!numericsSparseMatrix(A)->csr)
  {
    assert(A->matrix2);
    switch(A->matrix2->origin)
    {
    case NSM_TRIPLET:
    case NSM_UNKNOWN:
    {
      /*  triplet -> csr with allocation */
      A->matrix2->csr = NM_triplet_to_csr(NM_triplet(A));
      NSM_set_version(A->matrix2, NSM_CSR,
                      NSM_version(A->matrix2, NSM_TRIPLET));
      break;
    }
/* MB: there was a bug here */
/*    case NSM_CSR:*/
    case NSM_CSC:
    {
      A->matrix2->csr = NM_csc_to_csr(NM_csr(A));
      NSM_set_version(A->matrix2, NSM_CSC,
                      NSM_version(A->matrix2, NSM_CSC));
      break;
    }
    case NSM_HALF_TRIPLET:
    {
      numerics_error("NM_csr", "cannot get csr from half triplet");
      break;
    }
    default:
    {
      NSM_UNKNOWN_ERR("NM_csr", A->matrix2->origin);
      exit(EXIT_FAILURE);
    }
    }

    assert(A->matrix2->csr);
  }

  assert(NSM_version(A->matrix2, NSM_TRIPLET) <=
         NSM_version(A->matrix2, NSM_CSR));

  assert(NSM_version(A->matrix2, NSM_CSR) <=
         NSM_version(A->matrix2, NSM_CSR));

  return A->matrix2->csr;
}

/* Numerics Matrix wrapper  for y <- alpha A x + beta y */
void NM_gemv(const double alpha, NumericsMatrix* A, const double *x,
             const double beta, double *y)
{
  assert(A);
  assert(x);
  assert(y);

  switch(A->storageType)
  {
  case NM_DENSE:
  {
    cblas_dgemv(CblasColMajor, CblasNoTrans, A->size0, A->size1,
                alpha, A->matrix0, A->size0, x, 1, beta, y, 1);

    break;
  }

  case NM_SPARSE_BLOCK:
  {
    SBM_gemv(A->size1, A->size0, alpha, A->matrix1, x, beta, y);

    break;
  }

  case NM_SPARSE:
  {
    assert(A->storageType == NM_SPARSE);
    CHECK_RETURN(CSparseMatrix_aaxpby(alpha, NM_csc(A), x, beta, y));
    break;
  }
  default:
  {
    assert(0 && "NM_gemv unknown storageType");
  }
  }
}

/* Numerics Matrix wrapper  for y <- alpha trans(A) x + beta y */
void NM_tgemv(const double alpha, NumericsMatrix* A, const double *x,
              const double beta, double *y)
{
  switch(A->storageType)
  {
  case NM_DENSE:
  {
    cblas_dgemv(CblasColMajor, CblasTrans, A->size0, A->size1,
                alpha, A->matrix0, A->size0, x, 1, beta, y, 1);
    break;
  }
  case NM_SPARSE_BLOCK:
  case NM_SPARSE:
  {
    CHECK_RETURN(CSparseMatrix_aaxpby(alpha, NM_csc_trans(A), x, beta, y));
    break;
  }
  default:
  {
    assert(0 && "NM_tgemv unknown storageType");
  }
  }
}

/* Insert the submatrix B into the matrix A on the position defined in
 * (start_i, start_j) position.
 */
void NM_insert(NumericsMatrix* A, const NumericsMatrix* const B,
               const unsigned int start_i, const unsigned int start_j)
{
  DEBUG_BEGIN("NM_insert\n");

  DEBUG_EXPR(NM_display_storageType(A););
  DEBUG_EXPR(NM_display_storageType(B););

  /* validation */
  assert(A->size0 >= B->size0);
  assert(A->size1 >= B->size1);

  unsigned int end_i = start_i + B->size0;
  unsigned int end_j = start_j + B->size1;
  assert(start_i <= end_i);
  assert(start_j <= end_j);
  assert((int)end_i <= A->size0);
  assert((int)end_j <= A->size1);

  /* trivial case when size(A) == size(B) */
  if(A->size0 == B->size0 && A->size1 == B->size1)
  {
    /* version managed in NM_copy */
    NM_copy(B, A);
    DEBUG_END("NM_insert\n");
    return;
  }
  /* DEBUG_PRINTF("NM_insert -- A->storageType = %i\n", A->storageType); */
  /* check the case when A is sparse block matrix */
  switch(A->storageType)
  {
  case NM_SPARSE:
  {
    switch(A->matrix2->origin)
    {
    case NSM_TRIPLET:
    {
      break;
    }
    case NSM_CSC:
    {
      A->matrix2->triplet = NM_csc_to_triplet(A->matrix2->csc);
      NSM_set_version(A->matrix2, NSM_TRIPLET, NSM_version(A->matrix2,
                                                           NSM_CSC));
      break;
    }
    case NSM_CSR:
    {
      A->matrix2->triplet = NM_csr_to_triplet(A->matrix2->csr);
      NSM_set_version(A->matrix2, NSM_TRIPLET, NSM_version(A->matrix2,
                                                           NSM_CSR));

      break;
    }
    default:
    {
      numerics_error("NM_insert","unknown origin %d for matrix A\n", A->matrix2->origin);
    }
    }
    A->matrix2->origin = NSM_TRIPLET;
    break;
  }
  case NM_DENSE:
  {
    break;
  }
  default:
  {
    numerics_error("NM_insert", "unknown storageType %d for numerics matrix A\n", A->storageType);
  }
  }
  /* DEBUG_PRINTF("NM_insert -- B->storageType = %i\n", B->storageType); */

  /* We transform B into triplet to simplify: could be optimized */
  switch(B->storageType)
  {
  case NM_SPARSE:
  {
    switch(B->matrix2->origin)
    {
    case NSM_TRIPLET:
    {
      assert(B->matrix2->triplet);
      break;
    }
    case NSM_CSC:
    {
      B->matrix2->triplet = NM_csc_to_triplet(B->matrix2->csc);
      NSM_set_version(B->matrix2, NSM_TRIPLET, NSM_version(B->matrix2,
                                                           NSM_CSC));

      break;
    }
    case NSM_CSR:
    {
      B->matrix2->triplet = NM_csr_to_triplet(B->matrix2->csr);
      NSM_set_version(B->matrix2, NSM_TRIPLET, NSM_version(B->matrix2,
                                                           NSM_CSR));
      break;
    }
    default:
    {
      numerics_error("NM_insert","unknown origin %d for matrix B\n", B->matrix2->origin);
    }
    }

    B->matrix2->origin = NSM_TRIPLET;
    assert(NSM_max_version(B->matrix2) ==
           NSM_version(B->matrix2, B->matrix2->origin));

    CS_INT * Bi =   B->matrix2->triplet->i;
    CS_INT * Bp =   B->matrix2->triplet->p;
    double * Bx =   B->matrix2->triplet->x;
    // loop over the values of B
    for(int idx = 0 ; idx < B->matrix2->triplet->nz  ; idx++)
    {
      NM_entry(A, Bi[idx] + start_i, Bp[idx] + start_j, Bx[idx]);
    }
    break;
  }
  case NM_SPARSE_BLOCK:
  case NM_DENSE:
  {
    /* could be optimized */
    double val;
    for(unsigned int i = start_i; i < end_i; ++i)
    {
      for(unsigned int j = start_j; j < end_j; ++j)
      {
        val = NM_get_value(B, i - start_i, j - start_j);
        NM_entry(A, i, j, val);
      }
    }
    break;
  }
  default:
  {
    numerics_error("NM_insert","unknown storageType %d for numerics matrix B\n", B->storageType);
  }
  }
  DEBUG_END("NM_insert\n");
  return;
}


NumericsMatrix * NM_multiply(NumericsMatrix* A, NumericsMatrix* B)
{
  DEBUG_BEGIN("NM_multiply(...) \n")
  NM_types storageType;

  NumericsMatrix * C = NM_new();

  /* At the time of writing, we are able to transform anything into NM_SPARSE,
   * hence we use this format whenever possible */
  if(A->storageType == NM_SPARSE || B->storageType == NM_SPARSE || C->storageType == NM_SPARSE)
  {
    storageType = NM_SPARSE;
  }
  else
  {
    storageType = A->storageType;
  }
  switch(storageType)
  {
  case NM_DENSE:
  {
    assert(A->matrix0);
    assert(B->matrix0);

    C->size0 = A->size0;
    C->size1 = B->size1;

    assert(!C->matrix0);
    C->matrix0 = (double *)malloc(C->size0*C->size1*sizeof(double));
    assert(C->matrix0);

    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, A->size0, B->size1, B->size0,
                1.0, A->matrix0, A->size0, B->matrix0, B->size0, 0.0, C->matrix0, A->size0);
    NM_clearSparseBlock(C);
    NM_clearSparseStorage(C);
    C->storageType=storageType;
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    assert(A->matrix1);
    assert(B->matrix1);

    /* New version taht follows the principle of sparse matrices in Csparse*/
    SparseBlockStructuredMatrix * C_SBM = SBM_multiply(A->matrix1, B->matrix1);
    NM_clearSparseBlock(C);
    NM_clearDense(C);
    NM_clearSparseStorage(C);
    C->matrix1 = C_SBM;
    C->size0 = A->size0;
    C->size1 = B->size1;
    C->storageType=storageType;
    break;
  }
  case NM_SPARSE:
  {
    /* We need to free the allocated data here, hence we have to save the
     * matrix pointer. Otherwise, we have a memleak */
#ifdef WITH_MKL_SPBLAS
    if(check_mkl_lib() && (fabs(beta -1.) < 100*DBL_EPSILON))
    {
      if(!B->matrix2) NM_triplet(B);
      NumericsSparseMatrix* result = NM_MKL_spblas_gemm(0, A->matrix2, B->matrix2);
      assert(result);
      int size0 = C->size0;
      int size1 = C->size1;
      NM_clear(C);
      NM_null(C);
      NM_fill(C, NM_SPARSE, size0, size1, result);
      NM_MKL_to_sparse_matrix(C);
      return;
    }
#endif
    DEBUG_EXPR(NM_display(A));
    DEBUG_EXPR(NM_display(B));
    DEBUG_EXPR(cs_print((const cs *) NM_csc(A),0););
    DEBUG_EXPR(cs_print((const cs *) NM_csc(B),0););
    assert(A->size1 == B->size0 && "NM_gemm :: A->size1 != B->size0 ");
    CSparseMatrix* C_csc = cs_multiply(NM_csc(A), NM_csc(B));
    DEBUG_EXPR(cs_print((const cs *) C_csc,0););
    assert(C_csc && "NM_gemm :: cs_multiply failed");
    NSM_fix_csc(C_csc);

    NM_clearDense(C);
    NM_clearSparseBlock(C);
    NM_clearSparseStorage(C);
    C->storageType=storageType;
    numericsSparseMatrix(C)->csc = C_csc;
    C->size0 = (int)C->matrix2->csc->m;
    C->size1 = (int)C->matrix2->csc->n;
    numericsSparseMatrix(C)->origin = NSM_CSC;
    break;
  }
  default:
  {
    assert(0 && "NM_multiply unknown storageType");
  }
  }

  NM_MPI_copy(A, C);
  NM_MUMPS_copy(A, C);

  if (B->storageType == NM_SPARSE)
  {
    /* anything * sparse -> sparse */
    NM_version_copy(B, C);
  }
  else
  {
    NM_version_copy(A, C);
  }
  return C;
  DEBUG_END("NM_multiply(...) \n")
}

void NM_gemm(const double alpha, NumericsMatrix* A, NumericsMatrix* B,
             const double beta, NumericsMatrix* C)
{
  NM_types storageType;

  /* At the time of writing, we are able to transform anything into NM_SPARSE,
   * hence we use this format whenever possible */
  if(A->storageType == NM_SPARSE || B->storageType == NM_SPARSE || C->storageType == NM_SPARSE)
  {
    storageType = NM_SPARSE;
  }
  else
  {
    storageType = A->storageType;
  }
  switch(storageType)
  {
  case NM_DENSE:
  {
    assert(A->matrix0);
    assert(B->matrix0);
    assert(C->matrix0);
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, A->size0, B->size1, B->size0,
                alpha, A->matrix0, A->size0, B->matrix0, B->size0, beta, C->matrix0, A->size0);

    NM_clearSparseBlock(C);
    NM_clearSparseStorage(C);
    C->storageType=storageType;
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    assert(A->matrix1);
    assert(B->matrix1);
    assert(C->matrix1);

    /* old version that cannot work for beta != 0.0*/
    /* SBM_alloc_for_gemm(A->matrix1, B->matrix1, C->matrix1); */
    /* SBM_gemm(alpha, A->matrix1, B->matrix1, beta, C->matrix1); */
    /* NM_clearDense(C); */
    /* NM_clearSparseStorage(C); */
    /* C->storageType=storageType; */

    /* New version that follows the principle of sparse matrices in Csparse*/
    SparseBlockStructuredMatrix * C_tmp = SBM_multiply(A->matrix1, B->matrix1);
    SparseBlockStructuredMatrix * result = SBM_add(C_tmp, C->matrix1, alpha, beta);
    NM_clearSparseBlock(C);
    NM_clearDense(C);
    NM_clearSparseStorage(C);
    C->matrix1 = result;
    C->size0 = A->size0;
    C->size1 = B->size1;
    C->storageType=storageType;
    break;
  }
  case NM_SPARSE:
  {
    /* We need to free the allocated data here, hence we have to save the
     * matrix pointer. Otherwise, we have a memleak */
#ifdef WITH_MKL_SPBLAS
    if(check_mkl_lib() && (fabs(beta -1.) < 100*DBL_EPSILON))
    {
      if(!B->matrix2) NM_triplet(B);
      if(!C->matrix2) NM_triplet(C);
      NumericsSparseMatrix* tmp_matrix = NM_MKL_spblas_gemm(0, A->matrix2, B->matrix2);
      assert(tmp_matrix);
      NumericsSparseMatrix* result = NM_MKL_spblas_add(0, alpha, tmp_matrix, C->matrix2);
      int size0 = C->size0;
      int size1 = C->size1;
      NM_clear(C);
      NM_null(C);
      NM_fill(C, NM_SPARSE, size0, size1, result);
      NM_MKL_to_sparse_matrix(C);
      return;
    }
#endif
    DEBUG_EXPR(NM_display(A));
    DEBUG_EXPR(NM_display(B));
    DEBUG_EXPR(cs_print((const cs *) NM_csc(A),0););
    DEBUG_EXPR(cs_print((const cs *) NM_csc(B),0););
    assert(A->size1 == B->size0 && "NM_gemm :: A->size1 != B->size0 ");
    CSparseMatrix* tmp_matrix = cs_multiply(NM_csc(A), NM_csc(B));
    DEBUG_EXPR(cs_print((const cs *) tmp_matrix,0););
    assert(tmp_matrix && "NM_gemm :: cs_multiply failed");
    NSM_fix_csc(tmp_matrix);

    CSparseMatrix* result = cs_add(tmp_matrix, NM_csc(C), alpha, beta);
    assert(result && "NM_gemm :: cs_add failed");
    NSM_fix_csc(result);

    cs_spfree(tmp_matrix);
    NM_clearDense(C);
    NM_clearSparseBlock(C);
    NM_clearSparseStorage(C);
    C->storageType=storageType;
    numericsSparseMatrix(C)->csc = result;
    C->size0 = (int)C->matrix2->csc->m;
    C->size1 = (int)C->matrix2->csc->n;
    numericsSparseMatrix(C)->origin = NSM_CSC;
    break;
  }
  default:
  {
    assert(0 && "NM_gemm unknown storageType");
  }
  }

  NM_MPI_copy(A, C);
  NM_MUMPS_copy(A, C);

}

NumericsMatrixInternalData* NM_internalData(NumericsMatrix* A)
{
  if(!A->internalData)
  {
    NM_internalData_new(A);
  }
  return A->internalData;
}




void* NM_iWork(NumericsMatrix* A, size_t size, size_t sizeof_elt)
{
  size_t bit_size = size * sizeof_elt;
  NM_internalData(A)->sizeof_elt =   sizeof_elt;

  if(!NM_internalData(A)->iWork)
  {
    assert(A->internalData->iWorkSize == 0);
    A->internalData->iWork = malloc(bit_size);
    A->internalData->iWorkSize = bit_size;

  }
  else
  {
    assert(A->internalData);

    if(bit_size > A->internalData->iWorkSize)
    {
      A->internalData->iWork = realloc(A->internalData->iWork, bit_size);
      A->internalData->iWorkSize = bit_size;
    }
  }

  assert(A->internalData->iWork);
  assert(A->internalData->iWorkSize >= bit_size);

  return A->internalData->iWork;
}

double* NM_dWork(NumericsMatrix* A, int size)
{
  if(!NM_internalData(A)->dWork)
  {
    assert(A->internalData);

    assert(A->internalData->dWorkSize == 0);
    A->internalData->dWork = (double *) malloc(size * sizeof(double));
    A->internalData->dWorkSize = size;
  }
  else
  {
    assert(A->internalData);

    if((size_t)
        size > A->internalData->dWorkSize)
    {
      A->internalData->dWork = (double *) realloc(A->internalData->dWork, size * sizeof(double));
      A->internalData->dWorkSize = size;
    }
  }

  assert(A->internalData->dWork);
  assert(A->internalData->dWorkSize >= (size_t)size);

  return A->internalData->dWork;
}

int NM_LU_factorize(NumericsMatrix* Ao)
{
  DEBUG_BEGIN(" NM_LU_factorize(NumericsMatrix* Ao)\n");
  lapack_int info = 0;
  assert(Ao->destructible); /* by default Ao->destructible == Ao */
  NumericsMatrix* A = Ao->destructible;

  if (!NM_LU_factorized(Ao))
  {
    DEBUG_PRINT("NM_LU_factorize. The Factorization is performed\n");
#ifdef FACTORIZATION_DEBUG
    if (NM_internalData(Ao)->values_sha1_count > 0)
    {
      if(NM_check_values_sha1(Ao))
      {
        numerics_warning("NM_LU_factorize", "an attempt to  factorize this matrix has already been done");
        info = 1;
        return info;
      }
    }

    NM_set_values_sha1(Ao);
#endif

    switch (A->storageType)
    {
    case NM_DENSE:
    {
      assert(A->matrix0);

      numerics_printf_verbose(2,"NM_LU_factorize, using LAPACK (DGETRF)" );

      lapack_int* ipiv = (lapack_int*)NM_iWork(A, A->size0, sizeof(lapack_int));
      DEBUG_PRINTF("iwork are initialized with size %i and %i\n",A->size0*A->size1,A->size0 );

      numerics_printf_verbose(2,"NM_LU_factorize, we compute factors and keep them in place" );
      DEBUG_PRINT("Start to call DGETRF for NM_DENSE storage\n");
      //cblas_dcopy_msan(A->size0*A->size1, A->matrix0, 1, wkspace, 1);
      DGETRF(A->size0, A->size1, A->matrix0, A->size0, ipiv, &info);
      DEBUG_PRINT("end of call DGETRF for NM_DENSE storage\n");
      if (info > 0)
      {
        fprintf(stderr,"NM_LU_factorize: LU factorisation DGETRF failed. The %d-th diagonal element is 0\n", info);
      }
      else if (info < 0)
      {
        fprintf(stderr, "NM_LU_factorize: LU factorisation DGETRF failed. The %d-th argument has an illegal value, stopping\n", -info);
      }
      if (info)
      {
        NM_internalData_free(A);
        assert(!NM_internalData(A)->isLUfactorized);
      }
    }
    break;
    case NM_SPARSE_BLOCK: /* sparse block -> triplet -> csc */
    case NM_SPARSE:
    {
      NSM_linear_solver_params* p = NSM_linearSolverParams(A);
      assert(!NM_internalData(A)->isLUfactorized);
      switch (p->solver)
      {
      case NSM_CSPARSE:
        numerics_printf_verbose(2, "NM_LU_factorize, using CSparse");

        if (!p->linear_solver_data)
        {
          assert(!NSM_linear_solver_data(p));
          assert(!p->solver_free_hook);

          p->solver_free_hook = &NSM_clear_p;
        };

        CSparseMatrix_factors* cs_lu_A = (CSparseMatrix_factors*) malloc(sizeof(CSparseMatrix_factors));
        numerics_printf_verbose(2,"NM_LU_factorize, we compute factors and keep them" );
        //DEBUG_EXPR(cs_print(NM_csc(A),0));
        info = !CSparseMatrix_lu_factorization(1, NM_csc(A), DBL_EPSILON, cs_lu_A);
        if (info)
        {
          numerics_printf_verbose(2, "NM_LU_factorize: csparse factorization failed.");
        }
        assert(!p->linear_solver_data);
        p->linear_solver_data = cs_lu_A;
        break;
#ifdef WITH_MUMPS
      case NSM_MUMPS:
      {
        if(verbose >= 2)
        {
          printf("NM_LU_factorize: using MUMPS\n");
        }
        if(!NM_MUMPS_id(A)->job || (NM_MUMPS_id(A)->job == -2))
        {
          /* the mumps instance is initialized (call with job=-1) */
          NM_MUMPS_set_control_params(A);
          NM_MUMPS(A, -1);
          if ((NM_MUMPS_icntl(A, 1) == -1 ||
               NM_MUMPS_icntl(A, 2) == -1 ||
               NM_MUMPS_icntl(A, 3) == -1 ||
               verbose) ||
              (NM_MUMPS_icntl(A, 1) != -1 ||
               NM_MUMPS_icntl(A, 2) != -1 ||
               NM_MUMPS_icntl(A, 3) != -1 ||
               !verbose))
          {
            NM_MUMPS_set_verbosity(A, verbose);
          }
          NM_MUMPS_set_icntl(A, 24, 1); // Null pivot row detection
          NM_MUMPS_set_cntl(A, 5, 1.e20); // Fixation, recommended value
        }

        NM_MUMPS_set_matrix(A);

        NM_MUMPS(A, 4); /* analyzis,factorization */

        DMUMPS_STRUC_C* mumps_id = NM_MUMPS_id(A);

        info = mumps_id->info[0];

        /* MUMPS can return info codes with negative value */
        if(info)
        {
          if(verbose > 0)
          {
            fprintf(stderr,"NM_LU_factorize: MUMPS fails : info(1)=%d, info(2)=%d\n", info, mumps_id->info[1]);
          }
        }

        /* we should not do that here */
        if(!p->solver_free_hook)
        {
          p->solver_free_hook = &NM_MUMPS_free;
        }
        break;
      }
#endif /* WITH_MUMPS */
      default:
      {
        numerics_printf_verbose(0,"NM_LU_factorize, Unknown solver in NM_SPARSE case." );
        info = 1;
      }
      }
      break;
    }
    default:
      assert (0 && "NM_LU_factors: unknown storageType");
    }

    if (!info)
    {
      NM_internalData(A)->isLUfactorized = true;
    }
    else
    {
      assert(NM_internalData(A)->isLUfactorized == false);
    }
  }

  assert (NM_LU_factorized(Ao) == NM_LU_factorized(A));
  DEBUG_END(" NM_LU_factorize(NumericsMatrix* Ao)\n");
  return info;
}

int NM_LU_solve(NumericsMatrix* Ao, double *b, unsigned int nrhs)
{

  lapack_int info = 1;

  /* factorization is done on destructible part only if
   * !A->internalData->isLUfactorized */
  NM_LU_factorize(Ao);

  /* get the destructible part of the matrix */
  NumericsMatrix *A = Ao->destructible;

  if (NM_LU_factorized(A))
  {

    DEBUG_BEGIN("NM_LU_solve(NumericsMatrix* A, double *b, unsigned int nrhs)\n");
    assert(A->size0 == A->size1);

    switch (A->storageType)
    {
    case NM_DENSE:
    {
      assert(A->matrix0);

      numerics_printf_verbose(2, "NM_LU_solve, using LAPACK" );

      numerics_printf_verbose(2, "NM_LU_solve, factorization in-place" );

      /* dgetrf is called in NM_LU_factorize */
      DEBUG_PRINT("Start to call DGETRS for NM_DENSE storage\n");

      numerics_printf_verbose(2,"NM_LU_solve, we solve with given factors" );
      lapack_int* ipiv = (lapack_int*)NM_iWork(A, A->size0, sizeof(lapack_int));

      DGETRS(LA_NOTRANS, A->size0, nrhs, A->matrix0, A->size0, ipiv, b, A->size0, &info);

      DEBUG_PRINT("End of call DGETRS for NM_DENSE storage\n");

      if (info < 0)
      {
        numerics_printf_verbose(2,"NM_LU_solve: dense LU solve DGETRS failed. The %d-th argument has an illegal value\n", -info);
      }
      break;
    }

    case NM_SPARSE_BLOCK: /* sparse block -> triplet -> csc */
    case NM_SPARSE:
    {
      NSM_linear_solver_params* p = NSM_linearSolverParams(A);
      switch (p->solver)
      {
      case NSM_CSPARSE:
      {

        if (!p->dWork)
        {
          assert(!NSM_workspace(p));
          p->dWork = (double*) malloc(A->size1 * sizeof(double));
          p->dWorkSize = A->size1;
        };


        numerics_printf_verbose(2,"NM_LU_solve, using CSparse" );
        numerics_printf_verbose(2,"NM_LU_solve, we solve with given factors" );
        for(unsigned int j=0; j < nrhs ; j++ )
        {
          info = !CSparseMatrix_solve((CSparseMatrix_factors *)NSM_linear_solver_data(p), NSM_workspace(p), &b[j*A->size1]);
        }
        if (info < 0)
        {
          numerics_printf_verbose(2,"NM_LU_solve: Csparse solver failed with info = %i \n", info);
        }
        else
          numerics_printf_verbose(2,"NM_LU_solve: Csparse  with info = %i \n", info);
        break;
      }
#ifdef WITH_MUMPS
      case NSM_MUMPS:
      {
        if(verbose >= 2)
        {
          printf("NM_LU_solve: using MUMPS\n");
        }

        assert (NM_MUMPS_id(A)->job); /* this means that least a
                                       * factorization has already been
                                       * done */

        DMUMPS_STRUC_C* mumps_id = NM_MUMPS_id(A);

        NM_MUMPS_set_dense_rhs(A, nrhs, b);

        NM_MUMPS(A, 3); /* solve */
        info = mumps_id->info[0];

        /* MUMPS can return info codes with negative value */
        if(info)
        {
          if(verbose > 0)
          {
            fprintf(stderr,"NM_LU_solve: MUMPS fails : info(1)=%d, info(2)=%d\n", info, mumps_id->info[1]);
          }
        }
        break;
      }
#endif /* WITH_MUMPS */
      default:
      {
        fprintf(stderr, "NM_LU_solve: unknown sparse linearsolver %d\n", p->solver);
        exit(EXIT_FAILURE);
      }
      break;
      }
      break;
    }
    default:
      assert (0 && "NM_LU_solve unknown storageType");
    }


    /* WARNING: cs returns 0 (false) for failed and 1 (true) for ok
       CHECK_RETURN is ok for cs, but not for MUMPS and others */
    /* some time we cannot find a solution to a linear system, and its fine, for
     * instance with the minFBLSA. Therefore, we should not check here for
     * problems, but the calling function has to check the return code.*/
//  CHECK_RETURN(info);
    DEBUG_END("NM_LU_solve(NumericsMatrix* A, double *b, unsigned keep)\n");
  }

  return info;
}
int NM_LU_solve_matrix_rhs(NumericsMatrix* Ao, NumericsMatrix* B)
{

  lapack_int info = 1;

  /* factorization is done on destructible part only if
   * !A->internalData->isLUfactorized */
  NM_LU_factorize(Ao);

  /* get the destructible part of the matrix */
  NumericsMatrix *A = Ao->destructible;

  if (NM_LU_factorized(A))
  {

    DEBUG_BEGIN("NM_LU_solve(NumericsMatrix* A, double *b, unsigned int nrhs)\n");
    assert(A->size0 == A->size1);

    if (B->storageType == NM_DENSE)
    {
      assert(B->matrix0);
      info = NM_LU_solve(A, B->matrix0, B->size1);
    }
    else if ((B->storageType == NM_SPARSE) || (B->storageType == NM_SPARSE_BLOCK))
    {
      switch (A->storageType)
      {
      case NM_DENSE:
      {
        numerics_error("NM_LU_solve_matrix_rhs", "Solving Linear system with a dense matrix and a sparse matrix rhs is not implemented since it requires to copy the rhs into dense matrix." );
        break;
      }

      case NM_SPARSE_BLOCK: /* sparse block -> triplet -> csc */
      case NM_SPARSE:
      {
        NSM_linear_solver_params* p = NSM_linearSolverParams(A);
        switch (p->solver)
        {
        case NSM_CSPARSE:
        {
          numerics_printf_verbose(2,"NM_LU_solve_matrix_rhs, using CSparse" );
          numerics_printf_verbose(2,"NM_LU_solve_matrix_rhs, we solve with given factors" );

          if (!p->dWork)
          {
            assert(!NSM_workspace(p));
            p->dWork = (double*) malloc(A->size1 * sizeof(double));
            p->dWorkSize = A->size1;
          };


          CSparseMatrix *X = cs_spalloc(NM_csc(B)->m, NM_csc(B)->n, NM_csc(B)->nzmax, 1,0); /* csc format */

          info = !CSparseMatrix_spsolve((CSparseMatrix_factors *)NSM_linear_solver_data(p), X, NM_csc(B));
          //Invalidation
          B->matrix2->origin = NSM_CSC;
          NM_clearCSR(B);
          NM_clearCSCTranspose(B);
          NM_clearTriplet(B);
          NM_clearHalfTriplet(B);

          break;
        }
#ifdef WITH_MUMPS
        case NSM_MUMPS:
        {
          numerics_printf_verbose(2,"NM_LU_solve: using MUMPS\n");

          assert (NM_MUMPS_id(A)->job); /* this means that at least a
                                         * factorization has already been
                                         * done */

          DMUMPS_STRUC_C* mumps_id = NM_MUMPS_id(A);

          NM_MUMPS_set_sparse_rhs(A, B);

          NM_MUMPS(A, 3); /* solve */
          info = mumps_id->info[0];

          /* MUMPS can return info codes with negative value */
          if(info)
          {
            if(verbose > 0)
            {
              fprintf(stderr,"NM_LU_solve: MUMPS fails : info(1)=%d, info(2)=%d\n", info, mumps_id->info[1]);
            }
          }

          /* solution is returned in the DENSE format in the B matrix */
          NM_clear(B);
          B->storageType = NM_DENSE;
          B->size0 = A->size0;
          B->size1 = mumps_id->nrhs;
          B->matrix0 = (double *) malloc(
            A->size0*mumps_id->nrhs*sizeof(double));
          memcpy(B->matrix0, mumps_id->rhs,
                 A->size0*mumps_id->nrhs*sizeof(double));
          break;
        }
#endif /* WITH_MUMPS */
        default:
        {
          fprintf(stderr, "NM_LU_solve: unknown sparse linearsolver %d\n", p->solver);
          exit(EXIT_FAILURE);
        }
        break;
        }
        break;
      }
      default:
        assert (0 && "NM_LU_solve unknown storageType");
      }


      /* WARNING: cs returns 0 (false) for failed and 1 (true) for ok
         CHECK_RETURN is ok for cs, but not for MUMPS and others */
      /* some time we cannot find a solution to a linear system, and its fine, for
       * instance with the minFBLSA. Therefore, we should not check here for
       * problems, but the calling function has to check the return code.*/

    }
  }
  DEBUG_END("NM_LU_solve(NumericsMatrix* A, double *b, unsigned keep)\n");
  return info;
}

NumericsMatrix* NM_LU_inv(NumericsMatrix* A)
{

  DEBUG_BEGIN("NM_LU_inv(NumericsMatrix* A, double *b, unsigned keep)\n");
  assert(A->size0 == A->size1);
  double * b = (double *) malloc(A->size0*sizeof(double));
  for(int i = 0; i < A->size0; ++i)
  {
    b[i]=0.0;
  }


  NumericsMatrix* Atmp = NM_new();
  NM_copy(A,Atmp);

  NumericsMatrix * Ainv  = NM_new();
  Ainv->size0 =  A->size0;
  Ainv->size1 =  A->size1;

  int info =-1;

  switch(A->storageType)
  {
  case NM_DENSE:
  {
    Ainv->storageType = NM_DENSE;
    Ainv->matrix0 = (double *)malloc(A->size0*A->size1*sizeof(double));
    for(int col_rhs =0; col_rhs < A->size1; col_rhs++)
    {
      for(int i = 0; i < A->size0; ++i)
      {
        b[i]=0.0;
      }
      b[col_rhs] = 1.0;
      DEBUG_EXPR(NV_display(b,A->size1););
      info = NM_LU_solve(Atmp, b, 1);
      DEBUG_EXPR(NV_display(b,A->size1););
      if(info)
      {
        numerics_warning("NM_LU_inv", "problem in NM_LU_solve");
      }
      for(int i = 0; i < A->size0; ++i)
      {
        Ainv->matrix0[i+col_rhs*A->size0]  = b[i];
      }
    }
    break;
  }
  case NM_SPARSE_BLOCK: /* sparse block -> triplet -> csc */
  case NM_SPARSE:
  {

    Ainv->storageType = NM_SPARSE;
    NM_triplet_alloc(Ainv,  A->size0);
    Ainv->matrix2->origin = NSM_TRIPLET;

    for(int col_rhs =0; col_rhs < A->size1; col_rhs++)
    {
      for(int i = 0; i < A->size0; ++i)
      {
        b[i]=0.0;
      }
      b[col_rhs] = 1.0;
      DEBUG_EXPR(NV_display(b,A->size1););
      info = NM_LU_solve(Atmp, b, 1);
      if(info)
      {
        numerics_warning("NM_LU_inv", "problem in NM_LU_solve");
      }
      for(int i = 0; i < A->size0; ++i)
      {
        CHECK_RETURN(CSparseMatrix_entry(Ainv->matrix2->triplet, i, col_rhs, b[i]));
      }
    }
    break;
  }
  default:
    assert(0 && "NM_LU_inv :  unknown storageType");
  }

  NM_clear(Atmp);
  free(Atmp);
  free(b);
  DEBUG_END("NM_LU_inv(NumericsMatrix* A, double *b, unsigned keep)\n");
  return Ainv;

}

int NM_gesv_expert(NumericsMatrix* A, double *b, unsigned keep)
{

  DEBUG_BEGIN("NM_gesv_expert(NumericsMatrix* A, double *b, unsigned keep)\n");
  assert(A->size0 == A->size1);

  lapack_int info = 1;

  switch(A->storageType)
  {
  case NM_DENSE:
  {
    assert(A->matrix0);

    if(keep == NM_KEEP_FACTORS)
    {
      numerics_printf_verbose(2,"NM_gesv_expert, using LAPACK");
      //double* wkspace = NM_dWork(A, A->size0*A->size1);
      lapack_int* ipiv = (lapack_int*)NM_iWork(A, A->size0, sizeof(lapack_int));
      DEBUG_PRINTF("iwork and dwork are initialized with size %i and %i\n",A->size0*A->size1,A->size0);

      if(!NM_internalData(A)->isLUfactorized)
      {
        numerics_printf_verbose(2,"NM_gesv_expert, we compute factors and keep them");
        DEBUG_PRINT("Start to call DGETRF for NM_DENSE storage\n");
        //cblas_dcopy_msan(A->size0*A->size1, A->matrix0, 1, wkspace, 1);
        DGETRF(A->size0, A->size1, A->matrix0, A->size0, ipiv, &info);
        DEBUG_PRINT("end of call DGETRF for NM_DENSE storage\n");
        if(info > 0)
        {
          if(verbose >= 2)
          {
            printf("NM_gesv: LU factorisation DGETRF failed. The %d-th diagonal element is 0\n", info);
          }
        }
        else if(info < 0)
        {
          fprintf(stderr, "NM_gesv: LU factorisation DGETRF failed. The %d-th argument has an illegal value, stopping\n", -info);
        }
        else
        {
          if(verbose >= 2)
          {
            printf("NM_gesv: LU factorisation DGETRF suceeded.\n");
          }
        }
        if(info)
        {
          NM_internalData_free(A);
          return info;
        }

        NM_internalData(A)->isLUfactorized = true;
      }
      DEBUG_PRINT("Start to call DGETRS for NM_DENSE storage\n");
      numerics_printf_verbose(2,"NM_gesv_expert, we solve with given factors");
      DGETRS(LA_NOTRANS, A->size0, 1, A->matrix0, A->size0, ipiv, b, A->size0, &info);
      DEBUG_PRINT("End of call DGETRS for NM_DENSE storage\n");
      if(info < 0)
      {
        if(verbose >= 2)
        {
          printf("NM_gesv: dense LU solve DGETRS failed. The %d-th argument has an illegal value, stopping\n", -info);
        }
      }
    }
    else
    {
      double* mat;
      if(keep == NM_PRESERVE)
      {
        mat = NM_dWork(A, A->size0*A->size1);
        cblas_dcopy_msan(A->size0*A->size1, A->matrix0, 1, mat, 1);
      }
      else
      {
        mat = A->matrix0;
      }
      DGESV(A->size0, 1, mat, A->size0, (lapack_int*)NM_iWork(A, A->size0, sizeof(lapack_int)), b,
            A->size0, &info);
    }
    break;
  }

  case NM_SPARSE_BLOCK: /* sparse block -> triplet -> csc */
  case NM_SPARSE:
  {
    NSM_linear_solver_params* p = NSM_linearSolverParams(A);
    switch(p->solver)
    {
    case NSM_CSPARSE:
      numerics_printf_verbose(2,"NM_gesv, using CSparse");

      if(keep == NM_KEEP_FACTORS)
      {
        if(!(p->dWork && p->linear_solver_data))
        {
          assert(!NSM_workspace(p));
          assert(!NSM_linear_solver_data(p));
          assert(!p->solver_free_hook);

          p->solver_free_hook = &NSM_clear_p;
          p->dWork = (double*) malloc(A->size1 * sizeof(double));
          p->dWorkSize = A->size1;
          CSparseMatrix_factors* cs_lu_A = (CSparseMatrix_factors*) malloc(sizeof(CSparseMatrix_factors));
          numerics_printf_verbose(2,"NM_gesv_expert, we compute factors and keep them");
          CHECK_RETURN(CSparseMatrix_lu_factorization(1, NM_csc(A), DBL_EPSILON, cs_lu_A));
          p->linear_solver_data = cs_lu_A;
        }

        numerics_printf_verbose(2,"NM_gesv, we solve with given factors");
        info = !CSparseMatrix_solve((CSparseMatrix_factors *)NSM_linear_solver_data(p), NSM_workspace(p), b);
      }
      else
      {
        info = !cs_lusol(1, NM_csc(A), b, DBL_EPSILON);
      }
      break;

#ifdef WITH_MUMPS
    case NSM_MUMPS:
    {
      if(verbose >= 2)
      {
        printf("NM_gesv: using MUMPS\n");
      }
      if(!NM_MUMPS_id(A)->job || (NM_MUMPS_id(A)->job == -2))
      {
        /* the mumps instance is initialized (call with job=-1) */
        NM_MUMPS_set_control_params(A);
        NM_MUMPS(A, -1);
        if((NM_MUMPS_icntl(A, 1) == -1 ||
            NM_MUMPS_icntl(A, 2) == -1 ||
            NM_MUMPS_icntl(A, 3) == -1 ||
            verbose) ||
            (NM_MUMPS_icntl(A, 1) != -1 ||
             NM_MUMPS_icntl(A, 2) != -1 ||
             NM_MUMPS_icntl(A, 3) != -1 ||
             !verbose))
        {
          NM_MUMPS_set_verbosity(A, verbose);
        }
        NM_MUMPS_set_icntl(A, 24, 1); // Null pivot row detection
        NM_MUMPS_set_cntl(A, 5, 1.e20); // Fixation, recommended value
      }
      NM_MUMPS_set_matrix(A);
      NM_MUMPS_set_dense_rhs(A, 1, b);

      DMUMPS_STRUC_C* mumps_id = NM_MUMPS_id(A);

      if(keep != NM_KEEP_FACTORS|| mumps_id->job == -1)
      {
        NM_MUMPS(A, 6); /* analyzis,factorization,solve*/
      }
      else
      {
        NM_MUMPS(A, 3); /* solve */
      }

      info = mumps_id->info[0];

      /* MUMPS can return info codes with negative value */
      if(info)
      {
        if(verbose > 0)
        {
          fprintf(stderr,"NM_gesv: MUMPS fails : info(1)=%d, info(2)=%d\n", info, mumps_id->info[1]);
        }
      }
      if(keep != NM_KEEP_FACTORS)
      {
        NM_MUMPS(A, -2);
      }
      if(!p->solver_free_hook)
      {
        p->solver_free_hook = &NM_MUMPS_free;
      }
      break;
    }
#endif /* WITH_MUMPS */

#ifdef WITH_UMFPACK
    case NSM_UMFPACK:
    {
      if(verbose >= 2)
      {
        printf("NM_gesv: using UMFPACK\n");
      }

      NM_UMFPACK_WS* umfpack_ws = NM_UMFPACK_factorize(A);

      if(!umfpack_ws)
      {
        if(verbose > 1)
          fprintf(stderr, "NM_gesv: cannot factorize the matrix with UMFPACK\n");

        NM_UMFPACK_free(p);
        return -1;
      }

      CSparseMatrix* C = NM_csc(A);
      info = (int)UMFPACK_FN(wsolve)(UMFPACK_A, C->p, C->i, C->x, umfpack_ws->x, b, umfpack_ws->numeric, umfpack_ws->control, umfpack_ws->info, umfpack_ws->wi, umfpack_ws->wd);

      if(info)
      {
        UMFPACK_FN(report_status)(umfpack_ws->control, (CS_INT)info);
      }
      else
      {
        cblas_dcopy(C->n, umfpack_ws->x, 1, b, 1);
      }

      if(keep != NM_KEEP_FACTORS)
      {
        NM_UMFPACK_free(p);
      }
      else if(!p->solver_free_hook)
      {
        p->solver_free_hook = &NM_UMFPACK_free;
      }

      break;
    }

#endif /* WITH_UMFPACK */

#ifdef WITH_SUPERLU
    case NSM_SUPERLU:
    {
      if(verbose >= 2)
      {
        printf("NM_gesv: using SuperLU\n");
      }

      NM_SuperLU_WS* superlu_ws = NM_SuperLU_factorize(A);

      if(!superlu_ws)
      {
        if(verbose > 1)
          fprintf(stderr, "NM_gesv: cannot factorize the matrix with SuperLU\n");

        NM_SuperLU_free(p);
        return -1;
      }

      info = NM_SuperLU_solve(A, b, superlu_ws);

      if(info)
      {
        fprintf(stderr, "NM_gesv: cannot solve the system with SuperLU\n");
//        SuperLU_FN(report_status) (superlu_ws->control, (CS_INT)info);
      }

      if(keep != NM_KEEP_FACTORS)
      {
        NM_SuperLU_free(p);
      }
      else if(!p->solver_free_hook)
      {
        p->solver_free_hook = &NM_SuperLU_free;
      }

      break;
    }

#endif /* WITH_SUPERLU */

#ifdef WITH_SUPERLU_MT
    case NSM_SUPERLU_MT:
    {
      if(verbose >= 2)
      {
        printf("NM_gesv: using SuperLU_MT\n");
      }

      NM_SuperLU_MT_WS* superlu_mt_ws = NM_SuperLU_MT_factorize(A);

      if(!superlu_mt_ws)
      {
        if(verbose > 1)
          fprintf(stderr, "NM_gesv: cannot factorize the matrix with SuperLU_MT\n");

        NM_SuperLU_MT_free(p);
        return -1;
      }

      info = NM_SuperLU_MT_solve(A, b, superlu_mt_ws);

      if(info)
      {
        fprintf(stderr, "NM_gesv: cannot solve the system with SuperLU_MT\n");
//        SuperLU_MT_FN(report_status) (superlu_ws->control, (CS_INT)info);
      }

      if(keep != NM_KEEP_FACTORS)
      {
        NM_SuperLU_MT_free(p);
      }
      else if(!p->solver_free_hook)
      {
        p->solver_free_hook = &NM_SuperLU_MT_free;
      }

      break;
    }

#endif /* WITH_SUPERLU_MT */

#ifdef WITH_MKL_PARDISO
    case NSM_MKL_PARDISO:
    {
      if(verbose >= 2)
      {
        printf("NM_gesv: using MKL_PARDISO\n");
      }

      NM_MKL_pardiso_WS* pardiso_ws = NM_MKL_pardiso_factorize(A);

      if(!pardiso_ws)
      {
        if(verbose > 1)
          fprintf(stderr, "NM_gesv: cannot factorize the matrix with MKL_PARDISO\n");

        NM_MKL_pardiso_free(p);
        return -1;
      }

      info = NM_MKL_pardiso_solve(A, b, pardiso_ws);

      if(keep != NM_KEEP_FACTORS)
      {
        NM_MKL_pardiso_free(p);
      }
      else if(!p->solver_free_hook)
      {
        p->solver_free_hook = &NM_MKL_pardiso_free;
      }

      break;
    }

#endif /* WITH_MKL_pardiso */
    default:
    {
      fprintf(stderr, "NM_gesv: unknown sparse linearsolver %d\n", p->solver);
      exit(EXIT_FAILURE);
    }
    }
    break;
  }

  default:
    assert(0 && "NM_gesv unknown storageType");
  }

  /* some time we cannot find a solution to a linear system, and its fine, for
   * instance with the minFBLSA. Therefore, we should not check here for
   * problems, but the calling function has to check the return code.*/
//  CHECK_RETURN(info);
  DEBUG_END("NM_gesv_expert(NumericsMatrix* A, double *b, unsigned keep)\n");
  return (int)info;
}

int NM_posv_expert(NumericsMatrix* A, double *b, unsigned keep)
{

  DEBUG_BEGIN("NM_posv_expert(NumericsMatrix* A, double *b, unsigned keep)\n");
  assert(A->size0 == A->size1);

  lapack_int info = 1;
  /* verbose=2; */
  switch(A->storageType)
  {
  case NM_DENSE:
  {
    assert(A->matrix0);

    if(keep == NM_KEEP_FACTORS)
    {
      numerics_printf_verbose(2,"NM_posv_expert, using LAPACK");
      //double* wkspace = NM_dWork(A, A->size0*A->size1);

      DEBUG_PRINTF("iwork and dwork are initialized with size %i and %i\n",A->size0*A->size1,A->size0);

      if(!NM_internalData(A)->isLUfactorized)
      {
        numerics_printf_verbose(2,"NM_posv_expert, we compute factors and keep them");
        DEBUG_PRINT("Start to call DPOTRF for NM_DENSE storage\n");
        //cblas_dcopy_msan(A->size0*A->size1, A->matrix0, 1, wkspace, 1);
        DPOTRF(LA_UP, A->size1, A->matrix0, A->size0, &info);
        DEBUG_PRINT("end of call DPOTRF for NM_DENSE storage\n");
        if(info > 0)
        {
          if(verbose >= 2)
          {
            printf("NM_posv: Cholesky factorisation DPOTRF failed. The %d-th diagonal element is 0\n", info);
          }
        }
        else if(info < 0)
        {
          fprintf(stderr, "NM_posv: Cholesky factorisation DPOTRF failed. The %d-th argument has an illegal value, stopping\n", -info);
        }

        if(info)
        {
          NM_internalData_free(A);
          return info;
        }

        NM_internalData(A)->isLUfactorized = true;
      }
      DEBUG_PRINT("Start to call DPOTRS for NM_DENSE storage\n");
      numerics_printf_verbose(2,"NM_posv_expert, we solve with given factors");
      DEBUG_EXPR(NV_display(b,A->size0));
      DPOTRS(LA_UP, A->size0, 1, A->matrix0, A->size0, b, A->size0, &info);
      DEBUG_EXPR(NV_display(b,A->size0));
      DEBUG_EXPR(NM_display(A));
      DEBUG_PRINT("End of call DPOTRS for NM_DENSE storage\n");
      if(info < 0)
      {
        if(verbose >= 2)
        {
          printf("NM_posv: dense LU solve DGETRS failed. The %d-th argument has an illegal value, stopping\n", -info);
        }
      }
    }
    else
    {
      double* mat;
      if(keep == NM_PRESERVE)
      {
        mat = NM_dWork(A, A->size0*A->size1);
        cblas_dcopy_msan(A->size0*A->size1, A->matrix0, 1, mat, 1);
      }
      else
      {
        mat = A->matrix0;
      }
      DPOSV(LA_UP, A->size0, 1, mat, A->size0, b,
            A->size0, &info);

      if(info > 0)
      {
        if(verbose >= 2)
        {
          printf("NM_posv: Cholesky solver DPOSV failed. The %d-th diagonal element is 0\n", info);
        }
      }
      else if(info < 0)
      {
        fprintf(stderr, "NM_posv: Cholesky solver DPOSV failed. The %d-th argument has an illegal value, stopping\n", -info);
      }
    }
    break;
  }

  case NM_SPARSE_BLOCK: /* sparse block -> triplet -> csc */
  case NM_SPARSE:
  {
    NSM_linear_solver_params* p = NSM_linearSolverParams(A);
    switch(p->solver)
    {
    case NSM_CSPARSE:
      numerics_printf_verbose(2,"NM_posv, using CSparse cholsol");

      if(keep == NM_KEEP_FACTORS)
      {
        if(!(p->dWork && p->linear_solver_data))
        {
          assert(!NSM_workspace(p));
          assert(!NSM_linear_solver_data(p));
          assert(!p->solver_free_hook);

          p->solver_free_hook = &NSM_clear_p;
          p->dWork = (double*) malloc(A->size1 * sizeof(double));
          p->dWorkSize = A->size1;

          CSparseMatrix_factors* cs_chol_A = (CSparseMatrix_factors*) malloc(sizeof(CSparseMatrix_factors));

          numerics_printf_verbose(2,"NM_posv_expert, we compute factors and keep them");
          CHECK_RETURN(CSparseMatrix_chol_factorization(1, NM_csc(A),  cs_chol_A));

          p->linear_solver_data = cs_chol_A;
        }

        numerics_printf_verbose(2,"NM_posv, we solve with given factors");
        info = !CSparseMatrix_chol_solve((CSparseMatrix_factors *)NSM_linear_solver_data(p), NSM_workspace(p), b);
      }
      else
      {
        numerics_printf_verbose(2,"NM_posv");
        info = !cs_cholsol(1, NM_csc(A), b);
        if(info > 0)
        {
          numerics_printf("NM_posv: cs_cholsol failed. info = %i\n", info);
        }
        //DEBUG_EXPR(NV_display)
      }
      break;

#ifdef WITH_MUMPS
    case NSM_MUMPS:
    {
      if(verbose >= 2)
      {
        if(!NM_MUMPS_id(A)->job || (NM_MUMPS_id(A)->job == -2))
        {
          printf("NM_posv: using MUMPS\n");
        }
      }
      /* construction of lower triangular matrix */
      if(!NM_MUMPS_id(A)->job || (NM_MUMPS_id(A)->job == -2))
      {
        /* the mumps instance is initialized (call with job=-1) */
        NM_MUMPS_set_control_params(A);
        NM_MUMPS_set_sym(A, 2); /* general symmetric */
        NM_MUMPS(A, -1);
        if((NM_MUMPS_icntl(A, 1) == -1 ||
            NM_MUMPS_icntl(A, 2) == -1 ||
            NM_MUMPS_icntl(A, 3) == -1 ||
            verbose) ||
            (NM_MUMPS_icntl(A, 1) != -1 ||
             NM_MUMPS_icntl(A, 2) != -1 ||
             NM_MUMPS_icntl(A, 3) != -1 ||
             !verbose))
        {
          NM_MUMPS_set_verbosity(A, verbose);
        }

        NM_MUMPS_set_icntl(A, 24, 1); // Null pivot row detection
        NM_MUMPS_set_cntl(A, 5, 1.e20); // Fixation, recommended value
      }

      NM_MUMPS_set_matrix(A);
      NM_MUMPS_set_dense_rhs(A, 1, b);

      DMUMPS_STRUC_C* mumps_id = NM_MUMPS_id(A);

      if(keep != NM_KEEP_FACTORS|| mumps_id->job == -1)
      {
        NM_MUMPS(A, 6); /* analyzis,factorization,solve*/
      }
      else
      {
        NM_MUMPS(A, 3); /* solve */
      }

      info = mumps_id->info[0];

      /* MUMPS can return info codes with negative value */
      if(info)
      {
        if(verbose > 0)
        {
          fprintf(stderr, "NM_posv: MUMPS fails : info(1)=%d, info(2)=%d\n", info, mumps_id->info[1]);
        }
      }
      if(keep != NM_KEEP_FACTORS)
      {
        NM_MUMPS(A, -2);
      }
      if(!p->solver_free_hook)
      {
        p->solver_free_hook = &NM_MUMPS_free;
      }
      break;
    }
#endif /* WITH_MUMPS */

    default:
    {
      fprintf(stderr, "NM_posv: unknown sparse linearsolver %d\n", p->solver);
      exit(EXIT_FAILURE);
    }
    }
    break;
  }

  default:
    assert(0 && "NM_posv unknown storageType");
  }

  DEBUG_END("NM_posv_expert(NumericsMatrix* A, double *b, unsigned keep)\n");
  return (int)info;
}

int NM_gesv_expert_multiple_rhs(NumericsMatrix* A, double *b, unsigned int n_rhs, unsigned keep)
{
  assert(A->size0 == A->size1);
  int info = 0;
  for(unsigned int i = 0; i < n_rhs; ++i)
  {
    info = NM_gesv_expert(A, &b[A->size0*i],  keep);
    if(info) break;
  }
  return info;
}
NumericsMatrix* NM_gesv_inv(NumericsMatrix* A)
{

  DEBUG_BEGIN("NM_inv(NumericsMatrix* A, double *b, unsigned keep)\n");
  assert(A->size0 == A->size1);
  double * b = (double *) malloc(A->size0*sizeof(double));
  for(int i = 0; i < A->size0; ++i)
  {
    b[i]=0.0;
  }


  NumericsMatrix* Atmp = NM_new();
  NM_copy(A,Atmp);

  NumericsMatrix * Ainv  = NM_new();
  Ainv->size0 =  A->size0;
  Ainv->size1 =  A->size1;

  int info =-1;

  switch(A->storageType)
  {
  case NM_DENSE:
  {
    Ainv->storageType = NM_DENSE;
    Ainv->matrix0 = (double *)malloc(A->size0*A->size1*sizeof(double));
    for(int col_rhs =0; col_rhs < A->size1; col_rhs++)
    {
      for(int i = 0; i < A->size0; ++i)
      {
        b[i]=0.0;
      }
      b[col_rhs] = 1.0;
      DEBUG_EXPR(NV_display(b,A->size1););
      //info = NM_gesv_expert(Atmp, b, NM_PRESERVE);
      info = NM_gesv_expert(Atmp, b, NM_KEEP_FACTORS);
      DEBUG_EXPR(NV_display(b,A->size1););
      if(info)
      {
        numerics_warning("NM_inv", "problem in NM_gesv_expert");
      }
      for(int i = 0; i < A->size0; ++i)
      {
        Ainv->matrix0[i+col_rhs*A->size0]  = b[i];
      }
    }
    break;
  }
  case NM_SPARSE_BLOCK: /* sparse block -> triplet -> csc */
  case NM_SPARSE:
  {

    Ainv->storageType = NM_SPARSE;
    NM_triplet_alloc(Ainv,  A->size0);
    Ainv->matrix2->origin = NSM_TRIPLET;

    for(int col_rhs =0; col_rhs < A->size1; col_rhs++)
    {
      for(int i = 0; i < A->size0; ++i)
      {
        b[i]=0.0;
      }
      b[col_rhs] = 1.0;
      DEBUG_EXPR(NV_display(b,A->size1););
      //info = NM_gesv_expert(Atmp, b, NM_PRESERVE);
      info = NM_gesv_expert(Atmp, b, NM_KEEP_FACTORS);
      if(info)
      {
        numerics_warning("NM_inv", "problem in NM_gesv_expert");
      }
      for(int i = 0; i < A->size0; ++i)
      {
        CHECK_RETURN(CSparseMatrix_entry(Ainv->matrix2->triplet, i, col_rhs, b[i]));
      }
    }
    break;
  }
  default:
    assert(0 && "NM_inv :  unknown storageType");
  }

  NM_clear(Atmp);
  free(Atmp);
  free(b);
  DEBUG_END("NM_inv(NumericsMatrix* A, double *b, unsigned keep)\n");
  return Ainv;

}

int NM_inverse_diagonal_block_matrix_in_place(NumericsMatrix* A)
{

  DEBUG_BEGIN("NM_inverse_diagonal_block_matrix_in_place(NumericsMatrix* A)\n");
  assert(A->size0 == A->size1);
  int info =-1;

  // get internal data (allocation of needed)
  NM_internalData(A);

  switch(A->storageType)
  {
  case NM_SPARSE_BLOCK:
  {
    // get internal data (allocation of needed)
    lapack_int* ipiv = (lapack_int*)NM_iWork(A, A->size0, sizeof(lapack_int));
    assert(A->matrix1);
    info = SBM_inverse_diagonal_block_matrix_in_place(A->matrix1, ipiv);
    NM_internalData(A)->isInversed = true;
    break;
  }
  default:
    assert(0 && "NM_inverse_diagonal_block_matrix_in_place :  unknown storageType");
  }


  DEBUG_BEGIN("NM_inverse_diagonal_block_matrix_in_place(NumericsMatrix* A)\n");
  return (int)info;
}


void NM_update_size(NumericsMatrix* A)
{
  switch(A->storageType)
  {
  case NM_DENSE:
  {
    /* Can't do anything here ...  */
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    assert(A->matrix1);
    assert(A->matrix1->blocknumber0 > 0);
    assert(A->matrix1->blocknumber1 > 0);
    A->size0 = A->matrix1->blocksize0[A->matrix1->blocknumber0-1];
    A->size1 = A->matrix1->blocksize1[A->matrix1->blocknumber1-1];
    DEBUG_PRINT("NM_update_size :: to be implemented for NM_SPARSE_BLOCK");
    break;
  }
  case NM_SPARSE:
  {
    assert(A->matrix2);
    switch(A->matrix2->origin)
    {
    case NSM_CSC:
    {
      A->size0 = (int)A->matrix2->csc->m;
      A->size1 = (int)A->matrix2->csc->n;
      break;
    }
    case NSM_CSR:
    {
      A->size0 = (int)A->matrix2->csc->m;
      A->size1 = (int)A->matrix2->csc->n;
      break;
    }
    case NSM_TRIPLET:
    {
      A->size0 = (int)A->matrix2->triplet->m;
      A->size1 = (int)A->matrix2->triplet->n;
      break;
    }
    case NSM_HALF_TRIPLET:
    {
      A->size0 = (int)A->matrix2->half_triplet->m;
      A->size1 = (int)A->matrix2->half_triplet->n;
      break;
    }
    default:
    {
      assert(0 && "NM_update_size :: sparse matrice but neither csc nor triplet are != NULL");
    }
    }
    break;
  }
  default:
    DEBUG_PRINT("NM_update_size :: default case");
  }
}



void NM_csc_alloc(NumericsMatrix* A, CS_INT nzmax)
{
  numericsSparseMatrix(A)->csc = cs_spalloc(A->size0, A->size1, nzmax, 1, 0);
}
void NM_csc_empty_alloc(NumericsMatrix* A, CS_INT nzmax)
{
  NM_csc_alloc(A, nzmax);
  CS_INT* Ap = numericsSparseMatrix(A)->csc->p;
  for(int i =0 ; i < A->size1+1 ; i++)  Ap[i]= 0;
  //CS_INT* Ai = numericsSparseMatrix(A)->csc->i;
  //for (int i =0 ; i < nzmax ; i++)  Ai[i]= 0;
}




void NM_triplet_alloc(NumericsMatrix* A, CS_INT nzmax)
{
  numericsSparseMatrix(A)->triplet = cs_spalloc(A->size0, A->size1, nzmax, 1, 1);
  numericsSparseMatrix(A)->origin = NSM_TRIPLET;
}

void NM_setSparseSolver(NumericsMatrix* A, NSM_linear_solver solver_id)
{
  NSM_linearSolverParams(A)->solver = solver_id;
}



int NM_check(const NumericsMatrix* const A)
{
  int info = 0;
  if(!A->matrix2) return info;

  if(A->matrix2->csc) info = CSparseMatrix_check_csc(A->matrix2->csc);
  if(A->matrix2->csr) info = info ? info : CSparseMatrix_check_csc(A->matrix2->csr);
  if(A->matrix2->triplet) info = info ? info : CSparseMatrix_check_triplet(A->matrix2->triplet);
  return info;
}


size_t NM_nnz(const NumericsMatrix* M)
{
  switch(M->storageType)
  {
  case NM_DENSE:
    return M->size0 * M->size1;
  case NM_SPARSE:
  {
    assert(M->matrix2);
    return NSM_nnz(NSM_get_origin(M->matrix2));
  }
  default:
    numerics_error("NM_nnz", "Unsupported matrix type %d in %s", M->storageType);
    return SIZE_MAX;
  }
}

double NM_norm_1(NumericsMatrix* A)
{
  assert(A);


  switch(A->storageType)
  {
  case NM_DENSE:
  {
    assert(A->storageType == NM_DENSE);

    return cs_norm(NM_csc(A));
  }
  case NM_SPARSE_BLOCK:
  {
    assert(A->storageType == NM_SPARSE_BLOCK);

    return cs_norm(NM_csc(A));
  }
  case NM_SPARSE:
  {
    assert(A->storageType == NM_SPARSE);

    return cs_norm(NM_csc(A));
  }
  default:
  {
    assert(0 && "NM_norm unknown storageType");
  }
  }
  return NAN;
}

double NM_norm_inf(NumericsMatrix* A)
{
  assert(A);
  CSparseMatrix* Acsct = cs_transpose(NM_csc(A), 1);
  double norm = cs_norm(Acsct);
  assert(norm >=0);
  cs_spfree(Acsct);
  return norm;
}

int NM_is_symmetric(NumericsMatrix* A)
{
  assert(A);

  NumericsMatrix * Atrans = NM_transpose(A);
  NumericsMatrix * AplusATrans = NM_add(1/2.0, A, -1/2.0, Atrans);
  double norm_inf = NM_norm_inf(AplusATrans);

  NM_clear(Atrans);
  free(Atrans);
  NM_clear(AplusATrans);
  free(AplusATrans);

  if(norm_inf <= DBL_EPSILON*10)
  {
    return 1;
  }
  /* int n = A->size0; */
  /* int m = A->size1; */
  /* for (int i =0; i < n ; i++) */
  /* { */
  /*   for (int j =0 ; j < m ; j++) */
  /*   { */
  /*     if (fabs(NM_get_value(A,i,j)-NM_get_value(A,j,i)) >= DBL_EPSILON) */
  /*       return 0; */
  /*   } */
  /* }   */
  return 0;
}

double NM_symmetry_discrepancy(NumericsMatrix* A)
{
  int n = A->size0;
  int m = A->size1;

  double d = 0.0;
  for(int i =0; i < n ; i++)
  {
    for(int j =0 ; j < m ; j++)
    {
      d= fmax(d,fabs(NM_get_value(A,i,j)-NM_get_value(A,j,i)));
    }
  }
  return d;
}
#include "time.h"
double NM_iterated_power_method(NumericsMatrix* A, double tol, int itermax)
{
  int n = A->size0;
  int m = A->size1;
  DEBUG_EXPR(
    FILE* foutput = fopen("toto.py", "w");
    NM_write_in_file_python(A, foutput);
    fclose(foutput);
  );

  double eig = 0.0, eig_old = 2*tol;

  double * q = (double *) malloc(n*sizeof(double));
  double * z = (double *) malloc(n*sizeof(double));

  double * x1 = (double *) malloc(n*sizeof(double));
  double * z2= (double *) malloc(n*sizeof(double));

  double * q2 = (double *) malloc(n*sizeof(double));



  srand(time(NULL));
  for(int i = 0; i < n ; i++)
  {
    q[i] = (rand()/(double)RAND_MAX);
    DEBUG_PRINTF("q[%i] = %e \t",i, q[i]);
  }


  cblas_dcopy(n, q, 1, q2, 1);



  double norm = cblas_dnrm2(n, q, 1);
  cblas_dscal(m, 1.0/norm, q, 1);
  DEBUG_PRINT("\n");

  NM_gemv(1, A, q, 0.0, z);

  int k =0;

  double criteria = 1.0;

  while((criteria > tol) && k < itermax)
  {
    norm = cblas_dnrm2(n, z, 1);
    cblas_dscal(m, 1.0/norm, z, 1);
    cblas_dcopy(n, z, 1, q, 1);

    NM_gemv(1.0, A, q, 0.0, z);

    eig_old=eig;
    eig = cblas_ddot(n, q, 1, z, 1);


    cblas_dcopy(n, q, 1, x1, 1);
    NM_tgemv(1.0, A, q2, 0.0, z2);
    double norm2 = cblas_dnrm2(n, z2, 1);
    cblas_dscal(m, 1.0/norm2, z2, 1);
    cblas_dcopy(n, z2, 1, q2, 1);

    double costheta =  fabs(cblas_ddot(n, q2, 1, x1, 1));
    if(costheta <= 5e-02)
    {
      numerics_printf("NM_iterated_power_method. failed Multiple eigenvalue\n");
      break;
    }
    /*numerics_printf_verbose(1,"eig[%i] = %32.24e \t error = %e, costheta = %e \n",k, eig, fabs(eig-eig_old)/eig_old, costheta );*/
    k++;
    if(fabs(eig_old) > DBL_EPSILON)
      criteria = fabs((eig-eig_old)/eig_old);
    else
      criteria = fabs((eig-eig_old));
  }


  free(q);
  free(z);
  free(x1);
  free(z2);
  free(q2);

  return eig;
}

int NM_max_by_columns(NumericsMatrix *A, double * max)
{

  assert(A);
  assert(max);

  switch(A->storageType)
  {
  case NM_DENSE:
  case NM_SPARSE_BLOCK:
  case NM_SPARSE:
  {
    return CSparseMatrix_max_by_columns(NM_csc(A), max);
  }
  default:
  {
    assert(0 && "NM_max_by_columns unknown storageType");
  }
  }
  return 0;
}
int NM_max_by_rows(NumericsMatrix *A, double * max)
{

  assert(A);
  assert(max);

  CSparseMatrix* Acsct = cs_transpose(NM_csc(A), 1);
  int info  = CSparseMatrix_max_by_columns(Acsct, max);
  cs_spfree(Acsct);
  return info;
}
int NM_max_abs_by_columns(NumericsMatrix *A, double * max)
{

  assert(A);
  assert(max);

  switch(A->storageType)
  {
  case NM_DENSE:
  case NM_SPARSE_BLOCK:
  case NM_SPARSE:
  {
    return CSparseMatrix_max_abs_by_columns(NM_csc(A), max);
  }
  default:
  {
    assert(0 && "NM_max_abs_by_columns unknown storageType");
  }
  }
  return 0;
}
int NM_max_abs_by_rows(NumericsMatrix *A, double * max)
{

  assert(A);
  assert(max);

  CSparseMatrix* Acsct = cs_transpose(NM_csc(A), 1);
  int info  = CSparseMatrix_max_abs_by_columns(Acsct, max);
  cs_spfree(Acsct);
  return info;
}



BalancingMatrices * NM_BalancingMatrices_new(NumericsMatrix* A)
{
  BalancingMatrices * B = (BalancingMatrices *)malloc(sizeof(BalancingMatrices));
  B->size0 =  A->size0;
  B->size1 =  A->size1;
  B->D1 = NM_eye(B->size0);
  B->D2 = NM_eye(B->size1);
  B->A = NM_create(NM_SPARSE,A->size0,A->size1);
  NM_copy(A, B->A);
  return B;
}

BalancingMatrices * NM_BalancingMatrices_free(BalancingMatrices * B)
{
  if (B->D1)
    NM_free(B->D1);
  if (B->D2)
    NM_free(B->D2);
  if (B->A)
    NM_free(B->A);
  free(B);
  return NULL;
}

int NM_compute_balancing_matrices(NumericsMatrix* A, double tol, int itermax, BalancingMatrices * B)
{

  NumericsMatrix* D1_k = B->D1;
  NumericsMatrix* D2_k = B->D2;

  double * D1_k_x= D1_k->matrix2->triplet->x;
  double * D2_k_x= D2_k->matrix2->triplet->x;

  unsigned int size0 = B->size0;
  unsigned int size1 = B->size1;

  NumericsMatrix* D_R = NM_eye(size0);
  NumericsMatrix* D_C = NM_eye(size1);

  double * D_C_x= D_C->matrix2->triplet->x;
  double * D_R_x= D_R->matrix2->triplet->x;

  double * tmp = (double *) malloc(size0*size1 * sizeof(double));

  NumericsMatrix* D1_tmp = NM_eye(size0);
  NumericsMatrix* D2_tmp = NM_eye(size1);

  NumericsMatrix* A_k = B->A;
  NumericsMatrix* A_tmp = NM_create(NM_SPARSE, size0, size1);
  NM_copy(A, A_tmp);

  double error = tol + 1.0;
  int k =0;
  NM_max_abs_by_columns(A_k, D_C_x);
  NM_max_abs_by_rows(A_k, D_R_x);


  while((error > tol) && (k < itermax))
  {
    numerics_printf_verbose(2,"NM_compute_balancing_matrices iteration : %i ", k);

    NM_clearCSC(D_C);
    NM_clearCSC(D_R);

    /* inverse balancing matrices */
    for(unsigned int i=0 ; i < size0; i++)
    {
      D_R_x[i] =1.0/sqrt(D_R_x[i]);
    }
    for(unsigned int i=0 ; i < size1; i++)
    {
      D_C_x[i] =1.0/sqrt(D_C_x[i]);
    }

    /* Update balancing matrix */
    /* NM_gemm(1.0, D1_k, D_R, 0.0, D1_tmp); */
    /* NM_copy(D1_tmp, D1_k); */

    /* NM_gemm(1.0, D2_k, D_C, 0.0, D2_tmp); */
    /* NM_copy(D2_tmp, D2_k); */

    for(unsigned int i=0 ; i < size0; i++)
    {
      D1_k_x[i] = D1_k_x[i] * D_R_x[i];
    }
    for(unsigned int i=0 ; i < size1; i++)
    {
      D2_k_x[i] = D2_k_x[i] * D_C_x[i];
    }

    /* NM_display(D1_k); */
    /* DEBUG_PRINTF("D1_k ");NV_display(NM_triplet(D1_k)->x, size); */
    /* DEBUG_PRINTF("D2_k ");NV_display(NM_triplet(D2_k)->x, size); */
    /* printf("\n\n\n\n"); */

    /* Compute new balanced matrix */
    NM_gemm(1.0, A_k, D_C, 0.0, A_tmp);
    NM_gemm(1.0, D_R, A_tmp, 0.0, A_k);
    /* printf("D1_k ");NV_display(NM_triplet(D1_k)->x, size0); */
    /* printf("D2_k ");NV_display(NM_triplet(D2_k)->x, size1); */

    /* DEBUG_PRINTF("inv D_R ");NV_display(D_R_x, size); */
    /* DEBUG_PRINTF("inv D_C ");NV_display(D_C_x, size); */
    /* Compute error */

    NM_max_abs_by_rows(A_k, D_R_x);
    NM_max_abs_by_columns(A_k, D_C_x);
    /* printf("D_R ");NV_display(D_R_x, size0);  */
    /* printf("D_C ");NV_display(D_C_x, size1); */


    for(unsigned int i=0 ; i < size0; i++)
    {
      tmp[i] =(1.0 - D_R_x[i]);
    }
    error = fabs(NV_max(tmp,size0));

    for(unsigned int i=0 ; i < size1; i++)
    {
      tmp[i] =(1.0 - D_C_x[i]);
    }
    error += fabs(NV_max(tmp,size1));


    //error = fabs(1.0-NV_max(D_C_x,size)) +
    //  fabs(1.0-NV_max(D_R_x,size));

    numerics_printf_verbose(2,"NM_compute_balancing_matrices error =%e", error);
    /* printf("D_R ");NV_display(D_R_x, size); */
    /* printf("D_C ");NV_display(D_C_x, size); */
    /* printf("\n\n\n\n"); */

    k = k+1;
  }
  /* NM_clearCSC(D1_k); */
  /* NM_clearCSC(D2_k); */

  numerics_printf_verbose(1,"NM_compute_balancing_matrices final error =%e\n", error);
  free(tmp);
  NM_clear(D_R);
  free(D_R);
  NM_clear(D_C);
  free(D_C);
  NM_clear(D1_tmp);
  free(D1_tmp);
  NM_clear(D2_tmp);
  free(D2_tmp);
  NM_clear(A_tmp);
  free(A_tmp);

  if(error > tol)
  {
    return 0;
  }
  else
    return 1;
}

#ifdef WITH_OPENSSL
void NM_compute_values_sha1(NumericsMatrix* A, unsigned char* digest)
{
  switch(A->storageType)
  {
    /* ! A->matrix0 fortran layout != A->matrix2->triplet->x C layout */
  case NM_DENSE:
  case NM_SPARSE_BLOCK:
  {
    /* so triplet will be updated */
    NM_clearTriplet(A);
  }
  case NM_SPARSE:
  {
    SHA1((char*) NM_triplet(A)->x, NM_triplet(A)->nz*sizeof(double), digest);
    break;
  }
  default:
  {
    numerics_error("NM_compute_values_sha1",
                   "Unsupported matrix type %d in %s", A->storageType);
  }
  }
}

unsigned char* NM_values_sha1(NumericsMatrix* A)
{
  return NM_internalData(A)->values_sha1;
}

void NM_set_values_sha1(NumericsMatrix* A)
{
  NM_compute_values_sha1(A, NM_values_sha1(A));
  NM_internalData(A)->values_sha1_count += 1;
}

void NM_clear_values_sha1(NumericsMatrix* A)
{
  NM_internalData(A)->values_sha1_count = 0;
}


bool NM_check_values_sha1(NumericsMatrix* A)
{
  char current_digest[SHA_DIGEST_LENGTH];
  char* digest = NM_values_sha1(A);

  NM_compute_values_sha1(A, current_digest);

  if (memcmp(digest, current_digest, SHA_DIGEST_LENGTH) == 0)
  {
    return true;
  }
  else
  {
    return false;
  }
}

bool NM_equal_values_sha1(NumericsMatrix* A, NumericsMatrix* B)
{
  if (memcmp(NM_values_sha1(A), NM_values_sha1(B),
             SHA_DIGEST_LENGTH) == 0)
  {
    return true;
  }
  else
  {
    return false;
  }
}


#endif


int NM_Cholesky_factorize(NumericsMatrix* Ao)
{
  DEBUG_BEGIN("int NM_Cholesky_factorize(NumericsMatrix* Ao) \n");
  lapack_int info = 0;
  assert(Ao->destructible); /* by default Ao->destructible == Ao */
  NumericsMatrix* A = Ao->destructible;

  if (!NM_Cholesky_factorized(Ao))
  {

#ifdef FACTORIZATION_DEBUG
    if (NM_internalData(Ao)->values_sha1_count > 0)
    {
      if(NM_check_values_sha1(Ao))
      {
        numerics_error("NM_Cholesky_factorize", "this matrix is already factorized");
      }
    }

    NM_set_values_sha1(Ao);
#endif

    switch (A->storageType)
    {
    case NM_DENSE:
    {
      assert(A->matrix0);

      numerics_printf_verbose(2,"NM_Cholesky_factorize, using LAPACK (POTRF)" );

      DEBUG_PRINTF("iwork are initialized with size %i and %i\n",A->size0*A->size1,A->size0 );

      numerics_printf_verbose(2,"NM_Cholesky_factorize, we compute factors and keep them in place" );
      DEBUG_PRINT("Start to call DPOTRF for NM_DENSE storage\n");
      DPOTRF(LA_UP, A->size1, A->matrix0, A->size0, &info);
      DEBUG_PRINT("end of call DPOTRF for NM_DENSE storage\n");

      if (info > 0)
      {
        fprintf(stderr,"NM_Cholesky_factorize: Cholesky factorisation DPOTRF failed. The %d-th diagonal element is 0\n", info);
      }
      else if (info < 0)
      {
        fprintf(stderr, "NM_Cholesky_factorize: Cholesky factorisation DPOTRF failed. The %d-th argument has an illegal value, stopping\n", -info);
      }
      if (info)
      {
        NM_internalData_free(A);
        assert(!NM_internalData(A)->isCholeskyfactorized);
      }
    }
    break;
    case NM_SPARSE_BLOCK: /* sparse block -> triplet -> csc */
    case NM_SPARSE:
    {
      NSM_linear_solver_params* p = NSM_linearSolverParams(A);
      assert(!NM_internalData(A)->isCholeskyfactorized);
      switch (p->solver)
      {
      case NSM_CSPARSE:
        numerics_printf_verbose(2, "NM_Cholesky_factorize, using CSparse (chol_factorization)");

        if (!(p->dWork && p->linear_solver_data))
        {
          assert(!NSM_workspace(p));
          assert(!NSM_linear_solver_data(p));
          assert(!p->solver_free_hook);

          p->solver_free_hook = &NSM_clear_p;
          p->dWork = (double*) malloc(A->size1 * sizeof(double));
          p->dWorkSize = A->size1;
        };

        CSparseMatrix_factors* cs_chol_A = (CSparseMatrix_factors*) malloc(sizeof(CSparseMatrix_factors));

        info = !CSparseMatrix_chol_factorization(1, NM_csc(A),  cs_chol_A);

        if (info)
        {
          numerics_printf_verbose(2, "NM_Cholesky_factorize: csparse factorization failed.");
        }
        assert(!p->linear_solver_data);
        p->linear_solver_data = cs_chol_A;
        break;
#ifdef WITH_MUMPS
      case NSM_MUMPS:
      {
        if(verbose >= 2)
        {
          printf("NM_Cholesky_factorize: using MUMPS\n");
        }
        if(!NM_MUMPS_id(A)->job || (NM_MUMPS_id(A)->job == -2))
        {
          /* the mumps instance is initialized (call with job=-1) */
          NM_MUMPS_set_control_params(A);
          NM_MUMPS_set_sym(A, 1); /*  symmetric positive definite */
          NM_MUMPS(A, -1);
          if ((NM_MUMPS_icntl(A, 1) == -1 ||
               NM_MUMPS_icntl(A, 2) == -1 ||
               NM_MUMPS_icntl(A, 3) == -1 ||
               verbose) ||
              (NM_MUMPS_icntl(A, 1) != -1 ||
               NM_MUMPS_icntl(A, 2) != -1 ||
               NM_MUMPS_icntl(A, 3) != -1 ||
               !verbose))
          {
            NM_MUMPS_set_verbosity(A, verbose);
          }
          /* NM_MUMPS_set_icntl(A, 24, 1); // Null pivot row detection */
          /* NM_MUMPS_set_cntl(A, 5, 1.e20); // Fixation, recommended value */
        }

        NM_MUMPS_set_matrix(A);

        NM_MUMPS(A, 4); /* analyzis,factorization */

        DMUMPS_STRUC_C* mumps_id = NM_MUMPS_id(A);

        info = mumps_id->info[0];

        /* MUMPS can return info codes with negative value */
        if(info)
        {
          if(verbose > 0)
          {
            fprintf(stderr,"NM_Cholesky_factorize: MUMPS fails : info(1)=%d, info(2)=%d\n", info, mumps_id->info[1]);
          }
        }

        /* we should not do that here */
        if(!p->solver_free_hook)
        {
          p->solver_free_hook = &NM_MUMPS_free;
        }
        break;
      }
#endif /* WITH_MUMPS */
      default:
      {
        numerics_printf_verbose(0,"NM_Cholesky_factorize, Unknown solver in NM_SPARSE case." );
        info = 1;
      }
      }
      break;
    }
    default:
      assert (0 && "NM_Cholesky_factors: unknown storageType");
    }

    if (!info)
    {
      NM_internalData(A)->isCholeskyfactorized = true;
    }
    else
    {
      assert(NM_internalData(A)->isCholeskyfactorized == false);
    }
  }

  assert (NM_Cholesky_factorized(Ao) == NM_Cholesky_factorized(A));
  DEBUG_END("int NM_Cholesky_factorize(NumericsMatrix* Ao) \n");
  return info;
}


int NM_Cholesky_solve(NumericsMatrix* Ao, double *b, unsigned int nrhs)
{

  lapack_int info = 1;
  /* factorization is done on destructible part only if
   * !A->internalData->isLUfactorized */
  NM_Cholesky_factorize(Ao);

  /* get the destructible part of the matrix */
  NumericsMatrix *A = Ao->destructible;

  if (NM_Cholesky_factorized(A))
  {

    DEBUG_BEGIN("NM_Cholesky_solve(NumericsMatrix* A, double *b, unsigned int nrhs)\n");
    assert(A->size0 == A->size1);

    switch (A->storageType)
    {
    case NM_DENSE:
    {
      assert(A->matrix0);

      numerics_printf_verbose(2, "NM_Cholesky_solve, using LAPACK (DPOTRS)" );

      /* dpotrf is called in NM_Cholesky_factorize */
      DEBUG_PRINT("Start to call DPOTRS for NM_DENSE storage\n");

      numerics_printf_verbose(2,"NM_Cholesky_solve, we solve with given factors" );
      DEBUG_PRINT("Start to call DPOTRS for NM_DENSE storage\n");
      DEBUG_EXPR(NV_display(b,A->size0));
      DPOTRS(LA_UP, A->size0, nrhs, A->matrix0, A->size0, b, A->size0, &info);
      DEBUG_EXPR(NV_display(b,A->size0));
      DEBUG_PRINT("End of call DPOTRS for NM_DENSE storage\n");

      if (info < 0)
      {
        numerics_printf_verbose(2,"NM_Cholesky_solve: dense Cholesky solve DPOTRS failed. The %d-th argument has an illegal value\n", -info);
      }
      break;
    }

    case NM_SPARSE_BLOCK: /* sparse block -> triplet -> csc */
    case NM_SPARSE:
    {
      NSM_linear_solver_params* p = NSM_linearSolverParams(A);
      switch (p->solver)
      {
      case NSM_CSPARSE:
      {
        numerics_printf_verbose(2,"NM_Cholesky_solve, using CSparse" );

        numerics_printf_verbose(2,"NM_Cholesky_solve, we solve with given factors" );
        for(unsigned int j=0; j < nrhs ; j++ )
        {
          info = !CSparseMatrix_chol_solve((CSparseMatrix_factors *)NSM_linear_solver_data(p), NSM_workspace(p), &b[j*A->size1]);
        }
        break;
      }
#ifdef WITH_MUMPS
      case NSM_MUMPS:
      {
        if(verbose >= 2)
        {
          printf("NM_Cholesky_solve: using MUMPS\n");
        }

        assert (NM_MUMPS_id(A)->job); /* this means that least a
                                       * factorization has already been
                                       * done */

        DMUMPS_STRUC_C* mumps_id = NM_MUMPS_id(A);

        NM_MUMPS_set_dense_rhs(A, nrhs, b);

        NM_MUMPS(A, 3); /* solve */
        info = mumps_id->info[0];

        /* MUMPS can return info codes with negative value */
        if(info)
        {
          if(verbose > 0)
          {
            fprintf(stderr,"NM_Cholesky_solve: MUMPS fails : info(1)=%d, info(2)=%d\n", info, mumps_id->info[1]);
          }
        }
        break;
      }
#endif /* WITH_MUMPS */
      default:
      {
        fprintf(stderr, "NM_Cholesky_solve: unknown sparse linearsolver %d\n", p->solver);
        exit(EXIT_FAILURE);
      }
      break;
      }
      break;
    }
    default:
      assert (0 && "NM_Cholesky_solve unknown storageType");
    }


    /* WARNING: cs returns 0 (false) for failed and 1 (true) for ok
       CHECK_RETURN is ok for cs, but not for MUMPS and others */
    /* some time we cannot find a solution to a linear system, and its fine, for
     * instance with the minFBLSA. Therefore, we should not check here for
     * problems, but the calling function has to check the return code.*/
//  CHECK_RETURN(info);
    DEBUG_END("NM_Cholesky_solve(NumericsMatrix* A, double *b, unsigned keep)\n");
  }

  return info;
}

int NM_Cholesky_solve_matrix_rhs(NumericsMatrix* Ao, NumericsMatrix* B)
{

  lapack_int info = 1;

  /* factorization is done on destructible part only if
   * !A->internalData->isCHOLESKYfactorized */
  NM_Cholesky_factorize(Ao);

  /* get the destructible part of the matrix */
  NumericsMatrix *A = Ao->destructible;

  if (NM_Cholesky_factorized(A))
  {

    DEBUG_BEGIN("NM_Cholesky_solve_matrix_rhs(NumericsMatrix* Ao, NumericsMatrix* B)\n");
    assert(A->size0 == A->size1);

    if (B->storageType == NM_DENSE)
    {
      assert(B->matrix0);
      info = NM_Cholesky_solve(A, B->matrix0, B->size1);
    }
    else if ((B->storageType == NM_SPARSE) || (B->storageType == NM_SPARSE_BLOCK))
    {
      switch (A->storageType)
      {
      case NM_DENSE:
      {
        numerics_error("NM_Cholesky_solve_matrix_rhs", "Solving Linear system with a dense matrix and a sparse matrix rhs is not implemented since it requires to copy the rhs into dense matrix." );
        break;
      }

      case NM_SPARSE_BLOCK: /* sparse block -> triplet -> csc */
      case NM_SPARSE:
      {
        NSM_linear_solver_params* p = NSM_linearSolverParams(A);
        switch (p->solver)
        {
        case NSM_CSPARSE:
        {
          numerics_printf_verbose(2,"NM_Cholesky_solve_matrix_rhs, using CSparse" );
          numerics_printf_verbose(2,"NM_Cholesky_solve_matrix_rhs, we solve with given factors" );

          if (!p->dWork)
          {
            assert(!NSM_workspace(p));
            p->dWork = (double*) malloc(A->size1 * sizeof(double));
            p->dWorkSize = A->size1;
          };


          CSparseMatrix *X = cs_spalloc(NM_csc(B)->m, NM_csc(B)->n, NM_csc(B)->nzmax, 1,0); /* csc format */

          info = !CSparseMatrix_chol_spsolve((CSparseMatrix_factors *)NSM_linear_solver_data(p), X, NM_csc(B));
          //Invalidation
          B->matrix2->origin = NSM_CSC;
          NM_clearCSR(B);
          NM_clearCSCTranspose(B);
          NM_clearTriplet(B);
          NM_clearHalfTriplet(B);

          break;
        }
#ifdef WITH_MUMPS
        case NSM_MUMPS:
        {
          numerics_printf_verbose(2,"NM_Cholesky_solve: using MUMPS\n");

          assert (NM_MUMPS_id(A)->job); /* this means that at least a
                                         * factorization has already been
                                         * done */

          DMUMPS_STRUC_C* mumps_id = NM_MUMPS_id(A);

          NM_MUMPS_set_sparse_rhs(A, B);

          NM_MUMPS(A, 3); /* solve */
          info = mumps_id->info[0];

          /* MUMPS can return info codes with negative value */
          if(info)
          {
            if(verbose > 0)
            {
              fprintf(stderr,"NM_Cholesky_solve: MUMPS fails : info(1)=%d, info(2)=%d\n", info, mumps_id->info[1]);
            }
          }

          /* solution is returned in the DENSE format in the B matrix */
          NM_clear(B);
          B->storageType = NM_DENSE;
          B->size0 = A->size0;
          B->size1 = mumps_id->nrhs;
          B->matrix0 = (double *) malloc(
            A->size0*mumps_id->nrhs*sizeof(double));
          memcpy(B->matrix0, mumps_id->rhs,
                 A->size0*mumps_id->nrhs*sizeof(double));
          break;
        }
#endif /* WITH_MUMPS */
        default:
        {
          fprintf(stderr, "NM_Cholesky_solve: unknown sparse linearsolver %d\n", p->solver);
          exit(EXIT_FAILURE);
        }
        break;
        }
        break;
      }
      default:
        assert (0 && "NM_Cholesky_solve_matrix_rhs unknown storageType");
      }


      /* WARNING: cs returns 0 (false) for failed and 1 (true) for ok
         CHECK_RETURN is ok for cs, but not for MUMPS and others */
      /* some time we cannot find a solution to a linear system, and its fine, for
       * instance with the minFBLSA. Therefore, we should not check here for
       * problems, but the calling function has to check the return code.*/

    }
  }
  DEBUG_END("NM_Cholesky_solve_matrix_rhs(NumericsMatrix* Ao, NumericsMatrix* B)\n");
  return info;
}


int NM_LDLT_factorize(NumericsMatrix* Ao)
{
  DEBUG_BEGIN("int NM_LDLT_factorize(NumericsMatrix* Ao) \n");
  lapack_int info = 0;
  assert(Ao->destructible); /* by default Ao->destructible == Ao */
  NumericsMatrix* A = Ao->destructible;

  if (!NM_LDLT_factorized(Ao))
  {

#ifdef FACTORIZATION_DEBUG
    if (NM_internalData(Ao)->values_sha1_count > 0)
    {
      if(NM_check_values_sha1(Ao))
      {
        numerics_error("NM_LDLT_factorize", "this matrix is already factorized");
      }
    }

    NM_set_values_sha1(Ao);
#endif

    switch (A->storageType)
    {
    case NM_DENSE:
    {
      assert(A->matrix0);

      numerics_printf_verbose(2,"NM_LDLT_factorize, using LAPACK (SYTRF)" );
      lapack_int* ipiv = (lapack_int*)NM_iWork(A, A->size0, sizeof(lapack_int));
      DEBUG_PRINTF("iwork are initialized with size %i and %i\n",A->size0*A->size1,A->size0 );

      numerics_printf_verbose(2,"NM_LDLT_factorize, we compute factors and keep them in place" );
      DEBUG_PRINT("Start to call DSYTRF for NM_DENSE storage\n");
      DSYTRF(LA_UP, A->size1, A->matrix0, A->size0, ipiv, &info);
      DEBUG_PRINT("end of call DSYTRF for NM_DENSE storage\n");

      if (info > 0)
      {
        fprintf(stderr,"NM_LDLT_factorize: LDLT factorisation DSYTRF failed. The %d-th diagonal element is 0\n", info);
      }
      else if (info < 0)
      {
        fprintf(stderr, "NM_LDLT_factorize: LDLT factorisation DSYTRF failed. The %d-th argument has an illegal value, stopping\n", -info);
      }
      if (info)
      {
        NM_internalData_free(A);
        assert(!NM_internalData(A)->isLDLTfactorized);
      }
    }
    break;
    case NM_SPARSE_BLOCK: /* sparse block -> triplet -> csc */
    case NM_SPARSE:
    {
      NSM_linear_solver_params* p = NSM_linearSolverParams(A);
      assert(!NM_internalData(A)->isLDLTfactorized);
      switch (p->LDLT_solver)
      {
      case NSM_CSPARSE:
      {        numerics_printf_verbose(2, "NM_LDLT_factorize, using SuiteSparse (LDL )");

        if (!(p->dWork && p->linear_solver_data))
        {
          assert(!NSM_workspace(p));
          assert(!NSM_linear_solver_data(p));
          assert(!p->solver_free_hook);

          p->solver_free_hook = &NSM_clear_p;
          p->dWork = (double*) malloc(A->size1 * sizeof(double));
          p->dWorkSize = A->size1;
        };

        CSparseMatrix_factors* cs_ldlt_A = (CSparseMatrix_factors*) malloc(sizeof(CSparseMatrix_factors));

        info = !CSparseMatrix_ldlt_factorization(1, NM_csc(A),  cs_ldlt_A);

        if (info)
        {
          numerics_printf_verbose(2, "NM_LDLT_factorize: LDL (SuiteSparse) factorization failed.");
        }
        assert(!p->linear_solver_data);
        p->linear_solver_data = cs_ldlt_A;
        break;
      }
#ifdef WITH_MA57
      case NSM_HSL:
      {
        int  n = A->size0;
        CSparseMatrix * A_halftriplet= NM_half_triplet(A);
        CS_INT nz = A_halftriplet->nz;

        FILE * logfile = fopen("lbl.log", "w");
        int lblsolver = 1;  // MA27=0, MA57=1.

        // ... Initialize LBL data structure.
        LBL_Data * lbl = LBL_Initialize(nz, n, logfile, lblsolver);

        for(int i = 0; i < nz; i++)
        {
          lbl->irn[i]=A_halftriplet->i[i]+1;
          lbl->jcn[i]=A_halftriplet->p[i]+1;
        }

        p->linear_solver_data = lbl;

        if(!p->solver_free_hook)
        {
          p->solver_free_hook = &NM_MA57_free;
        }

        // Optionally, we may set some specific parameters (MA57 only).
        //LBL_set_int_parm(lbl, LBL_I_SCALING, 0);        // No scaling.
        //LBL_set_int_parm(lbl, LBL_I_PIV_NUMERICAL, 1);  // Do pivoting.
        //LBL_set_real_parm(lbl, LBL_D_PIV_THRESH, 0.5);  // Pivot threshold.

        // Analyze.
        info  = LBL_Analyze(lbl, 0);  // iflag=0: automatic pivot choice.
        if(info) {
            fprintf(stderr, "NM_LDLT_factorize: LBL_Analyze Error return from Analyze: %d\n", info);
        }

        // Factorize.
        info = LBL_Factorize(lbl, A_halftriplet->x);
        if(info) {
          fprintf(stderr, "NM_LDLT_factorize: LBL_Factorize. Error return from Factorize: %d\n", info);
        }
        // Close logging stream.
        fclose(logfile);

      break;
      }
#endif
#ifdef WITH_MUMPS
      case NSM_MUMPS:
      {
        if(verbose >= 2)
        {
          printf("NM_LDLT_factorize: using MUMPS\n");
        }
        if(!NM_MUMPS_id(A)->job || (NM_MUMPS_id(A)->job == -2))
        {
          /* the mumps instance is initialized (call with job=-1) */
          NM_MUMPS_set_control_params(A);
          NM_MUMPS_set_sym(A, 2); /*  general symmetric (LDLT)  */
          NM_MUMPS(A, -1);
          if ((NM_MUMPS_icntl(A, 1) == -1 ||
               NM_MUMPS_icntl(A, 2) == -1 ||
               NM_MUMPS_icntl(A, 3) == -1 ||
               verbose) ||
              (NM_MUMPS_icntl(A, 1) != -1 ||
               NM_MUMPS_icntl(A, 2) != -1 ||
               NM_MUMPS_icntl(A, 3) != -1 ||
               !verbose))
          {
            NM_MUMPS_set_verbosity(A, verbose);
          }
          /* NM_MUMPS_set_icntl(A, 24, 1); // Null pivot row detection */
          /* NM_MUMPS_set_cntl(A, 5, 1.e20); // Fixation, recommended value */
        }

        NM_MUMPS_set_matrix(A);

        NM_MUMPS(A, 4); /* analyzis,factorization */

        DMUMPS_STRUC_C* mumps_id = NM_MUMPS_id(A);

        info = mumps_id->info[0];

        /* MUMPS can return info codes with negative value */
        if(info)
        {
          if(verbose > 0)
          {
            fprintf(stderr,"NM_LDLT_factorize: MUMPS fails : info(1)=%d, info(2)=%d\n", info, mumps_id->info[1]);
          }
        }

        /* we should not do that here */
        if(!p->solver_free_hook)
        {
          p->solver_free_hook = &NM_MUMPS_free;
        }
        break;
      }
#endif /* WITH_MUMPS */
      default:
      {
        numerics_printf_verbose(0,"NM_LDLT_factorize, Unknown solver in NM_SPARSE case." );
        info = 1;
      }
      }
      break;
    }
    default:
      assert (0 && "NM_LDLT_factors: unknown storageType");
    }

    if (!info)
    {
      NM_internalData(A)->isLDLTfactorized = true;
    }
    else
    {
      assert(NM_internalData(A)->isLDLTfactorized == false);
    }
  }

  assert (NM_LDLT_factorized(Ao) == NM_LDLT_factorized(A));
  DEBUG_END("int NM_LDLT_factorize(NumericsMatrix* Ao) \n");
  return info;
}


int NM_LDLT_solve(NumericsMatrix* Ao, double *b, unsigned int nrhs)
{

  lapack_int info = 1;
  /* factorization is done on destructible part only if
   * !A->internalData->isLUfactorized */
  NM_LDLT_factorize(Ao);

  /* get the destructible part of the matrix */
  NumericsMatrix *A = Ao->destructible;

  if (NM_LDLT_factorized(A))
  {

    DEBUG_BEGIN("NM_LDLT_solve(NumericsMatrix* A, double *b, unsigned int nrhs)\n");
    assert(A->size0 == A->size1);

    switch (A->storageType)
    {
    case NM_DENSE:
    {
      assert(A->matrix0);

      numerics_printf_verbose(2, "NM_LDLT_solve, using LAPACK (DSYTRS)" );

      /* dpotrf is called in NM_LDLT_factorize */
      DEBUG_PRINT("Start to call DSYTRS for NM_DENSE storage\n");

      numerics_printf_verbose(2,"NM_LDLT_solve, we solve with given factors" );
      lapack_int* ipiv = (lapack_int*)NM_iWork(A, A->size0, sizeof(lapack_int));

      DEBUG_PRINT("Start to call DSYTRS for NM_DENSE storage\n");
      DEBUG_EXPR(NV_display(b,A->size0));
      DSYTRS(LA_UP, A->size0, nrhs, A->matrix0, A->size0,  ipiv, b, A->size0, &info);
      DEBUG_EXPR(NV_display(b,A->size0));
      DEBUG_PRINT("End of call DSYTRS for NM_DENSE storage\n");

      if (info < 0)
      {
        numerics_printf_verbose(2,"NM_LDLT_solve: dense LDLT solve DPOTRS failed. The %d-th argument has an illegal value\n", -info);
      }
      break;
    }

    case NM_SPARSE_BLOCK: /* sparse block -> triplet -> csc */
    case NM_SPARSE:
    {
      NSM_linear_solver_params* p = NSM_linearSolverParams(A);
      switch (p->LDLT_solver)
      {
      case NSM_CSPARSE:
      {
        numerics_printf_verbose(2,"NM_LDLT_solve, using SuiteSparse" );

        numerics_printf_verbose(2,"NM_LDLT_solve, we solve with given factors" );
        for(unsigned int j=0; j < nrhs ; j++ )
        {
          info = !CSparseMatrix_ldlt_solve((CSparseMatrix_factors *)NSM_linear_solver_data(p), NSM_workspace(p), &b[j*A->size1]);
        }
        break;
      }
#ifdef WITH_MA57
      case NSM_HSL:
      {
        LBL_Data * lbl = (LBL_Data *)p->linear_solver_data;
        // Solve.
        for (int irhs=0; irhs <nrhs ; irhs++)
        {
          info = LBL_Solve(lbl, &b[irhs*A->size1]); // MA57 is able to accept multiple rhs but the C wrapper lbl not.
          if(info)
          {
            fprintf(stderr, "NM_LDLT_solve. LBL_Solve error return from Solve: %d\n", info);
          }
        }
        break;
      }
#endif
#ifdef WITH_MUMPS
      case NSM_MUMPS:
      {
        if(verbose >= 2)
        {
          printf("NM_LDLT_solve: using MUMPS\n");
        }

        assert (NM_MUMPS_id(A)->job); /* this means that least a
                                       * factorization has already been
                                       * done */

        DMUMPS_STRUC_C* mumps_id = NM_MUMPS_id(A);

        NM_MUMPS_set_dense_rhs(A, nrhs, b);

        NM_MUMPS(A, 3); /* solve */
        info = mumps_id->info[0];

        /* MUMPS can return info codes with negative value */
        if(info)
        {
          if(verbose > 0)
          {
            fprintf(stderr,"NM_LDLT_solve: MUMPS fails : info(1)=%d, info(2)=%d\n", info, mumps_id->info[1]);
          }
        }
        break;
      }
#endif /* WITH_MUMPS */
      default:
      {
        fprintf(stderr, "NM_LDLT_solve: unknown sparse linearsolver %d\n", p->LDLT_solver);
        exit(EXIT_FAILURE);
      }
      break;
      }
      break;
    }
    default:
      assert (0 && "NM_LDLT_solve unknown storageType");
    }


    /* WARNING: cs returns 0 (false) for failed and 1 (true) for ok
       CHECK_RETURN is ok for cs, but not for MUMPS and others */
    /* some time we cannot find a solution to a linear system, and its fine, for
     * instance with the minFBLSA. Therefore, we should not check here for
     * problems, but the calling function has to check the return code.*/
//  CHECK_RETURN(info);
    DEBUG_END("NM_LDLT_solve(NumericsMatrix* A, double *b, unsigned keep)\n");
  }

  return info;
}

int NM_LDLT_refine(NumericsMatrix* Ao, double *x , double *b, unsigned int nrhs, double tol, int maxitref, int job )
{

  lapack_int info = 1;
  /* factorization is done on destructible part only if
   * !A->internalData->isLUfactorized */
  NM_LDLT_factorize(Ao);

  /* get the destructible part of the matrix */
  NumericsMatrix *A = Ao->destructible;

  if (NM_LDLT_factorized(A))
  {

    DEBUG_BEGIN("NM_LDLT_refine(NumericsMatrix* A, double *b, unsigned int nrhs)\n");
    assert(A->size0 == A->size1);

    switch (A->storageType)
    {
    case NM_DENSE:  
    case NM_SPARSE_BLOCK: /* sparse block -> triplet -> csc */
    case NM_SPARSE:
    {
      NSM_linear_solver_params* p = NSM_linearSolverParams(A);
      switch (p->LDLT_solver)
      {
#ifdef WITH_MA57
      case NSM_HSL:
      {
        LBL_Data * lbl = (LBL_Data *)p->linear_solver_data;
        // Solve.
        for (int irhs=0; irhs <nrhs ; irhs++)
        {
          info = LBL_Refine(lbl, &x[irhs*A->size1], &b[irhs*A->size1], NM_half_triplet(A)->x,
                            tol, maxitref, job); // MA57 is able to accept multiple rhs but the C wrapper lbl not.
          if(info)
          {
            fprintf(stderr, "NM_LDLT_refine. LBL_Refine error return from Refine: %d\n", info);
          }
        }
        break;
      }
#endif
      default:
      {
        fprintf(stderr, "NM_LDLT_refine: unknown sparse linearrefiner %d\n", p->LDLT_solver);
        exit(EXIT_FAILURE);
      }
      break;
      }
      break;
    }
    default:
      assert (0 && "NM_LDLT_refine unknown storageType");
    }


    /* WARNING: cs returns 0 (false) for failed and 1 (true) for ok
       CHECK_RETURN is ok for cs, but not for MUMPS and others */
    /* some time we cannot find a solution to a linear system, and its fine, for
     * instance with the minFBLSA. Therefore, we should not check here for
     * problems, but the calling function has to check the return code.*/
//  CHECK_RETURN(info);
    DEBUG_END("NM_LDLT_refine(NumericsMatrix* A, double *b, unsigned keep)\n");
  }

  return info;
}
