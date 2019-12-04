/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#include "CSparseMatrix.h"
#include "NumericsMatrix_internal.h"
#include "NumericsSparseMatrix.h"
#include "SiconosCompat.h"
#include "SparseBlockMatrix.h"
#include "NM_MUMPS.h"
#include "NM_MPI.h"
#include "NM_conversions.h"
#include "SiconosLapack.h"
#include "numerics_verbose.h"
#include "sanitizer.h"
/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"

#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"
#endif

#ifdef WITH_MKL_SPBLAS
#include "MKL_common.h"
#include "NM_MKL_spblas.h"
#endif


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

  int storage = A->storageType;

  /* double* storage */
  switch (storage)
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

  int storage = A->storageType;

  /* double* storage */
  if (storage == 0)
  {
    int incx = A->size0, incy = 1;
    double* mat = A->matrix0;
    if (init == 0) /* y += subAx */
    {
      for (int row = 0; row < sizeY; row++)
        y[row] += cblas_ddot(sizeX, &mat[currentRowNumber + row], incx, x, incy);
    }
    else
    {
      for (int row = 0; row < sizeY; row++)
        y[row] = cblas_ddot(sizeX, &mat[currentRowNumber + row], incx, x, incy);
    }

  }
  /* SparseBlock storage */
  else if (storage == 1)
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
    if (xsave) xSave = xsave;
    else xSave = (double*) malloc(sizeY * sizeof(double));

    memcpy(xSave, &x[row_start], sizeY * sizeof(double));
    memset(&x[row_start], 0, sizeY * sizeof(double));

    double * M = A->matrix0 + row_start;
    int incx = A->size0;
    int incy = 1;
    if (init)
    {
      memset(y, 0, sizeY*sizeof(double));
    }
    for (size_t i = 0; i < sizeY; ++i)
    {
      y[i] += cblas_ddot(A->size0, M, incx, x, incy);
      ++M;
    }

    memcpy(&x[row_start], xSave, sizeY*sizeof(double));

    if (!xsave)
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
    if (xsave) xSave = xsave;
    else xSave = (double*) malloc(sizeY * sizeof(double));

    memcpy(xSave, &x[row_start], sizeY * sizeof(double));
    memset(&x[row_start], 0, sizeY * sizeof(double));

    if (init)
    {
      memset(y, 0, sizeY*sizeof(double));
    }

    CSparseMatrix* M;
    if (A->matrix2->origin == NSM_CSR)
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

    for (size_t i = 0, j = row_start; i < sizeY; ++i, ++j)
    {
      for (CS_INT p = Mp[j]; p < Mp[j+1]; ++p)
      {
        y[i] += Mx[p] * x[Mi[p]];
      }
    }

    memcpy(&x[row_start], xSave, sizeY*sizeof(double));

    if (!xsave)
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
    if (init)
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
    if (init)
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
    if (A->matrix2->origin == NSM_CSR)
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

    for (size_t i = 0, j = row_start; i < 3; ++i, ++j)
    {
      for (CS_INT p = Mp[j]; p < Mp[j+1]; ++p)
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
    if (init)
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
    if (init)
    {
      y[0] = 0.;
    }

    double rin = x[row_start] ;
    x[row_start] = 0.;

    CSparseMatrix* M;
    if (A->matrix2->origin == NSM_CSR)
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

    for (CS_INT p = Mp[row_start]; p < Mp[row_start+1]; ++p)
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
  if (m->internalData)
  {
    if (m->internalData->iWork)
    {
      assert (m->internalData->iWorkSize > 0);
      free(m->internalData->iWork);
      m->internalData->iWork = NULL;
    }
    if (m->internalData->dWork)
    {
      assert (m->internalData->dWorkSize > 0);
      free(m->internalData->dWork);
      m->internalData->dWork = NULL;
    }
    free(m->internalData);
    m->internalData = NULL;
  }
}
void NM_internalData_copy(const NumericsMatrix* const A, NumericsMatrix* B )
  {

    if (!A->internalData)/* no internal data in A */
    {
      if (B->internalData)
      {
        NM_internalData_free(B);
      }
    }
    else
    {
      if (!B->internalData)
      {
        NM_internalData_new(B);
      }


      size_t sizeof_elt = max(A->internalData->sizeof_elt,sizeof(int));

      if (A->internalData->iWork)
      {

        int size = A->internalData->iWorkSize / sizeof_elt;


        if (! B->internalData->iWork)
        {
          B->internalData->iWork = malloc(size*sizeof_elt);
          B->internalData->iWorkSize = A->internalData->iWorkSize;
        }
        else
        {
          if( A->internalData->iWorkSize != B->internalData->iWorkSize)
          {
            B->internalData->iWorkSize = A->internalData->iWorkSize;
            B->internalData->iWork = realloc(B->internalData->iWork,size*sizeof_elt);
          }
        }

        memcpy(B->internalData->iWork, A->internalData->iWork, A->internalData->iWorkSize);
      }
      else
      {
        if (B->internalData->iWork)
          free(B->internalData->iWork);
        B->internalData->iWorkSize=0;
      }



      if (A->internalData->dWork)
      {
        if (! B->internalData->dWork)
        {
          B->internalData->dWork = malloc(A->internalData->dWorkSize*sizeof(double));
          B->internalData->dWorkSize = A->internalData->dWorkSize;
        }
        else
        {
          if( A->internalData->dWorkSize != B->internalData->dWorkSize)
          {
            B->internalData->dWorkSize = A->internalData->dWorkSize;
            B->internalData->dWork = realloc(B->internalData->dWork,A->internalData->dWorkSize*sizeof(double));
          }
        }

        for (unsigned i =0; i < A->internalData->dWorkSize; i++)
          B->internalData->dWork[i] = A->internalData->dWork[i];
      }
      else
      {
        if (B->internalData->dWork)
          free(B->internalData->dWork);
        B->internalData->dWorkSize=0;
      }
      B->internalData->isLUfactorized = A->internalData->isLUfactorized;
      B->internalData->isInversed = A->internalData->isInversed;
    }

  }
void NM_free(NumericsMatrix* m)
{
  assert(m && "NM_free, m == NULL");

  NM_clearDense(m);
  NM_clearSparseBlock(m);
  NM_clearSparse(m);

  NM_internalData_free(m);
}

void NM_dense_display_matlab(double * m, int nRow, int nCol, int lDim)
{
  if (m)
  {
    int lin, col;
    if (lDim == 0)
      lDim = nRow;
    printf("Matrix of size\t%d\t x \t%d =\n[", nRow, nCol);
    if (nRow == 0)
    {
      printf("]\n");
    }
    if (nCol == 0)
    {
      printf("]\n");
    }

    for (lin = 0; lin < nRow; lin++)
    {
      for (col = 0; col < nCol; col++)
      {
        printf(" %.15e", m[lin + col * lDim]);
        if (col != nCol - 1)
          printf(",");
      }
      if (lin != nRow - 1)
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
  if (m)
  {
    int lin, col;
    if (lDim == 0)
      lDim = nRow;
    printf("Matrix of size\t%d\t x \t%d =\n[", nRow, nCol);
    if (nRow == 0)
    {
      printf("]\n");
    }
    if (nCol == 0)
    {
      printf("]\n");
    }

    for (lin = 0; lin < nRow; lin++)
    {
      printf("[");
      for (col = 0; col < nCol; col++)
      {
        printf(" %.15e", m[lin + col * lDim]);
        if (col != nCol - 1)
          printf(",");
      }
      if (lin != nRow - 1)
        printf("],\n");
      else
        printf("]\t ]\n");
    }
  }
  else
    printf("Matrix : NULL\n");

}

void NM_zentry(NumericsMatrix* M, int i, int j, double val)
{
  switch (M->storageType)
  {
  case NM_DENSE:
  {
    // column major
    M->matrix0[i+j*M->size0] = val;
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    CHECK_RETURN(SBM_zentry(M->matrix1, i, j, val));
    break;
  }
  case NM_SPARSE:
  {
    assert(M->matrix2);
    switch (M->matrix2->origin)
    {
    case NSM_TRIPLET:
    {
      assert(M->matrix2->triplet);
      CHECK_RETURN(CSparseMatrix_zentry(M->matrix2->triplet, i, j, val));
      break;
    }
    case NSM_HALF_TRIPLET:
    {
      assert(M->matrix2->half_triplet);
      CHECK_RETURN(CSparseMatrix_symmetric_zentry(M->matrix2->triplet, i, j, val));
      break;
    }
    case NSM_CSC:
    {
      assert(M->matrix2->csc);
      CHECK_RETURN(CSparseMatrix_zentry(NM_triplet(M), i, j, val));
      M->matrix2->origin= NSM_TRIPLET;
      NM_clearCSC(M);
      //NM_display(M);
      //NM_csc_alloc(M,M->matrix2->triplet->nz);
      //M->matrix2->origin= NSM_CSC;
      NM_csc(M);

      M->matrix2->origin= NSM_CSC;
      NM_clearTriplet(M);
      NM_clearHalfTriplet(M);

      //NM_display(M);
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

  if (i>M->size0-1)
  {
    M->size0 = i+1;
  }
  if (j>M->size1-1)
  {
    M->size1 = j+1;
  }
}


double NM_get_value(NumericsMatrix* M, int i, int j)
{
  assert(M);

  if ((i + 1 > M->size0) || (j + 1 > M->size1) )
  {
    fprintf(stderr, "NM_get_value :: out of range \n");
    exit(EXIT_FAILURE);
  }
  switch (M->storageType)
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
    switch (M->matrix2->origin)
    {
    case NSM_TRIPLET:
    {
      assert(M->matrix2->triplet);
      assert(i <M->matrix2->triplet->m );
      assert(j <M->matrix2->triplet->n );
      CS_INT * Mi =   M->matrix2->triplet->i;
      CS_INT * Mp =   M->matrix2->triplet->p;
      double * Mx =   M->matrix2->triplet->x;

      for (int idx = 0 ; idx < M->matrix2->triplet->nz  ; idx++ )
      {
        if (Mi[idx] ==i && Mp[idx] == j)
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
      assert(i <M->matrix2->triplet->m );
      assert(j <M->matrix2->triplet->n );
      CS_INT * Mi =   M->matrix2->triplet->i;
      CS_INT * Mp =   M->matrix2->triplet->p;
      double * Mx =   M->matrix2->triplet->x;

      for (int idx = 0 ; idx < M->matrix2->triplet->nz  ; idx++ )
      {
        if ((Mi[idx] ==i && Mp[idx] == j)||(Mi[idx] ==j && Mp[idx] == i))
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

      for (int row = Mp[j]; row < Mp[j+1] ; row++)
      {
        if (i == Mi[row])
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
  return NM_compare(A, B, DBL_EPSILON);
};

bool NM_compare(NumericsMatrix* A, NumericsMatrix* B, double tol)
{
  DEBUG_BEGIN("NM_compare(NumericsMatrix* A, NumericsMatrix* B, double tol)\n");
  assert(A);
  assert(B);
  if (A->size0 != B->size0)
  {
    return false;
  }
  if (A->size1 != B->size1)
  {
    return false;
  }
  for (int i =0; i< A->size0 ; i++)
  {
    for (int j =0; j< A->size1 ; j++)
    {
      /* DEBUG_PRINTF("error %i %i = %e\n",i,j, fabs(NM_get_value(A, i, j) - NM_get_value(B, i, j))); */
      if (fabs(NM_get_value(A, i, j) - NM_get_value(B, i, j)) >= tol)
      {

        DEBUG_PRINTF("A(%i,%i) = %e\t, B(%i,%i) = %e\t,  error = %e\n",
                     i,j, NM_get_value(A, i, j),
                     i,j, NM_get_value(B, i, j),
                     fabs(NM_get_value(A, i, j) - NM_get_value(B, i, j)));
        return false;
      }
    }
  }
  DEBUG_END("NM_compare(NumericsMatrix* A, NumericsMatrix* B, double tol)\n");
  return true;
};


void NM_vector_display(double * m, int nRow)
{
  int lin;
  printf("vector of size\t%d\t =\n[", nRow);
  if (nRow == 0)
  {
    printf("]\n");
  }
  for (lin = 0; lin < nRow; lin++)
  {
    printf(" %.15e", m[lin]);
    if (lin != nRow - 1)
      printf(", ");
    else
      printf("]\n");
  }

}

void NM_display(const NumericsMatrix* const m)
{
  if (! m)
  {
    fprintf(stderr, "Numerics, NumericsMatrix display failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  printf("========== Numerics Matrix\n");

  printf("========== size0 = %i, size1 = %i\n", m->size0, m->size1);

  switch (m->storageType)
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
    switch (m->matrix2->origin)
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
    if (m->matrix2->triplet)
    {
      printf("========== a matrix in format triplet is stored\n" );
      cs_print(m->matrix2->triplet, 0);
    }
    if (m->matrix2->csc)
    {
      printf("========== a matrix in format csc is stored\n" );
      cs_print(m->matrix2->csc, 0);
    }
    if (m->matrix2->trans_csc)
    {
      printf("========== a matrix in format trans_csc is stored\n" );
      cs_print(m->matrix2->trans_csc, 0);
    }
    /* else */
    /* { */
    /*   fprintf(stderr, "display for sparse matrix: no matrix found!\n"); */
    /* } */
    if (m->matrix2->diag_indx)
    {
      printf("========== m->matrix2->diag_indx = %p\n", m->matrix2->diag_indx );
      for (int  i = 0; i < m->size0; ++i) printf("diag_indices[%i] = %li\t ", i, m->matrix2->diag_indx[i]);
    }
    else
    {
      printf("========== m->matrix2->diag_indx --> NULL\n");
    }
    if (m->matrix2->linearSolverParams)
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

  if (m->internalData)
  {
    printf("========== internalData->iWorkSize = %lu\n", m->internalData->iWorkSize );
    printf("========== internalData->iWork = %p\n", m->internalData->iWork );

    printf("========== internalData->dWorkSize = %lu\n", m->internalData->iWorkSize );
    printf("========== internalData->dWork = %p\n", m->internalData->iWork );
  }
  else
  {
    printf("========== internalData = NULL\n" );
  }


}
void NM_display_row_by_row(const NumericsMatrix* const m)
{
  if (! m)
  {
    fprintf(stderr, "Numerics, NumericsMatrix display failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  int storageType = m->storageType;
  if (storageType == NM_DENSE)
  {
    printf("\n ========== Numerics Matrix of dim %dX%d\n", m->size0, m->size1);
    for (int lin = 0; lin < m->size0; lin++)
    {
      for (int col = 0; col < m->size1; col++)
        printf("%lf ", m->matrix0[lin + col * m->size1]);
      printf("\n");
    }
  }
  else if (storageType == NM_SPARSE_BLOCK)
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

  if (! m)
  {
    fprintf(stderr, "Numerics, NM_write_in_file failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }

  fprintf(file, "%d\n", m->storageType);
  fprintf(file, "%d\n", m->size0);
  fprintf(file, "%d\n", m->size1);
  DEBUG_PRINTF("\n ========== storageType = %i\n", m->storageType);

  switch (m->storageType)
  {
  case NM_DENSE:
  {
    fprintf(file, "%i\t%i\n", m->size0, m->size1);
    for (int i = 0; i < m->size1 * m->size0; i++)
    {
      fprintf(file, "%32.24e ", m->matrix0[i]);
      if ((i + 1) % m->size1 == 0)
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
  if (! m)
  {
    fprintf(stderr, "Numerics, NumericsMatrix_write_in_file_python  failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  fprintf(file, "storageType = %d ; \n", m->storageType);
  fprintf(file, "size0 = %d; \n", m->size0);
  fprintf(file, "size1 = %d; \n", m->size1);
  fprintf(file, "data= [");
  for (int i = 0; i < m->size0; i++)
  {
    fprintf(file, "[");
    for (int j = 0; j < m->size1; j++)
    {
      fprintf(file, "%32.24e,\t ", NM_get_value((NumericsMatrix*) m,i,j));
    }
    fprintf(file, "],\n");
  }
  fprintf(file, "]");
}

void NM_write_in_file_scilab(const NumericsMatrix* const m, FILE* file)
{
  if (! m)
  {
    fprintf(stderr, "Numerics, NumericsMatrix printInFile failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  fprintf(file, "storageType = %d ; \n", m->storageType);
  fprintf(file, "size0 = %d; \n", m->size0);
  fprintf(file, "size1 = %d; \n", m->size1);
  fprintf(file, "data= [");
  for (int i = 0; i < m->size0; i++)
  {
    fprintf(file, "[");
    for (int j = 0; j < m->size1; j++)
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
  if (finput == NULL)
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
  if (!m)
  {
    fprintf(stderr, "Numerics, NumericsMatrix NM_read_in_file failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  CHECK_IO(fscanf(file, "%d", &(m->storageType)));
  CHECK_IO(fscanf(file, "%d", &(m->size0)));
  CHECK_IO(fscanf(file, "%d", &(m->size1)));
  int storageType = m->storageType;

  if (storageType == NM_DENSE)
  {
    CHECK_IO(fscanf(file, "%d\t%d\n", &(m->size0), &(m->size1)));

    for (int i = 0; i < m->size1 * m->size0; i++)
    {
      CHECK_IO(fscanf(file, "%le ", &(m->matrix0[i])));
    }
  }
  else if (storageType == NM_SPARSE_BLOCK)
  {
    SBM_read_in_file(m->matrix1, file);
  }
  else if (storageType == NM_SPARSE)
  {
    NumericsMatrix * tmp = NM_new_from_file(file);
    NM_copy(tmp,m);
    NM_free(tmp);
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

  int storageType;
  size_t size0;
  size_t size1;
  int info = 0;
  void* data = NULL;

  CHECK_IO(fscanf(file, "%d", &storageType), &info);
  CHECK_IO(fscanf(file, SN_SIZE_T_F, &size0), &info);
  CHECK_IO(fscanf(file, SN_SIZE_T_F, &size1), &info);

  if (storageType == NM_DENSE)
  {
    CHECK_IO(fscanf(file, SN_SIZE_T_F "\t" SN_SIZE_T_F "\n", &size0, &size1), &info);

    data =  malloc(size1 * size0 * sizeof(double));
    double* data_d = (double *) data;

    for (size_t i = 0; i < size1 * size0; ++i)
    {
      CHECK_IO(fscanf(file, "%le ", &(data_d[i])), &info);
    }
  }
  else if (storageType == NM_SPARSE_BLOCK)
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

NumericsMatrix* NM_new_from_filename(char * filename)
{
  FILE* finput = fopen(filename, "r");
  if (finput == NULL)
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

NumericsMatrix* NM_create_from_filename(char * filename)
{
  FILE* finput = fopen(filename, "r");
  if (finput == NULL)
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
  int storageType = M->storageType;
  switch (storageType)
  {
  case NM_DENSE:
  {
    double* Mptr = M->matrix0 + (M->size0 + 1)*start_row;
    double* Bmat = *Block;
    /* The part of MM which corresponds to the current block is copied into MLocal */
    for (size_t i = 0; i < (size_t) size; ++i)
    {
      memcpy(Bmat, Mptr, size*sizeof(double));
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
  int storageType = M->storageType;
  switch (storageType)
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

void NM_extract_diag_block5(NumericsMatrix* M, int block_row_nb, double ** Block)
{
  int storageType = M->storageType;
  switch (storageType)
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
    int diagPos = SBM_diagonal_block_index(M->matrix1, block_row_nb);
    (*Block) = M->matrix1->block[diagPos];
    break;
  }
  case NM_SPARSE:
  {
    size_t start_row = (size_t)block_row_nb + block_row_nb + block_row_nb;
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
  int storageType = M->storageType;
  switch (storageType)
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
  switch (M->storageType)
  {
  case NM_DENSE:
  {
    for (size_t indx = 0; indx < n*n; indx += n+1) M->matrix0[indx] += alpha;
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    for (size_t ic = 0; ic < n/3; ++ic)
    {
      int diagPos = SBM_diagonal_block_index(M->matrix1, ic);
      M->matrix1->block[diagPos][0] += alpha;
      M->matrix1->block[diagPos][4] += alpha;
      M->matrix1->block[diagPos][8] += alpha;
    }
    break;
  }
  case NM_SPARSE:
  {
    CS_INT* diag_indices = NSM_diag_indices(M);

    DEBUG_EXPR(
      printf("diag_indices:\n");
      for (size_t i = 0; i < n; ++i) printf("diag_indices[%zu] = %li\t ", i, diag_indices[i]);
      );

    double* Mx = NSM_data(M->matrix2);
    for (size_t i = 0; i < n; ++i) Mx[diag_indices[i]] += alpha;

    break;
  }
  default:
    printf("NM_add_to_diag3 :: unsupported matrix storage %d", M->storageType);
    exit(EXIT_FAILURE);
  }
}

NumericsMatrix *  NM_add(double alpha, NumericsMatrix* A, double beta, NumericsMatrix* B)
{
  assert(A->size0 == B->size0 && "NM_add :: A->size0 != A->size0 ");
  assert(A->size1 == B->size1 && "NM_add :: A->size1 != A->size1 ");


  /* The storageType  for C inherits from A except for NM_SPARSE_BLOCK */
  NumericsMatrix *C = NM_create(A->storageType, A->size0, A->size1);

  /* should we copy the whole internal data ? */
  /*NM_internalData_copy(A, C);*/
  NM_MPI_copy(A, C);
  NM_MUMPS_copy(A, C);

  switch (A->storageType)
  {
  case NM_DENSE:
  {
    int nm= A->size0*A->size1;
    cblas_dcopy(nm, A->matrix0, 1, C->matrix0, 1);
    cblas_dscal(nm, alpha, C->matrix0,1);
    switch (B->storageType)
    {
    case NM_DENSE:
    {
      cblas_daxpy(nm, beta, B->matrix0, 1, C->matrix0,1 );
      break;
    }
    case NM_SPARSE_BLOCK:
    case NM_SPARSE:
    {
      NumericsMatrix* B_dense = NM_create(NM_DENSE, A->size0, A->size1);
      NM_to_dense(B, B_dense);
      cblas_daxpy(nm, beta, B_dense->matrix0, 1, C->matrix0,1 );
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

  switch (A->storageType)
  {
  case NM_DENSE:
  {
    int nm= A->size0*A->size1;
    cblas_dscal(nm, alpha, A->matrix0,1);
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    SBM_scal(alpha, A->matrix1);
    break;
  }
  case NM_SPARSE:
  {
    CSparseMatrix_scal(alpha, NM_csc(A));
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
  switch (mat->storageType)
  {
    case NM_DENSE:
      data = malloc(size0*size1*sizeof(double));
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
  NumericsMatrix* M = (NumericsMatrix*) malloc(sizeof(NumericsMatrix));
  M->storageType = -1;
  M->size0 = 0;
  M->size1 = 0;
  NM_null(M);

  return M;
}

NumericsMatrix* NM_eye(int size)
{
  NumericsMatrix* M = NM_create(NM_SPARSE, size, size);
  M->matrix2 = NSM_triplet_eye(size);
  return M;
}
NumericsMatrix* NM_create(int storageType, int size0, int size1)
{
  NumericsMatrix* M = NM_new();

  void* data;

  switch (storageType)
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


void NM_fill(NumericsMatrix* M, int storageType, int size0, int size1, void* data)
{

  assert(M);
  M->storageType = storageType;
  M->size0 = size0;
  M->size1 = size1;

  NM_null(M);

  if (data)
  {
    switch (storageType)
    {
      case NM_DENSE:
        M->matrix0 = (double*) data;
        break;
      case NM_SPARSE_BLOCK:
        M->matrix1 = (SparseBlockStructuredMatrix*) data;
        break;
      case NM_SPARSE:
        M->matrix2 = (NumericsSparseMatrix*) data;
        if (data)
        {
          if (M->matrix2->origin == NSM_UNKNOWN)
          {
            if (M->matrix2->triplet) { M->matrix2->origin = NSM_TRIPLET; }
            else if (M->matrix2->half_triplet) { M->matrix2->origin = NSM_HALF_TRIPLET; }
            else if (M->matrix2->csc) { M->matrix2->origin = NSM_CSC; }
            else if (M->matrix2->csr) { M->matrix2->origin = NSM_CSR; }
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
  switch (A->storageType)
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
    DEBUG_EXPR(NM_display(Atrans););
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


void NM_clearDense(NumericsMatrix* A)
{
  if (A->matrix0)
  {
    free(A->matrix0);
    A->matrix0 = NULL;
  }
}

void NM_clearSparseBlock(NumericsMatrix* A)
{
  if (A->matrix1)
  {
    SBM_free(A->matrix1);
    free(A->matrix1);
    A->matrix1 = NULL;
  }
}

void NM_clearSparse(NumericsMatrix* A)
{
  if (A->matrix2)
  {
    NSM_free(A->matrix2);
    free(A->matrix2);
    A->matrix2 = NULL;
  }
}

void NM_clearTriplet(NumericsMatrix* A)
{
  if (A->matrix2)
  {
    if (A->matrix2->triplet)
    {
      cs_spfree(A->matrix2->triplet);
      A->matrix2->triplet = NULL;
    }
  }
}

void NM_clearHalfTriplet(NumericsMatrix* A)
{
  if (A->matrix2)
  {
    if (A->matrix2->half_triplet)
    {
      cs_spfree(A->matrix2->half_triplet);
      A->matrix2->half_triplet = NULL;
    }
  }
}

void NM_clearCSC(NumericsMatrix* A)
{
  if (A->matrix2)
  {
    if (A->matrix2->csc)
    {
      cs_spfree(A->matrix2->csc);
      A->matrix2->csc = NULL;
    }
  }
}

void NM_clearCSCTranspose(NumericsMatrix* A)
{
  if (A->matrix2)
  {
    if (A->matrix2->trans_csc)
    {
      cs_spfree(A->matrix2->trans_csc);
      A->matrix2->trans_csc = NULL;
    }
  }
}

void NM_clearCSR(NumericsMatrix* A)
{
  if (A->matrix2)
  {
    if (A->matrix2->csr)
    {
      cs_spfree(A->matrix2->csr);
      A->matrix2->csr = NULL;
    }
  }
}

void NM_clearSparseStorage(NumericsMatrix *A)
{
  if (A->matrix2){
    A->matrix2->origin = NSM_UNKNOWN;
    if (A->matrix2->linearSolverParams)
      A->matrix2->linearSolverParams = NSM_linearSolverParams_free(A->matrix2->linearSolverParams);
  }
  NM_clearTriplet(A);
  NM_clearHalfTriplet(A);
  NM_clearCSC(A);
  NM_clearCSCTranspose(A);
  NM_clearCSR(A);

}


void NM_dense_to_sparse(const NumericsMatrix* const A, NumericsMatrix* B)
{
  assert(A->matrix0);
  assert(B->matrix2->triplet);
  for (int i = 0; i < A->size0; ++i)
  {
    for (int j = 0; j < A->size1; ++j)
    {
      CHECK_RETURN(CSparseMatrix_zentry(B->matrix2->triplet, i, j, A->matrix0[i + A->size0*j]));
    }
  }
}
int NM_to_dense(const NumericsMatrix* const A, NumericsMatrix* B)
{
  int info = 1;
  if (!B->matrix0)
  {
    B->matrix0 = (double *)calloc(A->size0*A->size1, sizeof(double));
  }
  else if (B->size0 != A->size0 || B->size0 != A->size0)
  {
    free(B->matrix0);
    B->matrix0 = (double *)calloc(A->size0*A->size1, sizeof(double));
  }

  assert(B->matrix0);

  B->size0 = A->size0;
  B->size1 = A->size1;
  B->storageType=NM_DENSE;

  switch (A->storageType)
  {
  case NM_DENSE:
  {
    NM_copy(A, B);
    info=0;
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    SBM_to_dense(A->matrix1, B->matrix0);
    info=0;
    break;
  }
  case NM_SPARSE:
  {
    assert(A->matrix2);
    info  = NSM_to_dense(A->matrix2, B->matrix0);
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

  return info;


}


void NM_copy_to_sparse(const NumericsMatrix* const A, NumericsMatrix* B)
{
  DEBUG_BEGIN("NM_copy_to_sparse(...)\n")
  assert(A);
  assert(B);
  B->size0 = A->size0;
  B->size1 = A->size1;

  assert(B->storageType == NM_SPARSE);
  if (!B->matrix2)
  {
    B->matrix2 = NSM_new();
  }

  switch (A->storageType)
  {
  case NM_DENSE:
  {
    B->matrix2->triplet = cs_spalloc(0,0,1,1,1);
    B->matrix2->origin = NSM_TRIPLET;
    NM_dense_to_sparse(A, B);
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
    break;
  }
  case NM_SPARSE:
  {
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
  switch (A->storageType)
  {
  case NM_DENSE:
  {
    if (B->matrix0)
    {
      if (sizeB < sizeA)
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

    break;
  }
  case NM_SPARSE_BLOCK:
  {
    if (!B->matrix1)
    {
      B->matrix1 = SBM_new();
    }

    SparseBlockStructuredMatrix* A_ = A->matrix1;
    SparseBlockStructuredMatrix* B_ = B->matrix1;

    SBM_copy(A_,B_,1);

    /* invalidations */
    NM_clearDense(B);
    NM_clearSparseStorage(B);

    break;
  }
  case NM_SPARSE:
  {
    CSparseMatrix* A_;
    CSparseMatrix* B_;

    if (!B->matrix2)
    {
      B->matrix2 = NSM_new();
    }

    B->matrix2->origin = A->matrix2->origin;

    switch (A->matrix2->origin)
    {
    case NSM_TRIPLET:
    {
      A_ = A->matrix2->triplet;

      if (!B->matrix2->triplet)
      {
        B->matrix2->triplet = cs_spalloc(A_->m, A_->n, A_->nzmax, 0, 1);
      }

      B_ = B->matrix2->triplet;
      break;
    }
    case NSM_HALF_TRIPLET:
    {
      A_ = A->matrix2->half_triplet;

      if (!B->matrix2->half_triplet)
      {
        B->matrix2->half_triplet = cs_spalloc(A_->m, A_->n, A_->nzmax, 0, 1);
      }

      B_ = B->matrix2->half_triplet;
      break;
    }
    case NSM_CSC:
    {
      assert (A->matrix2->csc);

      A_ = A->matrix2->csc;

      if (!B->matrix2->csc)
      {
        B->matrix2->csc = cs_spalloc(A_->m, A_->n, A_->nzmax, 0, 0);
      }

      B_ = B->matrix2->csc;
      break;
    }
    case NSM_CSR:
    {
      assert (A->matrix2->csr);

      A_ = A->matrix2->csr;

      if (!B->matrix2->csr)
      {
        NM_csr_alloc(B, A_->nzmax);
      }

      B_ = B->matrix2->csr;
      break;
    }
    default:
    {
      fprintf(stderr, "NM_copy :: error unknown origin %d for sparse matrix\n", A->matrix2->origin);
      exit(EXIT_FAILURE);
    }
    }
    CSparseMatrix_copy(A_, B_);


    /* invalidations */
    NM_clearDense(B);
    NM_clearSparseBlock(B);

    /* We remove diag_indx from B and  we copy it from A if it exists */
    if (numericsSparseMatrix(B)->diag_indx)
    {
      free(numericsSparseMatrix(B)->diag_indx);
      numericsSparseMatrix(B)->diag_indx=NULL;
    }
    if (A->matrix2->diag_indx)
    {
      numericsSparseMatrix(B)->diag_indx = (CS_INT*) malloc(A->size0 * sizeof(CS_INT));
      memcpy(numericsSparseMatrix(B)->diag_indx, A->matrix2->diag_indx, A->size0 * sizeof(CS_INT));
    }

    if (B_->nz >= 0)
    {
      NM_clearCSC(B);
      NM_clearCSCTranspose(B);
      NM_clearCSR(B);
    }
    else
    {
      NM_clearTriplet(B);
      NM_clearHalfTriplet(B);
      if (A->matrix2->origin == NSM_CSC) { NM_clearCSR(B); }
      else { NM_clearCSC(B); }
    }

    break;
  }
  }
  NM_internalData_copy(A, B);
  NM_MPI_copy(A, B);
  NM_MUMPS_copy(A, B);
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
  if (!numericsSparseMatrix(A)->triplet)
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

      if (A->matrix1)
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
                CHECK_RETURN(CSparseMatrix_zentry(A->matrix2->triplet, i + roffset, j + coffset,
                                                  A->matrix1->block[bn][i + j*inbr]));
              }
            }
          }
        }
      }
      else if (A->matrix0)
      {
        NM_dense_to_sparse(A, A);
      }
      else if (A->size0 > 0 || A->size1 > 0)
      {
        fprintf(stderr, "NM_triplet: sparse matrix cannot be constructed.\n");
        exit(EXIT_FAILURE);
      }
      break;
    }
    case NM_SPARSE:
    {
      switch (A->matrix2->origin)
      {
      case NSM_CSC:
      {
        assert(A->matrix2->csc);
        A->matrix2->triplet = NM_csc_to_triplet(A->matrix2->csc);
        break;
      }
      case NSM_CSR:
      {
        assert(A->matrix2->csr);
        A->matrix2->triplet = NM_csr_to_triplet(A->matrix2->csr);
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
  assert (A->matrix2->triplet);

  return A->matrix2->triplet;
}


CSparseMatrix* NM_half_triplet(NumericsMatrix* A)
{
  if (!numericsSparseMatrix(A)->half_triplet)
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

      if (A->matrix1)
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
                CHECK_RETURN(CSparseMatrix_symmetric_zentry(A->matrix2->half_triplet, i + roffset, j + coffset,
                                                            A->matrix1->block[bn][i + j*inbr]));
              }
            }
          }
        }
      }
      else if (A->matrix0)
      {
        fprintf(stderr, "NM_half_triplet: conversion is not implemented");
        exit(EXIT_FAILURE);
      }
      break;
    }
    case NM_SPARSE:
    {
      switch (A->matrix2->origin)
      {
      case NSM_TRIPLET:
      {
        A->matrix2->half_triplet = NM_csc_to_half_triplet(NM_csc(A));
        break;
      }
      case NSM_CSC:
      {
        assert(A->matrix2->csc);
        A->matrix2->half_triplet = NM_csc_to_half_triplet(A->matrix2->csc);
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
        NSM_UNKNOWN_ERR("NM_half_triplet", A->matrix2->origin);
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
  assert (A->matrix2->half_triplet);

  return A->matrix2->half_triplet;
}

CSparseMatrix* NM_csc(NumericsMatrix *A)
{
  DEBUG_BEGIN("NM_csc(NumericsMatrix *A)\n");
  assert(A);

  if(!numericsSparseMatrix(A)->csc)
  {
    assert(A->matrix2);
    switch (A->matrix2->origin)
    {
    case NSM_TRIPLET:
    case NSM_UNKNOWN:
    {
      /*  triplet -> csc with allocation */
      A->matrix2->csc = cs_compress(NM_triplet(A));
      break;
    }
    case NSM_CSR:
    {
      A->matrix2->csc = NM_csr_to_csc(NM_csr(A));
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
  DEBUG_END("NM_csc(NumericsMatrix *A)\n");
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
    switch (A->matrix2->origin)
    {
    case NSM_TRIPLET:
    case NSM_UNKNOWN:
    {
      /*  triplet -> csr with allocation */
      A->matrix2->csr = NM_triplet_to_csr(NM_triplet(A));
      break;
    }
    case NSM_CSR:
    {
      A->matrix2->csr = NM_csc_to_csr(NM_csr(A));
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
  return A->matrix2->csr;
}

/* Numerics Matrix wrapper  for y <- alpha A x + beta y */
void NM_gemv(const double alpha, NumericsMatrix* A, const double *x,
             const double beta, double *y)
{
  assert(A);
  assert(x);
  assert(y);

  switch (A->storageType)
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
  switch (A->storageType)
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

NumericsMatrix * NM_multiply(NumericsMatrix* A, NumericsMatrix* B)
{
  DEBUG_BEGIN("NM_multiply(...) \n")
  size_t storageType;

  NumericsMatrix * C = NM_new();

  /* should we copy the whole internal data ? */
  /*NM_internalData_copy(A, C);*/
  NM_MPI_copy(A, C);
  NM_MUMPS_copy(A, C);

  /* At the time of writing, we are able to transform anything into NM_SPARSE,
   * hence we use this format whenever possible */
  if (A->storageType == NM_SPARSE || B->storageType == NM_SPARSE || C->storageType == NM_SPARSE)
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
    if (check_mkl_lib() && (fabs(beta -1.) < 100*DBL_EPSILON))
    {
      if (!B->matrix2) NM_triplet(B);
      NumericsSparseMatrix* result = NM_MKL_spblas_gemm(0, A->matrix2, B->matrix2);
      assert(result);
      int size0 = C->size0;
      int size1 = C->size1;
      NM_free(C);
      NM_null(C);
      NM_fill(C, NM_SPARSE, size0, size1, result);
      NM_MKL_to_sparse_matrix(C);
      return;
    }
#endif
    DEBUG_EXPR(NM_display(A));
    DEBUG_EXPR(NM_display(B));
    DEBUG_EXPR(cs_print((const cs * ) NM_csc(A),0););
    DEBUG_EXPR(cs_print((const cs * ) NM_csc(B),0););
    assert(A->size1 == B->size0 && "NM_gemm :: A->size1 != B->size0 ");
    CSparseMatrix* C_csc = cs_multiply(NM_csc(A), NM_csc(B));
    DEBUG_EXPR(cs_print((const cs * ) C_csc,0););
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
  return C;
  DEBUG_END("NM_multiply(...) \n")
}

void NM_gemm(const double alpha, NumericsMatrix* A, NumericsMatrix* B,
             const double beta, NumericsMatrix* C)
{
  size_t storageType;

  /* At the time of writing, we are able to transform anything into NM_SPARSE,
   * hence we use this format whenever possible */
  if (A->storageType == NM_SPARSE || B->storageType == NM_SPARSE || C->storageType == NM_SPARSE)
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
    if (check_mkl_lib() && (fabs(beta -1.) < 100*DBL_EPSILON))
    {
      if (!B->matrix2) NM_triplet(B);
      if (!C->matrix2) NM_triplet(C);
      NumericsSparseMatrix* tmp_matrix = NM_MKL_spblas_gemm(0, A->matrix2, B->matrix2);
      assert(tmp_matrix);
      NumericsSparseMatrix* result = NM_MKL_spblas_add(0, alpha, tmp_matrix, C->matrix2);
      int size0 = C->size0;
      int size1 = C->size1;
      NM_free(C);
      NM_null(C);
      NM_fill(C, NM_SPARSE, size0, size1, result);
      NM_MKL_to_sparse_matrix(C);
      return;
    }
#endif
    DEBUG_EXPR(NM_display(A));
    DEBUG_EXPR(NM_display(B));
    DEBUG_EXPR(cs_print((const cs * ) NM_csc(A),0););
    DEBUG_EXPR(cs_print((const cs * ) NM_csc(B),0););
    assert(A->size1 == B->size0 && "NM_gemm :: A->size1 != B->size0 ");
    CSparseMatrix* tmp_matrix = cs_multiply(NM_csc(A), NM_csc(B));
    DEBUG_EXPR(cs_print((const cs * ) tmp_matrix,0););
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
  if (!A->internalData)
  {
    NM_internalData_new(A);
  }
  return A->internalData;
}




void* NM_iWork(NumericsMatrix* A, size_t size, size_t sizeof_elt)
{
  size_t bit_size = size * sizeof_elt;
  NM_internalData(A)->sizeof_elt =   sizeof_elt;

  if (!NM_internalData(A)->iWork)
  {
    assert(A->internalData->iWorkSize == 0);
    A->internalData->iWork = malloc(bit_size);
    A->internalData->iWorkSize = bit_size;

  }
  else
  {
    assert(A->internalData);

    if (bit_size > A->internalData->iWorkSize)
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
  if (!NM_internalData(A)->dWork)
  {
    assert(A->internalData);

    assert(A->internalData->dWorkSize == 0);
    A->internalData->dWork = (double *) malloc(size * sizeof(double));
    A->internalData->dWorkSize = size;
  }
  else
  {
    assert(A->internalData);

    if ( (size_t)
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

int NM_gesv_expert(NumericsMatrix* A, double *b, unsigned keep)
{

  DEBUG_BEGIN("NM_gesv_expert(NumericsMatrix* A, double *b, unsigned keep)\n");
  assert(A->size0 == A->size1);

  lapack_int info = 1;

  switch (A->storageType)
  {
  case NM_DENSE:
  {
    assert(A->matrix0);

    if (keep == NM_KEEP_FACTORS)
    {
      numerics_printf_verbose(2,"NM_gesv_expert, using LAPACK" );
      //double* wkspace = NM_dWork(A, A->size0*A->size1);
      lapack_int* ipiv = (lapack_int*)NM_iWork(A, A->size0, sizeof(lapack_int));
      DEBUG_PRINTF("iwork and dwork are initialized with size %i and %i\n",A->size0*A->size1,A->size0 );

      if (!NM_internalData(A)->isLUfactorized)
      {
        numerics_printf_verbose(2,"NM_gesv_expert, we compute factors and keep it" );
        DEBUG_PRINT("Start to call DGETRF for NM_DENSE storage\n");
        //cblas_dcopy_msan(A->size0*A->size1, A->matrix0, 1, wkspace, 1);
        DGETRF(A->size0, A->size1, A->matrix0, A->size0, ipiv, &info);
        DEBUG_PRINT("end of call DGETRF for NM_DENSE storage\n");
        if (info > 0)
        {
          if (verbose >= 2)
          {
            printf("NM_gesv: LU factorisation DGETRF failed. The %d-th diagonal element is 0\n", info);
          }
        }
        else if (info < 0)
        {
          fprintf(stderr, "NM_gesv: LU factorisation DGETRF failed. The %d-th argument has an illegal value, stopping\n", -info);
        }

        if (info) { NM_internalData_free(A); return info; }

        NM_internalData(A)->isLUfactorized = true;
      }
      DEBUG_PRINT("Start to call DGETRS for NM_DENSE storage\n");
      numerics_printf_verbose(2,"NM_gesv_expert, we solve with given factors" );
      DGETRS(LA_NOTRANS, A->size0, 1, A->matrix0, A->size0, ipiv, b, A->size0, &info);
      DEBUG_PRINT("End of call DGETRS for NM_DENSE storage\n");
      if (info < 0)
      {
        if (verbose >= 2)
        {
          printf("NM_gesv: dense LU solve DGETRS failed. The %d-th argument has an illegal value, stopping\n", -info);
        }
      }
    }
    else
    {
      double* mat;
      if (keep == NM_PRESERVE)
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
    switch (p->solver)
    {
    case NSM_CS_LUSOL:
      numerics_printf_verbose(2,"NM_gesv, using CSparse" );

      if (keep == NM_KEEP_FACTORS)
      {
        if (!(p->dWork && p->linear_solver_data))
        {
          assert(!NSM_workspace(p));
          assert(!NSM_linear_solver_data(p));
          assert(!p->solver_free_hook);

          p->solver_free_hook = &NSM_free_p;
          p->dWork = (double*) malloc(A->size1 * sizeof(double));
          p->dWorkSize = A->size1;
          CSparseMatrix_factors* cs_lu_A = (CSparseMatrix_factors*) malloc(sizeof(CSparseMatrix_factors));
          numerics_printf_verbose(2,"NM_gesv_expert, we compute factors and keep it" );
          CHECK_RETURN(CSparsematrix_lu_factorization(1, NM_csc(A), DBL_EPSILON, cs_lu_A));
          p->linear_solver_data = cs_lu_A;
        }

        numerics_printf_verbose(2,"NM_gesv, we solve with given factors" );
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
      if (verbose >= 2)
      {
        printf("NM_gesv: using MUMPS\n" );
      }
      if (!NM_MUMPS_id(A)->job || (NM_MUMPS_id(A)->job == -2))
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
      NM_MUMPS_set_problem(A, b);

      DMUMPS_STRUC_C* mumps_id = NM_MUMPS_id(A);

      if (keep != NM_KEEP_FACTORS|| mumps_id->job == -1)
      {
        NM_MUMPS(A, 6); /* analyzis,factorization,solve*/
      }
      else
      {
        NM_MUMPS(A, 3); /* solve */
      }

      info = mumps_id->info[0];

      /* MUMPS can return info codes with negative value */
      if (info)
      {
        fprintf(stderr,"NM_gesv: MUMPS fails : info(1)=%d, info(2)=%d\n", info, mumps_id->info[1]);
      }
      if (keep != NM_KEEP_FACTORS)
      {
        NM_MUMPS(A, -2);
      }
      if (!p->solver_free_hook)
      {
        p->solver_free_hook = &NM_MUMPS_free;
      }
      break;
    }
#endif /* WITH_MUMPS */

#ifdef WITH_UMFPACK
    case NSM_UMFPACK:
    {
      if (verbose >= 2)
      {
        printf("NM_gesv: using UMFPACK\n" );
      }

      NM_UMFPACK_WS* umfpack_ws = NM_UMFPACK_factorize(A);

      if (!umfpack_ws)
      {
        if (verbose > 1)
          fprintf(stderr, "NM_gesv: cannot factorize the matrix with UMFPACK\n");

        NM_UMFPACK_free(p);
        return -1;
      }

      CSparseMatrix* C = NM_csc(A);
      info = (int)UMFPACK_FN(wsolve) (UMFPACK_A, C->p, C->i, C->x, umfpack_ws->x, b, umfpack_ws->numeric, umfpack_ws->control, umfpack_ws->info, umfpack_ws->wi, umfpack_ws->wd);

      if (info)
      {
        UMFPACK_FN(report_status) (umfpack_ws->control, (CS_INT)info);
      }
      else
      {
        cblas_dcopy(C->n, umfpack_ws->x, 1, b, 1);
      }

      if (keep != NM_KEEP_FACTORS)
      {
        NM_UMFPACK_free(p);
      }
      else if (!p->solver_free_hook)
      {
        p->solver_free_hook = &NM_UMFPACK_free;
      }

      break;
    }

#endif /* WITH_UMFPACK */

#ifdef WITH_SUPERLU
    case NSM_SUPERLU:
    {
      if (verbose >= 2)
      {
        printf("NM_gesv: using SuperLU\n" );
      }

      NM_SuperLU_WS* superlu_ws = NM_SuperLU_factorize(A);

      if (!superlu_ws)
      {
        if (verbose > 1)
          fprintf(stderr, "NM_gesv: cannot factorize the matrix with SuperLU\n");

        NM_SuperLU_free(p);
        return -1;
      }

      info = NM_SuperLU_solve(A, b, superlu_ws);

      if (info)
      {
        fprintf(stderr, "NM_gesv: cannot solve the system with SuperLU\n");
//        SuperLU_FN(report_status) (superlu_ws->control, (CS_INT)info);
      }

      if (keep != NM_KEEP_FACTORS)
      {
        NM_SuperLU_free(p);
      }
      else if (!p->solver_free_hook)
      {
        p->solver_free_hook = &NM_SuperLU_free;
      }

      break;
    }

#endif /* WITH_SUPERLU */

#ifdef WITH_SUPERLU_MT
    case NSM_SUPERLU_MT:
    {
      if (verbose >= 2)
      {
        printf("NM_gesv: using SuperLU_MT\n" );
      }

      NM_SuperLU_MT_WS* superlu_mt_ws = NM_SuperLU_MT_factorize(A);

      if (!superlu_mt_ws)
      {
        if (verbose > 1)
          fprintf(stderr, "NM_gesv: cannot factorize the matrix with SuperLU_MT\n");

        NM_SuperLU_MT_free(p);
        return -1;
      }

      info = NM_SuperLU_MT_solve(A, b, superlu_mt_ws);

      if (info)
      {
        fprintf(stderr, "NM_gesv: cannot solve the system with SuperLU_MT\n");
//        SuperLU_MT_FN(report_status) (superlu_ws->control, (CS_INT)info);
      }

      if (keep != NM_KEEP_FACTORS)
      {
        NM_SuperLU_MT_free(p);
      }
      else if (!p->solver_free_hook)
      {
        p->solver_free_hook = &NM_SuperLU_MT_free;
      }

      break;
    }

#endif /* WITH_SUPERLU_MT */

#ifdef WITH_MKL_PARDISO
    case NSM_MKL_PARDISO:
    {
      if (verbose >= 2)
      {
        printf("NM_gesv: using MKL_PARDISO\n" );
      }

      NM_MKL_pardiso_WS* pardiso_ws = NM_MKL_pardiso_factorize(A);

      if (!pardiso_ws)
      {
        if (verbose > 1)
          fprintf(stderr, "NM_gesv: cannot factorize the matrix with MKL_PARDISO\n");

        NM_MKL_pardiso_free(p);
        return -1;
      }

      info = NM_MKL_pardiso_solve(A, b, pardiso_ws);

      if (keep != NM_KEEP_FACTORS)
      {
        NM_MKL_pardiso_free(p);
      }
      else if (!p->solver_free_hook)
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
    assert (0 && "NM_gesv unknown storageType");
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
  switch (A->storageType)
  {
  case NM_DENSE:
  {
    assert(A->matrix0);

    if (keep == NM_KEEP_FACTORS)
    {
      numerics_printf_verbose(2,"NM_posv_expert, using LAPACK" );
      //double* wkspace = NM_dWork(A, A->size0*A->size1);

      DEBUG_PRINTF("iwork and dwork are initialized with size %i and %i\n",A->size0*A->size1,A->size0 );

      if (!NM_internalData(A)->isLUfactorized)
      {
        numerics_printf_verbose(2,"NM_posv_expert, we compute factors and keep it" );
        DEBUG_PRINT("Start to call DPOTRF for NM_DENSE storage\n");
        //cblas_dcopy_msan(A->size0*A->size1, A->matrix0, 1, wkspace, 1);
        DPOTRF(LA_UP, A->size1, A->matrix0, A->size0, &info);
        DEBUG_PRINT("end of call DPOTRF for NM_DENSE storage\n");
        if (info > 0)
        {
          if (verbose >= 2)
          {
            printf("NM_posv: Cholesky factorisation DPOTRF failed. The %d-th diagonal element is 0\n", info);
          }
        }
        else if (info < 0)
        {
          fprintf(stderr, "NM_posv: Cholesky factorisation DPOTRF failed. The %d-th argument has an illegal value, stopping\n", -info);
        }

        if (info) { NM_internalData_free(A); return info; }

        NM_internalData(A)->isLUfactorized = true;
      }
      DEBUG_PRINT("Start to call DPOTRS for NM_DENSE storage\n");
      numerics_printf_verbose(2,"NM_posv_expert, we solve with given factors" );
      DPOTRS(LA_UP, A->size0, 1, A->matrix0, A->size0, b, A->size0, &info);
      DEBUG_PRINT("End of call DPOTRS for NM_DENSE storage\n");
      if (info < 0)
      {
        if (verbose >= 2)
        {
          printf("NM_posv: dense LU solve DGETRS failed. The %d-th argument has an illegal value, stopping\n", -info);
        }
      }
    }
    else
    {
      double* mat;
      if (keep == NM_PRESERVE)
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
      if (info > 0)
      {
        if (verbose >= 2)
        {
          printf("NM_posv: Cholesky solver DPOSV failed. The %d-th diagonal element is 0\n", info);
        }
      }
      else if (info < 0)
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
    switch (p->solver)
    {
    case NSM_CS_CHOLSOL:
      numerics_printf_verbose(2,"NM_posv, using CSparse cholsol" );

      if (keep == NM_KEEP_FACTORS)
      {
        if (!(p->dWork && p->linear_solver_data))
        {
          assert(!NSM_workspace(p));
          assert(!NSM_linear_solver_data(p));
          assert(!p->solver_free_hook);

          p->solver_free_hook = &NSM_free_p;
          p->dWork = (double*) malloc(A->size1 * sizeof(double));
          p->dWorkSize = A->size1;

          CSparseMatrix_factors* cs_chol_A = (CSparseMatrix_factors*) malloc(sizeof(CSparseMatrix_factors));

          numerics_printf_verbose(2,"NM_posv_expert, we compute factors and keep it" );
          CHECK_RETURN(CSparsematrix_chol_factorization(1, NM_csc(A),  cs_chol_A));

          p->linear_solver_data = cs_chol_A;
        }

        numerics_printf_verbose(2,"NM_posv, we solve with given factors" );
        info = !CSparseMatrix_chol_solve((CSparseMatrix_factors *)NSM_linear_solver_data(p), NSM_workspace(p), b);
      }
      else
      {
        numerics_printf_verbose(2,"NM_posv" );
        info = !cs_cholsol(1, NM_csc(A), b);
        if (info > 0)
        {
          numerics_printf("NM_posv: cs_cholsol failed. info = %i\n", info);
        }
        //DEBUG_EXPR(NV_display)
      }
      break;

#ifdef WITH_MUMPS
    case NSM_MUMPS:
    {
      if (verbose >= 2)
      {
        printf("NM_posv: using MUMPS\n" );
      }
      /* construction of lower triangular matrix */
      if (!NM_MUMPS_id(A)->job || (NM_MUMPS_id(A)->job == -2))
      {
        /* the mumps instance is initialized (call with job=-1) */
        NM_MUMPS_set_control_params(A);
        NM_MUMPS_set_sym(A, 2); /* general symmetric */
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

      NM_MUMPS_set_problem(A, b);

      DMUMPS_STRUC_C* mumps_id = NM_MUMPS_id(A);

      if (keep != NM_KEEP_FACTORS|| mumps_id->job == -1)
      {
        NM_MUMPS(A, 6); /* analyzis,factorization,solve*/
      }
      else
      {
        NM_MUMPS(A, 3); /* solve */
      }

      info = mumps_id->info[0];

      /* MUMPS can return info codes with negative value */
      if (info)
      {
        fprintf(stderr, "NM_posv: MUMPS fails : info(1)=%d, info(2)=%d\n", info, mumps_id->info[1]);
      }
      if (keep != NM_KEEP_FACTORS)
      {
        NM_MUMPS(A, -2);
      }
      if (!p->solver_free_hook)
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
    assert (0 && "NM_posv unknown storageType");
  }

  DEBUG_END("NM_posv_expert(NumericsMatrix* A, double *b, unsigned keep)\n");
  return (int)info;
}

int NM_gesv_expert_multiple_rhs(NumericsMatrix* A, double *b, unsigned int n_rhs, unsigned keep)
{
  assert(A->size0 == A->size1);
  int info = 0;
  for (unsigned int i = 0; i < n_rhs; ++i)
  {
    info = NM_gesv_expert(A, &b[A->size0*i],  keep);
  }
  return info;
}
NumericsMatrix* NM_inv(NumericsMatrix* A)
{

  DEBUG_BEGIN("NM_inv(NumericsMatrix* A, double *b, unsigned keep)\n");
  assert(A->size0 == A->size1);
  double * b = (double *) malloc(A->size0*sizeof(double));
  for (int i = 0; i < A->size0; ++i)
  {
    b[i]=0.0;
  }


  NumericsMatrix* Atmp = NM_new();
  NM_copy(A,Atmp);

  NumericsMatrix * Ainv  = NM_new();
  Ainv->size0 =  A->size0;
  Ainv->size1 =  A->size1;

  int info =-1;

  switch (A->storageType)
  {
  case NM_DENSE:
  {
    Ainv->storageType = NM_DENSE;
    Ainv->matrix0 = (double *)malloc(A->size0*A->size1*sizeof(double));
    for( int col_rhs =0; col_rhs < A->size1; col_rhs++ )
    {
      for (int i = 0; i < A->size0; ++i)
      {
        b[i]=0.0;
      }
      b[col_rhs] = 1.0;
      DEBUG_EXPR(NV_display(b,A->size1););
      //info = NM_gesv_expert(Atmp, b, NM_PRESERVE);
      info = NM_gesv_expert(Atmp, b, NM_KEEP_FACTORS);
      DEBUG_EXPR(NV_display(b,A->size1););
      if (info)
      {
        numerics_warning("NM_inv", "problem in NM_gesv_expert");
      }
      for (int i = 0; i < A->size0; ++i)
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

    for( int col_rhs =0; col_rhs < A->size1; col_rhs++ )
    {
      for (int i = 0; i < A->size0; ++i)
      {
        b[i]=0.0;
      }
      b[col_rhs] = 1.0;
      DEBUG_EXPR(NV_display(b,A->size1););
      //info = NM_gesv_expert(Atmp, b, NM_PRESERVE);
      info = NM_gesv_expert(Atmp, b, NM_KEEP_FACTORS);
      if (info)
      {
        numerics_warning("NM_inv", "problem in NM_gesv_expert");
      }
      for (int i = 0; i < A->size0; ++i)
      {
        CHECK_RETURN(CSparseMatrix_zentry(Ainv->matrix2->triplet, i, col_rhs, b[i]));
      }
    }
    break;
  }
  default:
    assert (0 && "NM_inv :  unknown storageType");
  }

  NM_free(Atmp);
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

  switch (A->storageType)
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
    assert (0 && "NM_inverse_diagonal_block_matrix_in_place :  unknown storageType");
  }


  DEBUG_BEGIN("NM_inverse_diagonal_block_matrix_in_place(NumericsMatrix* A)\n");
  return (int)info;
}


void NM_update_size(NumericsMatrix* A)
{
  switch (A->storageType)
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
    switch (A->matrix2->origin)
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
  for (int i =0 ; i < A->size1+1 ; i++)  Ap[i]= 0;
  //CS_INT* Ai = numericsSparseMatrix(A)->csc->i;
  //for (int i =0 ; i < nzmax ; i++)  Ai[i]= 0;
}


void NM_csr_alloc(NumericsMatrix* A, CS_INT nzmax)
{
  numericsSparseMatrix(A)->csr = cs_spalloc(A->size1, A->size0, nzmax, 1, 0);
  numericsSparseMatrix(A)->csr->nz = -2;
  numericsSparseMatrix(A)->csr->m =  A->size0;
  numericsSparseMatrix(A)->csr->n = A->size1;
}

void NM_triplet_alloc(NumericsMatrix* A, CS_INT nzmax)
{
  numericsSparseMatrix(A)->triplet = cs_spalloc(A->size0, A->size1, nzmax, 1, 1);
}

void NM_setSparseSolver(NumericsMatrix* A, unsigned solver_id)
{
  NSM_linearSolverParams(A)->solver = (NSM_linear_solver)solver_id;
}



int NM_check(const NumericsMatrix* const A)
{
  int info = 0;
  if (!A->matrix2) return info;

  if (A->matrix2->csc) info = CSparseMatrix_check_csc(A->matrix2->csc);
  if (A->matrix2->csr) info = info ? info : CSparseMatrix_check_csc(A->matrix2->csr);
  if (A->matrix2->triplet) info = info ? info : CSparseMatrix_check_triplet(A->matrix2->triplet);
  return info;
}


size_t NM_nnz(const NumericsMatrix* M)
{
  switch (M->storageType)
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


  switch (A->storageType)
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


  NumericsMatrix * Atrans = NM_transpose(A);
  NumericsMatrix * AplusATrans = NM_add(1/2.0, A, -1/2.0, Atrans);
  double norm_inf = NM_norm_inf(AplusATrans);
  NM_free(Atrans); free(Atrans);
  NM_free(AplusATrans);  free(AplusATrans);

  if (norm_inf <= DBL_EPSILON*10)
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
  for (int i =0; i < n ; i++)
  {
    for (int j =0 ; j < m ; j++)
    {
      d= fmax (d,fabs(NM_get_value(A,i,j)-NM_get_value(A,j,i)));
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
  srand(time(NULL));
  for (int i = 0; i < n ; i++)
  {
    q[i] = (rand()/(double)RAND_MAX);
    DEBUG_PRINTF("q[%i] = %e \t",i, q[i]);
  }
  double norm = cblas_dnrm2(n , q, 1);
  cblas_dscal(m, 1.0/norm, q, 1);
  DEBUG_PRINT("\n");

  NM_gemv(1, A, q, 0.0, z);

  int k =0;
  while ((fabs((eig-eig_old)/eig_old) > tol) && k < itermax)
  {
    norm = cblas_dnrm2(n , z, 1);
    cblas_dscal(m, 1.0/norm, z, 1);
    cblas_dcopy(n , z  , 1 , q, 1);

    NM_gemv(1.0, A, q, 0.0, z);

    eig_old=eig;
    eig = cblas_ddot(n, q, 1, z, 1);

    DEBUG_PRINTF("eig[%i] = %32.24e \t error = %e\n",k, eig, fabs(eig-eig_old)/eig_old );
    k++;
  }
  free(q);
  free(z);
  return eig;
}
