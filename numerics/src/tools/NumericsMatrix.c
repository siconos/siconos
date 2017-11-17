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
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#include "csparse.h"
#include "NumericsMatrix.h"
#include "NumericsSparseMatrix.h"
#include "SiconosCompat.h"
#include "SparseBlockMatrix.h"
#include "NumericsMatrix_private.h"
#include "NumericsMatrix.h"
#include "NM_conversions.h"
#include "SiconosLapack.h"
#include "numerics_verbose.h"
#include "sanitizer.h"
/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"

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
      cs_aaxpy(alpha, NM_csc(A), x, beta, y);
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
    if (A->matrix2->origin == NS_CSR)
    {
      M = NM_csr(A);
    }
    else
    {
      M = NM_csc_trans(A);
    }

    csi* Mp = M->p;
    csi* Mi = M->i;
    double* Mx = M->x;

    for (size_t i = 0, j = row_start; i < sizeY; ++i, ++j)
    {
      for (csi p = Mp[j]; p < Mp[j+1]; ++p)
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
    if (A->matrix2->origin == NS_CSR)
    {
      M = NM_csr(A);
    }
    else
    {
      M = NM_csc_trans(A);
    }

    csi* Mp = M->p;
    csi* Mi = M->i;
    double* Mx = M->x;

    for (size_t i = 0, j = row_start; i < 3; ++i, ++j)
    {
      for (csi p = Mp[j]; p < Mp[j+1]; ++p)
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
    if (A->matrix2->origin == NS_CSR)
    {
      M = NM_csr(A);
    }
    else
    {
      M = NM_csc_trans(A);
    }

    csi* Mp = M->p;
    csi* Mi = M->i;
    double* Mx = M->x;

    for (csi p = Mp[row_start]; p < Mp[row_start+1]; ++p)
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
  case NM_SPARSE:
  {
    assert(M->matrix2);
    switch (M->matrix2->origin)
    {
    case NS_TRIPLET:
    {
      assert(M->matrix2->triplet);
      CHECK_RETURN(cs_zentry(M->matrix2->triplet, i, j, val));
      break;
    }
    case NS_CSC:
    {
      assert(M->matrix2->csc);
      CHECK_RETURN(cs_zentry(NM_triplet(M), i, j, val));
      M->matrix2->origin= NS_TRIPLET;
      NM_clearCSC(M);
      //NM_display(M);
      //NM_csc_alloc(M,M->matrix2->triplet->nz);
      //M->matrix2->origin= NS_CSC;
      NM_csc(M);

      M->matrix2->origin= NS_CSC;
      NM_clearTriplet(M);

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
    case NS_TRIPLET:
    {
      assert(M->matrix2->triplet);
      assert(i <M->matrix2->triplet->m );
      assert(j <M->matrix2->triplet->n );
      csi * Mi =   M->matrix2->triplet->i;
      csi * Mp =   M->matrix2->triplet->p;
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
    case NS_CSC:
    {
      assert(M->matrix2->csc);
      csi * Mi =   M->matrix2->csc->i;
      csi * Mp =   M->matrix2->csc->p;
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
  }
  default:
    fprintf(stderr, "NM_get_value ::  unknown matrix storage = %d\n", M->storageType);
  }
return 0.0;

}

bool NM_equal(NumericsMatrix* A, NumericsMatrix* B)
{
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
      //printf("error = %e\n", fabs(NM_get_value(A, i, j) - NM_get_value(B, i, j)));
      if (fabs(NM_get_value(A, i, j) - NM_get_value(B, i, j)) >= DBL_EPSILON) return false;
    }
  }
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
    case NS_TRIPLET:
    {
      printf("========== origin =  NS_TRIPLET\n");
      break;
    }
    case NS_CSC:
    {
      printf("========== origin =  NS_CSC\n");
      break;
    }
    case NS_CSR:
    {
      printf("========== origin =  NS_CSR\n");
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
  DEBUG_PRINT("\n  ========== printInFile(const NumericsMatrix* const m, FILE* file) start\n");

  if (! m)
  {
    fprintf(stderr, "Numerics, NumericsMatrix printInFile failed, NULL input.\n");
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
    assert(m->matrix2);
    if (m->matrix2->triplet)
    {
      cs_printInFile(m->matrix2->triplet, 0, file);
    }
    else if (m->matrix2->csc)
    {
      cs_printInFile(m->matrix2->csc, 0, file);
    }
    else if (m->matrix2->trans_csc)
    {
      cs_printInFile(m->matrix2->trans_csc, 0, file);
    }
    else
    {
      fprintf(stderr, "display for sparse matrix: no matrix found!\n");
    }
    break;
  }
  default:
  {
    fprintf(stderr, "Numerics, NumericsMatrix printInFile failed, unknown storage type .\n");
    exit(EXIT_FAILURE);
  }
  }
  DEBUG_PRINT("\n  ========== printInFile(const NumericsMatrix* const m, FILE* file) end\n");

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

void NM_read_in_file(NumericsMatrix* const m, FILE *file)
{
  if (! m)
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
}


int NM_new_from_file(NumericsMatrix* const m, FILE *file)
{
  if (! m)
  {
    fprintf(stderr, "Numerics, NumericsMatrix NM_new_from_file failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }

  if (m->storageType >= 0)
  {
    NM_free(m);
  }

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

    data = malloc(size1 * size0 * sizeof(double));
    double* data_d = (double*) data;

    for (size_t i = 0; i < size1 * size0; ++i)
    {
      CHECK_IO(fscanf(file, "%le ", &(data_d[i])), &info);
    }
  }
  else if (storageType == NM_SPARSE_BLOCK)
  {
    data = SBM_new();
    SBM_new_from_file((SparseBlockStructuredMatrix*)data, file);
  }
  else
  {
    fprintf(stderr, "newFromFile :: storageType %d not supported yet", storageType);
  }

  NM_fill(m, storageType, (int)size0, (int)size1, data);

  return info;
}

void NM_read_in_filename(NumericsMatrix* const m, const char *filename)
{
  FILE* finput = fopen(filename, "r");
  NM_new_from_file(m, finput);
  fclose(finput);
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
    int diagPos = SBM_get_position_diagonal_block(M->matrix1, block_row_nb);
    (*Block) = M->matrix1->block[diagPos];
    break;
  }
  case NM_SPARSE:
  {
    NM_sparse_extract_block(M, *Block, start_row, start_row, size, size);
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
    int diagPos = SBM_get_position_diagonal_block(M->matrix1, block_row_nb);
    (*Block) = M->matrix1->block[diagPos];
    break;
  }
  case NM_SPARSE:
  {
    size_t start_row = (size_t)block_row_nb + block_row_nb + block_row_nb;
    NM_sparse_extract_block(M, *Block, start_row, start_row, 3, 3);
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
    int diagPos = SBM_get_position_diagonal_block(M->matrix1, block_row_nb);
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
    NM_sparse_extract_block(M, *Block, start_row, start_row, 3, 3);
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
      int diagPos = SBM_get_position_diagonal_block(M->matrix1, ic);
      M->matrix1->block[diagPos][0] += alpha;
      M->matrix1->block[diagPos][4] += alpha;
      M->matrix1->block[diagPos][8] += alpha;
    }
    break;
  }
  case NM_SPARSE:
  {
    csi* diag_indices = NM_sparse_diag_indices(M);
    double* Mx = NM_sparse_data(M->matrix2);
    for (size_t i = 0; i < n; ++i) Mx[diag_indices[i]] += alpha;

    break;
  }
  default:
    printf("NM_add_to_diag3 :: unsupported matrix storage %d", M->storageType);
    exit(EXIT_FAILURE);
  }
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
      data = newNumericsSparseMatrix();
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

NumericsMatrix* NM_create(int storageType, int size0, int size1)
{
  assert(size0 > 0);
  assert(size1 > 0);
  NumericsMatrix* M = NM_new();

  void* data;

  switch (storageType)
  {
    case NM_DENSE:
      data = malloc(size0*size1*sizeof(double));
      break;
    case NM_SPARSE_BLOCK:
      data = SBM_new();
      break;
    case NM_SPARSE:
      data = newNumericsSparseMatrix();
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
  assert(size0 > 0);
  assert(size1 > 0);
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
          if (M->matrix2->origin == NS_UNKNOWN)
          {
            if (M->matrix2->triplet) { M->matrix2->origin = NS_TRIPLET; }
            else if (M->matrix2->csc) { M->matrix2->origin = NS_CSC; }
            else if (M->matrix2->csr) { M->matrix2->origin = NS_CSR; }
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
    freeNumericsSparseMatrix(A->matrix2);
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
  if (A->matrix2) A->matrix2->origin = NS_UNKNOWN;
  NM_clearTriplet(A);
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
      CHECK_RETURN(cs_zentry(B->matrix2->triplet, i, j, A->matrix0[i + A->size0*j]));
    }
  }
}
void NM_to_dense(const NumericsMatrix* const A, NumericsMatrix* B)
{
  assert(B->matrix0);
  for (int i = 0; i < A->size0; ++i)
  {
    for (int j = 0; j < A->size1; ++j)
    {
      B->matrix0[i+ j*A->size0] = NM_get_value((NumericsMatrix*) A,i,j);
    }
  }
}


void NM_copy_sparse(const CSparseMatrix* const A, CSparseMatrix* B)
{
  assert (A);
  assert (B);

  if (B->nzmax < A->nzmax)
  {
    B->x = (double *) realloc(B->x, A->nzmax * sizeof(double));
    B->i = (csi *) realloc(B->i, A->nzmax * sizeof(csi));
  }
  else if (!(B->x))
  {
    B->x = (double *) malloc(A->nzmax * sizeof(double));
  }

  if (A->nz >= 0)
  {
    /* triplet */
    B->p = (csi *) realloc(B->p, A->nzmax * sizeof(csi));
  }
  else if ((A->nz == -1) && (B->n < A->n))
  {
    /* csc */
    B->p = (csi *) realloc(B->p, (A->n + 1) * sizeof(csi));
  }
  else if ((A->nz == -2) && (B->m < A->m))
  {
    /* csr */
    B->p = (csi *) realloc(B->p, (A->m + 1) * sizeof(csi));
  }


  B->nzmax = A->nzmax;
  B->nz = A->nz;
  B->m = A->m;
  B->n = A->n;

  memcpy(B->x, A->x, A->nzmax * sizeof(double));
  memcpy(B->i, A->i, A->nzmax * sizeof(csi));

  size_t size_cpy = -1;
  if (A->nz >= 0) { size_cpy = A->nzmax; }
  else if (A->nz == -1) { size_cpy = A->n + 1; }
  else if (A->nz == -2) { size_cpy = A->m + 1; }

  memcpy(B->p, A->p, size_cpy * sizeof(csi));
}

void NM_sparse_extract_block(NumericsMatrix* M, double* blockM, size_t pos_row, size_t pos_col,
                             size_t block_row_size, size_t block_col_size)
{
  assert(M);
  assert(M->storageType == NM_SPARSE);
  assert(blockM);
  assert(pos_row < (size_t)M->size0);
  assert(pos_col < (size_t)M->size1);
  assert(block_row_size > 0 && block_row_size + pos_row <= (unsigned long int)M->size0);
  assert(block_col_size > 0 && block_col_size + pos_col <= (unsigned long int)M->size1);

  assert(M->matrix2);

  /* Clear memory */
  memset(blockM, 0, block_row_size*block_col_size * sizeof(double));

//  switch (Msparse->origin)
  {
//  case NS_CSC:
  {
    CSparseMatrix* Mcsc = NM_csc(M);
    assert(Mcsc);
    csi* Mp = Mcsc->p;
    csi* Mi = Mcsc->i;
    double* Mx = Mcsc->x;
    for (size_t j = pos_col; j < pos_col + block_col_size; ++j)
    {
      for (csi p = Mp[j]; p < Mp[j+1]; ++p)
      {
        csi row_nb = Mi[p];
        if (row_nb >= (csi) pos_row)
        {
          if (row_nb >= (csi)(pos_row + block_row_size))
          {
            break;
          }
          else
          {
            blockM[(j-pos_col)*block_col_size + row_nb - pos_row] = Mx[p];
          }
        }
      }
    }
//    break;
  }
//  default:
//    printf("NM_sparse_extract_block :: unsupported matrix type %d\n", Msparse->origin);
//    exit(EXIT_FAILURE);
  }
}

void NM_copy_to_sparse(const NumericsMatrix* const A, NumericsMatrix* B)
{
  assert(A);
  assert(B);
  B->size0 = A->size0;
  B->size1 = A->size1;

  assert(B->storageType == NM_SPARSE);
  if (!B->matrix2)
  {
    B->matrix2 = newNumericsSparseMatrix();
  }

  switch (A->storageType)
  {
  case NM_DENSE:
  {
    B->matrix2->triplet = cs_spalloc(0,0,1,1,1);
    B->matrix2->origin = NS_TRIPLET;
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
}

void NM_copy(const NumericsMatrix* const A, NumericsMatrix* B)
{
  assert(A);
  assert(B);
  int sizeA = A->size0 * A->size1;
  int sizeB = B->size0 * B->size1;
  B->size0 = A->size0;
  B->size1 = A->size1;

  if (B->storageType >= 0 && B->storageType != A->storageType)
  {
    NM_internalData_free(B);
  }
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
    int need_blocks = 0;

    SparseBlockStructuredMatrix* A_ = A->matrix1;
    SparseBlockStructuredMatrix* B_ = B->matrix1;

    if (B_)
    {
      if (B_->nbblocks < A_->nbblocks)
      {
        need_blocks = 1;
        for (unsigned i=0; i<B_->nbblocks; ++i)
        {
          free(B_->block [i]);
          B_->block [i] = NULL;
        }
        B_->block = (double **) realloc(B_->block, A_->nbblocks * sizeof(double *));
      }
      B_->nbblocks = A_->nbblocks;

      if (B_->blocknumber0 < A_->blocknumber0)
      {
        B_->blocksize0 = (unsigned int*) realloc(B_->blocksize0, A_->blocknumber0 * sizeof(unsigned int));
      }
      B_->blocknumber0 = A_->blocknumber0;

      if (B_->blocknumber1 < A_->blocknumber1)
      {
        B_->blocksize1 = (unsigned int*) realloc(B_->blocksize1, A_->blocknumber1 * sizeof(unsigned int));
      }
      B_->blocknumber1 = A_->blocknumber1;

      if (B_->filled1 < A_->filled1)
      {
        B_->index1_data = (size_t*) realloc(B_->index1_data, A_->filled1 * sizeof(size_t));
      }
      B_->filled1 = A_->filled1;

      if (B_->filled2 < A_->filled2)
      {
        B_->index2_data = (size_t*) realloc(B_->index2_data, A_->filled2 * sizeof(size_t));
      }
      B_->filled2 = A_->filled2;
    }
    else
    {
      B->matrix1 = SBM_new();
      B_ = B->matrix1;

      B_->block = (double **) malloc(A_->nbblocks * sizeof(double *));
      B_->nbblocks = A_->nbblocks;

      B_->blocksize0 = (unsigned int*) malloc(A_->blocknumber0 * sizeof(unsigned int));
      B_->blocknumber0 = A_->blocknumber0;

      B_->blocksize1 = (unsigned int*) malloc(A_->blocknumber1 * sizeof(unsigned int));
      B_->blocknumber1 = A_->blocknumber1;

      B_->index1_data = (size_t*) malloc(A_->filled1 * sizeof(size_t));
      B_->filled1 = A_->filled1;

      B_->index2_data = (size_t*) malloc(A_->filled2 * sizeof(size_t));
      B_->filled2 = A_->filled2;
    }

    memcpy(B_->blocksize0, A_->blocksize0, A_->blocknumber0 * sizeof(unsigned int));
    memcpy(B_->blocksize1, A_->blocksize1, A_->blocknumber1 * sizeof(unsigned int));
    memcpy(B_->index1_data, A_->index1_data, A_->filled1 * sizeof(size_t));
    memcpy(B_->index2_data, A_->index2_data, A_->filled2 * sizeof(size_t));

    /* cf SBM_copy */
    unsigned int currentRowNumber ;
    size_t colNumber;
    unsigned int nbRows, nbColumns;
    for (currentRowNumber = 0 ; currentRowNumber < A_->filled1 - 1; ++currentRowNumber)
    {
      for (size_t blockNum = A_->index1_data[currentRowNumber];
           blockNum < A_->index1_data[currentRowNumber + 1]; ++blockNum)
      {
        assert(blockNum < A_->filled2);
        colNumber = A_->index2_data[blockNum];
        /* Get dim. of the current block */
        nbRows = A_->blocksize0[currentRowNumber];
        if (currentRowNumber != 0)
          nbRows -= A_->blocksize0[currentRowNumber - 1];
        nbColumns = A_->blocksize1[colNumber];

        if (colNumber != 0)
          nbColumns -= A_->blocksize1[colNumber - 1];

        if (need_blocks)
        {
          B_->block[blockNum] = (double*)malloc(nbRows * nbColumns * sizeof(double));
        }

        for (unsigned int i = 0; i < nbRows * nbColumns; i++)
        {
          B_->block[blockNum] [i] = A_->block[blockNum] [i] ;
        }
      }
    }

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
      B->matrix2 = newNumericsSparseMatrix();
    }

    B->matrix2->origin = A->matrix2->origin;

    switch (A->matrix2->origin)
    {
    case NS_TRIPLET:
    {
      A_ = A->matrix2->triplet;

      if (!B->matrix2->triplet)
      {
        B->matrix2->triplet = cs_spalloc(A_->m, A_->n, A_->nzmax, 0, 1);
      }

      B_ = B->matrix2->triplet;
      break;
    }
    case NS_CSC:
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
    case NS_CSR:
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

    NM_copy_sparse(A_, B_);

    /* invalidations */
    NM_clearDense(B);
    NM_clearSparseBlock(B);

    if (B_->nz >= 0)
    {
      NM_clearCSC(B);
      NM_clearCSCTranspose(B);
      NM_clearCSR(B);
    }
    else
    {
      NM_clearTriplet(B);
      if (A->matrix2->origin == NS_CSC) { NM_clearCSR(B); }
      else { NM_clearCSC(B); }
    }

    break;
  }
  }
}

NumericsSparseMatrix* NM_sparse(NumericsMatrix* A)
{
  if(!A->matrix2)
  {
    A->matrix2 = newNumericsSparseMatrix();
  }
  return A->matrix2;
}

NumericsSparseLinearSolverParams* NM_linearSolverParams(NumericsMatrix* A)
{
  if(!NM_sparse(A)->linearSolverParams)
  {
    NM_sparse(A)->linearSolverParams = newNumericsSparseLinearSolverParams();
  }
  return NM_sparse(A)->linearSolverParams;
}


CSparseMatrix* NM_triplet(NumericsMatrix* A)
{
  if (!NM_sparse(A)->triplet)
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

      A->matrix2->origin = NS_TRIPLET;

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
                CHECK_RETURN(cs_zentry(A->matrix2->triplet, i + roffset, j + coffset,
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
      else
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
      case NS_CSC:
      {
        assert(A->matrix2->csc);
        A->matrix2->triplet = NM_csc_to_triplet(A->matrix2->csc);
        break;
      }
      case NS_CSR:
      {
        assert(A->matrix2->csr);
        A->matrix2->triplet = NM_csr_to_triplet(A->matrix2->csr);
        break;
      }
      default:
      case NS_UNKNOWN:
      {
        NS_UNKNOWN_ERR("NM_triplet", A->matrix2->origin);
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

CSparseMatrix* NM_csc(NumericsMatrix *A)
{
  assert(A);

  if(!NM_sparse(A)->csc)
  {
    assert(A->matrix2);
    switch (A->matrix2->origin)
    {
    case NS_TRIPLET:
    case NS_UNKNOWN:
    {
      /*  triplet -> csc with allocation */
      A->matrix2->csc = cs_compress(NM_triplet(A));
      break;
    }
    case NS_CSR:
    {
      A->matrix2->csc = NM_csr_to_csc(NM_csr(A));
      break;
    }
    default:
    {
      NS_UNKNOWN_ERR("NM_csc", A->matrix2->origin);
      exit(EXIT_FAILURE);
    }
    }

    assert(A->matrix2->csc);
    NM_clearCSCTranspose(A);
  }
  return A->matrix2->csc;
}

CSparseMatrix* NM_csc_trans(NumericsMatrix* A)
{
  if(!NM_sparse(A)->trans_csc)
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

  if(!NM_sparse(A)->csr)
  {
    assert(A->matrix2);
    switch (A->matrix2->origin)
    {
    case NS_TRIPLET:
    case NS_UNKNOWN:
    {
      /*  triplet -> csr with allocation */
      A->matrix2->csr = NM_triplet_to_csr(NM_triplet(A));
      break;
    }
    case NS_CSR:
    {
      A->matrix2->csr = NM_csc_to_csr(NM_csr(A));
      break;
    }
    default:
    {
      NS_UNKNOWN_ERR("NM_csr", A->matrix2->origin);
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
    // if possible use the much simpler version provided by CSparse
    // Also at the time of writing, cs_aaxpy is bugged --xhub
    bool beta_check = false;
    if (fabs(beta) < 100*DBL_EPSILON)
    {
      memset(y, 0, A->size0*sizeof(double));
      beta_check = true;
    }
    if (fabs(alpha - 1.) < 100*DBL_EPSILON && (fabs(beta - 1.) < 100*DBL_EPSILON || beta_check))
    {
      CHECK_RETURN(cs_gaxpy(NM_csc(A), x, y));
    }
    else
    {
      CHECK_RETURN(cs_aaxpy(alpha, NM_csc(A), x, beta, y));
    }
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
      /* if possible use the much simpler version provided by CSparse
       Also at the time of writing, cs_aaxpy is bugged --xhub */
      bool beta_check = false;
      if (fabs(beta) < 100*DBL_EPSILON)
      {
        memset(y, 0, A->size0*sizeof(double));
        beta_check = true;
      }
      if (fabs(alpha - 1.) < 100*DBL_EPSILON && (fabs(beta - 1.) < 100*DBL_EPSILON || beta_check))
      {
        CHECK_RETURN(cs_gaxpy(NM_csc_trans(A), x, y));
      }
      else
      {
        CHECK_RETURN(cs_aaxpy(alpha, NM_csc_trans(A), x, beta, y));
      }
      break;
    }
  default:
    {
      assert(0 && "NM_tgemv unknown storageType");
    }
  }
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
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, A->size0, B->size1, B->size1, alpha, A->matrix0, A->size0, B->matrix0, B->size0, beta, C->matrix0, A->size0);

    NM_clearSparseBlock(C);
    NM_clearSparseStorage(C);

    break;
  }
  case NM_SPARSE_BLOCK:
  {
    assert(A->matrix1);
    assert(B->matrix1);
    assert(C->matrix1);

    SBM_alloc_for_gemm(A->matrix1, B->matrix1, C->matrix1);
    SBM_gemm(alpha, A->matrix1, B->matrix1, beta, C->matrix1);

    NM_clearDense(C);
    NM_clearSparseStorage(C);
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

    CSparseMatrix* tmp_matrix = cs_multiply(NM_csc(A), NM_csc(B));
    DEBUG_EXPR(cs_print((const cs * ) tmp_matrix,0););
    assert(tmp_matrix && "NM_gemm :: cs_multiply failed");
    NM_sparse_fix_csc(tmp_matrix);

    CSparseMatrix* result = cs_add(tmp_matrix, NM_csc(C), alpha, beta);
    assert(result && "NM_gemm :: cs_add failed");
    NM_sparse_fix_csc(result);

    cs_spfree(tmp_matrix);
    NM_clearDense(C);
    NM_clearSparseBlock(C);
    NM_clearSparseStorage(C);

    NM_sparse(C)->csc = result;
    C->size0 = (int)C->matrix2->csc->m;
    C->size1 = (int)C->matrix2->csc->n;
    NM_sparse(C)->origin = NS_CSC;
    break;
  }
  default:
  {
    assert(0 && "NM_gemm unknown storageType");
  }
  }
}

NumericsMatrixInternalData* NM_internalData(NumericsMatrix* A)
{
  if (!A->internalData)
  {
    NM_alloc_internalData(A);
    A->internalData->iWork = NULL;
    A->internalData->iWorkSize = 0;
    A->internalData->dWork = NULL;
    A->internalData->dWorkSize = 0;
    A->internalData->isLUfactorized = 0;
  }
  return A->internalData;
}

void* NM_iWork(NumericsMatrix* A, size_t size, size_t sizeof_elt)
{
  size_t bit_size = size * sizeof_elt;
  if (!NM_internalData(A)->iWork)
  {
    assert(A->internalData);

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

      //double* wkspace = NM_dWork(A, A->size0*A->size1);
      lapack_int* ipiv = (lapack_int*)NM_iWork(A, A->size0, sizeof(lapack_int));
      DEBUG_PRINTF("iwork and dwork are initialized with size %i and %i\n",A->size0*A->size1,A->size0 );

      if (!NM_internalData(A)->isLUfactorized)
      {
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
    NumericsSparseLinearSolverParams* p = NM_linearSolverParams(A);
    switch (p->solver)
    {
    case NS_CS_LUSOL:
      if (verbose >= 2)
      {
        printf("NM_gesv: using CSparse\n" );
      }
      if (keep == NM_KEEP_FACTORS)
      {
        if (!(p->dWork && p->solver_data))
        {
          assert(!NM_sparse_workspace(p));
          assert(!NM_sparse_solver_data(p));
          assert(!p->solver_free_hook);

          p->solver_free_hook = &NM_sparse_free;
          p->dWork = (double*) malloc(A->size1 * sizeof(double));
          p->dWorkSize = A->size1;
          cs_lu_factors* cs_lu_A = (cs_lu_factors*) malloc(sizeof(cs_lu_factors));
          CHECK_RETURN(cs_lu_factorization(1, NM_csc(A), DBL_EPSILON, cs_lu_A));
          p->solver_data = cs_lu_A;
        }

        info = !cs_solve((cs_lu_factors *)NM_sparse_solver_data(p), NM_sparse_workspace(p), b);
      }
      else
      {
        info = !cs_lusol(1, NM_csc(A), b, DBL_EPSILON);
      }
      break;

#ifdef WITH_MUMPS
    case NS_MUMPS:
    {
      if (verbose >= 2)
      {
        printf("NM_gesv: using MUMPS\n" );
      }
      /* the mumps instance is initialized (call with job=-1) */
      DMUMPS_STRUC_C* mumps_id = NM_MUMPS_id(A);

      mumps_id->rhs = b;

      if (keep != NM_KEEP_FACTORS|| mumps_id->job == -1)
      {
        mumps_id->job = 6;
      }
      else
      {
        mumps_id->job = 3;
      }


      /* compute the solution */
      dmumps_c(mumps_id);

      info = mumps_id->info[0];

      /* MUMPS can return info codes with negative value */
      if (info)
      {
        if (verbose > 0)
        {
          printf("NM_gesv: MUMPS fails : info(1)=%d, info(2)=%d\n", info, mumps_id->info[1]);
        }
      }
      if (keep != NM_KEEP_FACTORS)
      {
        NM_MUMPS_free(p);
      }
      else if (!p->solver_free_hook)
      {
        p->solver_free_hook = &NM_MUMPS_free;
      }

      break;
    }
#endif /* WITH_MUMPS */

#ifdef WITH_UMFPACK
    case NS_UMFPACK:
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
        UMFPACK_FN(report_status) (umfpack_ws->control, (csi)info);
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
    case NS_SUPERLU:
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
//        SuperLU_FN(report_status) (superlu_ws->control, (csi)info);
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
    case NS_SUPERLU_MT:
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
//        SuperLU_MT_FN(report_status) (superlu_ws->control, (csi)info);
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
    case NS_MKL_PARDISO:
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
int NM_inv(NumericsMatrix* A, NumericsMatrix* Ainv)
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

   int info =-1;

  switch (A->storageType)
  {
  case NM_DENSE:
  {
    assert (0 && "NM_inv :  not implemented");
    break;
  }
  case NM_SPARSE_BLOCK: /* sparse block -> triplet -> csc */
  case NM_SPARSE:
  {
    // We clear Ainv
    NM_clearSparse(Ainv);
    Ainv->storageType = NM_SPARSE;
    Ainv->size0 = A->size0;
    Ainv->size1 = A->size1;
    NM_triplet_alloc(Ainv,  A->size0);
    Ainv->matrix2->origin = NS_TRIPLET;

    for( int col_rhs =0; col_rhs < A->size1; col_rhs++ )
    {
      if (col_rhs >0) b[col_rhs-1] = 0.0;
      b[col_rhs] = 1.0;
      //info = NM_gesv_expert(Atmp, b, NM_PRESERVE);
      info = NM_gesv_expert(Atmp, b, NM_KEEP_FACTORS);

      for (int i = 0; i < A->size0; ++i)
      {
        CHECK_RETURN(cs_zentry(Ainv->matrix2->triplet, i, col_rhs, b[i]));
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
  return (int)info;

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
    case NS_CSC:
    {
      A->size0 = (int)A->matrix2->csc->m;
      A->size1 = (int)A->matrix2->csc->n;
      break;
    }
    case NS_CSR:
    {
      A->size0 = (int)A->matrix2->csc->m;
      A->size1 = (int)A->matrix2->csc->n;
      break;
    }
    case NS_TRIPLET:
    {
      A->size0 = (int)A->matrix2->triplet->m;
      A->size1 = (int)A->matrix2->triplet->n;
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

csi* NM_sparse_diag_indices(NumericsMatrix* M)
{
  NumericsSparseMatrix* A = M->matrix2;
  assert(A);
  if (A->diag_indx) return A->diag_indx;

  csi* indices = (csi*) malloc(M->size0 * sizeof(csi));
  A->diag_indx = indices;
  /* XXX hack --xhub  */
  if (A->origin == NS_TRIPLET) { NM_csc(M); A->origin = NS_CSC; }
  switch (A->origin)
  {
  case NS_CSC:
  {
    assert(A->csc);
    indices[0] = 0;
    CSparseMatrix* newMat = cs_spalloc(M->size0, M->size1, A->csc->p[M->size0]+M->size0, 1, 0);
    csi* Ai = A->csc->i;
    csi* Ap = A->csc->p;
    double* Ax = A->csc->x;
    csi* Ni = newMat->i;
    csi* Np = newMat->p;
    double* Nx = newMat->x;
    csi end = Ap[1];
    csi inc = 0;
    Np[0] = 0;
    if (Ai[0] == 0)
    {
      memcpy(Ni, Ai, end*sizeof(csi));
      Np[1] = Ap[1];
      memcpy(Nx, Ax, end*sizeof(double));
    }
    else
    {
      Ni[0] = 0;
      Np[1] = Ap[1] + 1;
      Nx[0] = 0.;
      memcpy(&Ni[1], Ai, end*sizeof(csi));
      memcpy(&Nx[1], Ax, end*sizeof(double));
      ++inc;
    }

    /* Could optimize further and copy everything using memcpy */
    for (size_t j = 1; j < (size_t)M->size0; ++j)
    {
      csi rem = 0;
      for (csi p = Ap[j]; (rem == 0) && (p < Ap[j+1]); ++p)
      {
        if (Ai[p] < (csi) j)
        {
          Ni[p+inc] = Ai[p];
          Nx[p+inc] = Ax[p];
        }
        else
        {
          if (Ai[p] > (csi) j)
          {
            Ni[p+inc] = j;
            Nx[p+inc] = 0.;
            indices[j] = p+inc;
            ++inc;
          }
          else
          {
            indices[j] = p+inc;
          }
          rem = p;
          Np[j] = Ap[j] + inc;
        }
        end = Ap[j+1] - rem;
        memcpy(&Ni[rem+inc], &Ai[rem], end*sizeof(csi));
        memcpy(&Nx[rem+inc], &Ax[rem], end*sizeof(double));
        assert(inc <= M->size0);
      }
    }
    Np[M->size0] = Ap[M->size0] + inc;
    NM_clearSparseStorage(M);
    A->origin = NS_CSC;
    A->csc = newMat;
    break;
  }
  case NS_TRIPLET:
  case NS_CSR:
  default:
    printf("NM_sparse_diag_indices :: unknown matrix origin %d", A->origin);
    exit(EXIT_FAILURE);
  }

  return indices;
}

void NM_csc_alloc(NumericsMatrix* A, csi nzmax)
{
  NM_sparse(A)->csc = cs_spalloc(A->size0, A->size1, nzmax, 1, 0);
}
void NM_csc_empty_alloc(NumericsMatrix* A, csi nzmax)
{
  NM_csc_alloc(A, nzmax);
  csi* Ap = NM_sparse(A)->csc->p;
  for (int i =0 ; i < A->size1+1 ; i++)  Ap[i]= 0;
  //csi* Ai = NM_sparse(A)->csc->i;
  //for (int i =0 ; i < nzmax ; i++)  Ai[i]= 0;
}


void NM_csr_alloc(NumericsMatrix* A, csi nzmax)
{
  NM_sparse(A)->csr = cs_spalloc(A->size1, A->size0, nzmax, 1, 0);
  NM_sparse(A)->csr->nz = -2;
  NM_sparse(A)->csr->m =  A->size0;
  NM_sparse(A)->csr->n = A->size1;
}

void NM_triplet_alloc(NumericsMatrix* A, csi nzmax)
{
  NM_sparse(A)->triplet = cs_spalloc(A->size0, A->size1, nzmax, 1, 1);
}

void NM_setSparseSolver(NumericsMatrix* A, unsigned solver_id)
{
  NM_linearSolverParams(A)->solver = (NumericsSparseLinearSolver)solver_id;
}

unsigned NM_sparse_origin(NumericsMatrix* M)
{
  assert(M);
  if (!M->matrix2) return -1;
  return M->matrix2->origin;
}

int NM_check(const NumericsMatrix* const A)
{
  int info = 0;
  if (!A->matrix2) return info;

  if (A->matrix2->csc) info = cs_check_csc(A->matrix2->csc);
  if (A->matrix2->csr) info = info ? info : cs_check_csc(A->matrix2->csr);
  if (A->matrix2->triplet) info = info ? info : cs_check_triplet(A->matrix2->triplet);
  return info;
}

CSparseMatrix* NM_sparse_get_origin(const NumericsMatrix* M)
{
  assert(M);
  if (!M->matrix2) { return NULL; }

  switch (M->matrix2->origin)
  {
  case NS_CSC:
    return M->matrix2->csc;
  case NS_TRIPLET:
    return M->matrix2->triplet;
  case NS_CSR:
    return M->matrix2->csr;
  default:
    numerics_error_nonfatal("NM_sparse_get_origin", "Unknown matrix origin %d", M->matrix2->origin);
    return NULL;
  }
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
    return NM_sparse_nnz(NM_sparse_get_origin(M));
  }
  default:
    numerics_error("NM_nnz", "Unsupported matrix type %d in %s", M->storageType);
    return SIZE_MAX;
  }
}
