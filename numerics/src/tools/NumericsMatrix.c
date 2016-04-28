/* Siconos-Numerics, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#include "NumericsMatrix.h"
#include "NumericsMatrix_private.h"
#include "NM_conversions.h"
#include "SiconosLapack.h"
#include "NumericsOptions.h"
#include "misc.h"
#include "sanitizer.h"
//#define DEBUG_MESSAGES
#include "debug.h"

void prodNumericsMatrix(int sizeX, int sizeY, double alpha, NumericsMatrix* A, const double* const x, double beta, double* y)
{

  assert(A);
  assert(x);
  assert(y);
  assert(A->size0 == sizeY);
  assert(A->size1 == sizeX);

  int storage = A->storageType;

  /* double* storage */
  switch (storage)
  {
    case NM_DENSE:
      cblas_dgemv(CblasColMajor, CblasNoTrans, sizeY, sizeX, alpha, A->matrix0, sizeY, x, 1, beta, y, 1);
    break;
  /* SparseBlock storage */
    case NM_SPARSE_BLOCK:
      prodSBM(sizeX, sizeY, alpha, A->matrix1, x, beta, y);
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

void prodNumericsMatrixNumericsMatrix(double alpha, const NumericsMatrix* const A, const NumericsMatrix* const B, double beta,  NumericsMatrix* C)
{

  assert(A);
  assert(B);
  assert(C);
  int astorage = A->storageType;
  int bstorage = B->storageType;
  int cstorage = C->storageType;
  assert(A->size1 == B->size0);
  assert(C->size0 == A->size0);
  assert(C->size1 == B->size1);


  /* double* storage */
  if ((astorage == NM_DENSE) && (bstorage == NM_DENSE) && (cstorage == NM_DENSE))
  {
    /*      cblas_dgemv(CblasColMajor,CblasNoTrans, sizeY, sizeX, alpha, A->matrix0, sizeY, x, 1, beta, y, 1); */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, A->size0, B->size1, A->size1, alpha, A->matrix0, A->size0, B->matrix0, B->size0, beta, C->matrix0, C->size0);
  }
  /* SparseBlock storage */
  else if ((astorage == NM_SPARSE_BLOCK) && (bstorage == NM_SPARSE_BLOCK) && (cstorage == NM_SPARSE_BLOCK))
  {
    prodSBMSBM(alpha, A->matrix1, B->matrix1, beta, C->matrix1);
  }
  else
  {
    fprintf(stderr, "Numerics, NumericsMatrix, product matrix - matrix prod(A,B,C) not yet implemented.\n");
    exit(EXIT_FAILURE);
  }
}

void subRowProd(int sizeX, int sizeY, int currentRowNumber, const NumericsMatrix* A, const double* const x, double* y, int init)
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
    subRowProdSBM(sizeX, sizeY, currentRowNumber, A->matrix1, x, y, init);
  else
  {
    fprintf(stderr, "Numerics, NumericsMatrix, product matrix - vector subRowProd(A,x,y) failed, unknown storage type for A.\n");
    exit(EXIT_FAILURE);
  }

}

void rowProdNoDiag(int sizeX, int sizeY, int currentRowNumber, const NumericsMatrix* A, const double* const x, double* y, int init)
{

  assert(A);
  assert(x);
  assert(y);
  assert(A->size0 >= sizeY);
  assert(A->size1 == sizeX);

  /* Checks storage type */
  int storage = A->storageType;


  /* double* storage */
  if (storage == 0)
  {
    double * xSave = (double*) malloc(sizeY * sizeof(double));
    double * xx = (double *)x; /*because of const*/
    for (int i = 0; i < sizeY; i++)
    {
      xSave[i] = x[currentRowNumber + i];
      xx[currentRowNumber + i] = 0;
    }
    double * MM = A->matrix0 + currentRowNumber;
    int incx = A->size0;
    int incy = 1;
    if (init)
    {
      for (int i = 0; i < sizeY; i++)
        y[i] = 0;
    }
    for (int i = 0; i < sizeY; i++)
    {
      y[i] += cblas_ddot(A->size0 , MM + i , incx , x , incy);
    }
    for (int i = 0; i < sizeY; i++)
    {
      xx[currentRowNumber + i] = xSave[i];
    }
    free(xSave);

  }
  else if (storage == 1)
    rowProdNoDiagSBM(sizeX, sizeY, currentRowNumber, A->matrix1, x, y, init);
  else
  {
    fprintf(stderr, "Numerics, NumericsMatrix, product matrix - vector rowProdNoDiag(A,x,y) failed, unknown storage type for A.\n");
    exit(EXIT_FAILURE);
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

void freeNumericsMatrix(NumericsMatrix* m)
{
  assert(m && "freeNumericsMatrix, m == NULL");

  NM_clearDense(m);
  NM_clearSparseBlock(m);
  NM_clearSparse(m);

  NM_internalData_free(m);
}


void displayMat(double * m, int nRow, int nCol, int lDim)
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

void display(const NumericsMatrix* const m)
{
  if (! m)
  {
    fprintf(stderr, "Numerics, NumericsMatrix display failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  printf("\n ========== Numerics Matrix\n");
  printf("\n ========== storageType = %i\n", m->storageType);

  switch (m->storageType)
  {
  case NM_DENSE:
  {

    displayMat(m->matrix0, m->size0, m->size1, m->size0);
    break;
  }
  case NM_SPARSE_BLOCK:
  {
    assert(m->matrix1);
    printSBM(m->matrix1);
    break;
  }
  case NM_SPARSE:
  {
    assert(m->matrix2);
    if (m->matrix2->triplet)
    {
      cs_print(m->matrix2->triplet, 0);
    }
    else if (m->matrix2->csc)
    {
      cs_print(m->matrix2->csc, 0);
    }
    else if (m->matrix2->trans_csc)
    {
      cs_print(m->matrix2->trans_csc, 0);
    }
    else
    {
      fprintf(stderr, "display for sparse matrix: no matrix found!\n");
    }
    break;
  }
  default:
  {
    fprintf(stderr, "display for NumericsMatrix: matrix type %d not supported!\n", m->storageType);
  }
  }
}
void displayRowbyRow(const NumericsMatrix* const m)
{
  if (! m)
  {
    fprintf(stderr, "Numerics, NumericsMatrix display failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  int storageType = m->storageType;
  if (storageType == 0)
  {
    printf("\n ========== Numerics Matrix of dim %dX%d\n", m->size0, m->size1);
    for (int lin = 0; lin < m->size0; lin++)
    {
      for (int col = 0; col < m->size1; col++)
        printf("%lf ", m->matrix0[lin + col * m->size1]);
      printf("\n");
    }
  }
  else if (storageType == 1)
    printSBM(m->matrix1);
}

static int cs_printInFile(const cs *A, int brief, FILE* file)
{
  csi m, n, nzmax, nz, p, j, *Ap, *Ai ;
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
    fprintf(file,"%ld-by-%ld, nzmax: %ld nnz: %ld, 1-norm: %g\n", m, n, nzmax,
           Ap [n], cs_norm(A)) ;
    for(j = 0 ; j < n ; j++)
    {
      fprintf(file,"    col %ld : locations %ld to %ld\n", j, Ap [j], Ap [j+1]-1);
      for(p = Ap [j] ; p < Ap [j+1] ; p++)
      {
        fprintf(file,"      %ld : %g\n", Ai [p], Ax ? Ax [p] : 1) ;
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
    fprintf(file,"triplet: %ld-by-%ld, nzmax: %ld nnz: %ld\n", m, n, nzmax, nz) ;
    for(p = 0 ; p < nz ; p++)
    {
      fprintf(file,"    %ld %ld : %g\n", Ai [p], Ap [p], Ax ? Ax [p] : 1) ;
      if(brief && p > 20)
      {
        fprintf(file,"  ...\n") ;
        return (1) ;
      }
    }
  }
  return (1) ;
}


void printInFile(const NumericsMatrix* const m, FILE* file)
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
    printInFileSBM(m->matrix1, file);
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

void printInFileForScilab(const NumericsMatrix* const m, FILE* file)
{
  if (! m)
  {
    fprintf(stderr, "Numerics, NumericsMatrix printInFile failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  int storageType = m->storageType;
  fprintf(file, "storageType = %d ; \n", m->storageType);
  fprintf(file, "size0 = %d; \n", m->size0);
  fprintf(file, "size1 = %d; \n", m->size1);

  if (storageType == 0)
  {
    fprintf(file, "data= [");
    for (int i = 0; i < m->size0; i++)
    {
      fprintf(file, "[");
      for (int j = 0; j < m->size1; j++)
      {
        fprintf(file, "%32.24e,\t ", m->matrix0[i + j * m->size0]);
      }
      fprintf(file, "];\n");
    }
    fprintf(file, "]");
  }
  else if (storageType == 1)
  {
    printInFileSBMForScilab(m->matrix1, file);
    /*       fprintf(stderr,"Numerics, NumericsMatrix printInFileForScilab. Not yet implemented fo storageType = %i.\n", storageType); */
    /*       exit(EXIT_FAILURE); */

  }
}

void printInFileName(const NumericsMatrix* const m, const char *filename)
{
  FILE* foutput = fopen(filename, "w");
  printInFile(m, foutput);
  fclose(foutput);
}

void readInFile(NumericsMatrix* const m, FILE *file)
{
  if (! m)
  {
    fprintf(stderr, "Numerics, NumericsMatrix readInFile failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  CHECK_IO(fscanf(file, "%d", &(m->storageType)));
  CHECK_IO(fscanf(file, "%d", &(m->size0)));
  CHECK_IO(fscanf(file, "%d", &(m->size1)));
  int storageType = m->storageType;

  if (storageType == 0)
  {
    CHECK_IO(fscanf(file, "%d\t%d\n", &(m->size0), &(m->size1)));

    for (int i = 0; i < m->size1 * m->size0; i++)
    {
      CHECK_IO(fscanf(file, "%le ", &(m->matrix0[i])));
    }


  }
  else if (storageType == 1)
  {
    readInFileSBM(m->matrix1, file);
  }
}


int newFromFile(NumericsMatrix* const m, FILE *file)
{
  if (! m)
  {
    fprintf(stderr, "Numerics, NumericsMatrix newFromFile failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }

  if (m->storageType >= 0)
  {
    freeNumericsMatrix(m);
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
    data = newSBM();
    newFromFileSBM((SparseBlockStructuredMatrix*)data, file);
  }
  else
  {
    fprintf(stderr, "newFromFile :: storageType %d not supported yet", storageType);
  }

  fillNumericsMatrix(m, storageType, (int)size0, (int)size1, data);

  return info;
}

void readInFileName(NumericsMatrix* const m, const char *filename)
{
  FILE* finput = fopen(filename, "r");
  newFromFile(m, finput);
  fclose(finput);
}

void readInFileForScilab(NumericsMatrix* const M, FILE *file)
{
  fprintf(stderr, "Numerics, NumericsMatrix,readInFileForScilab");
  exit(EXIT_FAILURE);
}

void getDiagonalBlock(NumericsMatrix* m, int numBlockRow, int numRow, int size, double ** Bout)
{
  int storageType = m->storageType;
  if (storageType == 0)
  {
    double * MM = m->matrix0;
    double * elem = 0;
    /* The part of MM which corresponds to the current block is copied into MLocal */
    for (int i = 0; i < size; i++)
    {
      elem = MM + numRow + (numRow + i) * (m->size0);
      for (int j = 0; j < size; j++)
      {
        (*Bout)[j + i * size] = *elem;
        elem++;
      }
    }
  }
  else if (storageType == 1)
  {
    int diagPos = getDiagonalBlockPos(m->matrix1, numBlockRow);
    (*Bout) = m->matrix1->block[diagPos];

  }
}

NumericsMatrix* createNumericsMatrixFromData(int storageType, int size0, int size1, void* data)
{
  NumericsMatrix* M = newNumericsMatrix();

  fillNumericsMatrix(M, storageType, size0, size1, data);

  return M;
}

NumericsMatrix* duplicateNumericsMatrix(NumericsMatrix* mat)
{
  NumericsMatrix* M = newNumericsMatrix();

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
      data = malloc(sizeof(SparseBlockStructuredMatrix));
      break;
    case NM_SPARSE:
      data = malloc(sizeof(CSparseMatrix));
      break;
    default:
      fprintf(stderr, "duplicateNumericsMatrix :: storageType value %d not implemented yet !", mat->storageType);
      exit(EXIT_FAILURE);
  }

  fillNumericsMatrix(M, mat->storageType, size0, size1, data);

  return M;
}

NumericsMatrix* newNumericsMatrix(void)
{
  NumericsMatrix* M = (NumericsMatrix*) malloc(sizeof(NumericsMatrix));
  M->storageType = -1;
  M->size0 = 0;
  M->size1 = 0;
  NM_null(M);

  return M;
}

NumericsMatrix* createNumericsMatrix(int storageType, int size0, int size1)
{
  NumericsMatrix* M = newNumericsMatrix();

  void* data;

  switch (storageType)
  {
    case NM_DENSE:
      data = malloc(size0*size1*sizeof(double));
      break;
    case NM_SPARSE_BLOCK:
      data = newSBM();
      break;
    case NM_SPARSE:
      data = newNumericsSparseMatrix();
      break;
    default:
      fprintf(stderr, "createNumericsMatrix :: storageType value %d not implemented yet !", storageType);
      exit(EXIT_FAILURE);
  }

  fillNumericsMatrix(M, storageType, size0, size1, data);

  return M;
}


void fillNumericsMatrix(NumericsMatrix* M, int storageType, int size0, int size1, void* data)
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
          if (M->matrix2->origin == NS_UNKNOWN)
          {
            if (M->matrix2->triplet) { M->matrix2->origin = NS_TRIPLET; }
            else if (M->matrix2->csc) { M->matrix2->origin = NS_CSC; }
            else if (M->matrix2->csr) { M->matrix2->origin = NS_CSR; }
          }
        }
        break;

      default:
        printf("fillNumericsMatrix :: storageType value %d not implemented yet !", storageType);
        exit(EXIT_FAILURE);
    }
  }
}

NumericsMatrix* newSparseNumericsMatrix(int size0, int size1, SparseBlockStructuredMatrix* m1)
{
  return createNumericsMatrixFromData(NM_SPARSE_BLOCK, size0, size1, (void*)m1);
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
    freeSBM(A->matrix1);
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
  NM_clearTriplet(A);
  NM_clearCSC(A);
  NM_clearCSCTranspose(A);
  NM_clearCSR(A);
}

static inline void NM_dense_to_sparse(const NumericsMatrix* const A, NumericsMatrix* B)
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
      B->matrix1 = newSBM();
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

    /* cf copySBM */
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
    prodSBM(A->size1, A->size0, alpha, A->matrix1, x, beta, y);

     break;
  }

  case NM_SPARSE:
  default:
  {
    assert(A->storageType == NM_SPARSE);
    // if possible use the much simpler version provided by CSparse
    // Also at the time of writing, cs_aaxpy is bugged --xhub
    if (fabs(alpha - 1.) < 100*DBL_EPSILON && fabs(beta - 1.) < 100*DBL_EPSILON)
    {
      CHECK_RETURN(cs_gaxpy(NM_csc(A), x, y));
    }
    else
    {
      CHECK_RETURN(cs_aaxpy(alpha, NM_csc(A), x, beta, y));
    }
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
        if (fabs(alpha - 1.) < 100*DBL_EPSILON && fabs(beta - 1.) < 100*DBL_EPSILON)
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
    prodNumericsMatrixNumericsMatrix(alpha, A, B, beta, C);

    NM_clearDense(C);
    NM_clearSparseStorage(C);
    break;
  }
  case NM_SPARSE:
  {
    /* We need to free the allocated data here, hence we have to save the
     * matrix pointer. Otherwise, we have a memleak */
    CSparseMatrix* tmp_matrix = cs_multiply(NM_csc(A), NM_csc(B));
    assert(tmp_matrix && "NM_gemm :: cs_multiply failed");
    CSparseMatrix* result = cs_add(tmp_matrix,
                                   NM_csc(C), alpha, beta);
    assert(result && "NM_gemm :: cs_add failed");

    cs_spfree(tmp_matrix);
    NM_clearDense(C);
    NM_clearSparseBlock(C);
    NM_clearSparseStorage(C);

    NM_sparse(C)->csc = result;
    C->size0 = (int)C->matrix2->csc->m;
    C->size1 = (int)C->matrix2->csc->n;
    break;
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
  }
  return A->internalData;
}

int* NM_iWork(NumericsMatrix* A, int size)
{
  if (!NM_internalData(A)->iWork)
  {
    assert(A->internalData);

    assert(A->internalData->iWorkSize == 0);
    A->internalData->iWork = (int *) malloc(size * sizeof(int));
    A->internalData->iWorkSize = size;
  }
  else
  {
    assert(A->internalData);

    if (size > A->internalData->iWorkSize)
    {
      A->internalData->iWork = (int *) realloc(A->internalData->iWork, size * sizeof(int));
      A->internalData->iWorkSize = size;
    }
  }

  assert(A->internalData->iWork);
  assert(A->internalData->iWorkSize >= size);

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

    if (size > A->internalData->dWorkSize)
    {
      A->internalData->dWork = (double *) realloc(A->internalData->dWork, size * sizeof(double));
      A->internalData->dWorkSize = size;
    }
  }

  assert(A->internalData->dWork);
  assert(A->internalData->dWorkSize >= size);

  return A->internalData->dWork;
}

int NM_gesv_expert(NumericsMatrix* A, double *b, bool keep)
{
  assert(A->size0 == A->size1);

  int info = 1;

  switch (A->storageType)
  {
  case NM_DENSE:
  {
    assert(A->matrix0);

    if (keep)
    {
      double* wkspace = NM_dWork(A, A->size0*A->size1);
      int* ipiv = NM_iWork(A, A->size0);
      if (!NM_internalData(A)->isLUfactorized)
      {
        cblas_dcopy_msan(A->size0*A->size1, A->matrix0, 1, wkspace, 1);
        DGETRF(A->size0, A->size1, wkspace, A->size0, ipiv, &info);
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

      DGETRS(LA_NOTRANS, A->size0, 1, wkspace, A->size0, ipiv, b, A->size0, &info);
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
      DGESV(A->size0, 1, A->matrix0, A->size0, NM_iWork(A, A->size0), b,
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
      if (keep)
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

      if (!keep || mumps_id->job == -1)
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
      if (!keep)
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

      if (!keep)
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
    default:
    {
      fprintf(stderr, "NM_gesv: unknown sparse linearsolver : %d\n", p->solver);
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
  return info;
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

void NM_csc_alloc(NumericsMatrix* A, csi nzmax)
{
  NM_sparse(A)->csc = cs_spalloc(A->size0, A->size1, nzmax, 1, 0);
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
    fprintf(stderr, "NM_copy :: error unknown type %d for CSparse matrix\n", m->nz);
    exit(EXIT_FAILURE);
  }

  return out;
}

