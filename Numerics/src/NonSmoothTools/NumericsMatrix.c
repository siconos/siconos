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

#include "NumericsMatrix.h"
#include "SiconosLapack.h"
#include "misc.h"
#include "GlobalFrictionContact3D_AlartCurnier.h"
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


void freeNumericsMatrix(NumericsMatrix* m)
{
  assert(m && "freeNumericsMatrix, m == NULL");
  if (m->storageType == 0)
  {
    if (m->matrix0)
    {
      free(m->matrix0);
      m->matrix0 = NULL;
    }
  }
  else
  {
    if (m->matrix1)
    {
      freeSBM(m->matrix1);
      free(m->matrix1);
      m->matrix1 = NULL;
    }
    if (m->matrix2)
    {
      freeNumericsSparseMatrix(m->matrix2);
      free(m->matrix2);
      m->matrix2 = NULL;
    }
  }
  if (m->internalData)
  {
    if (m->internalData->iWork)
    {
      assert (m->internalData->iWorkSize > 0);
      free(m->internalData->iWork);
      m->internalData->iWork = NULL;
    }
    free(m->internalData);
    m->internalData = NULL;
  }
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
int cs_printInFile(const cs *A, int brief, FILE* file);
/* print a sparse matrix */
int cs_printInFile(const cs *A, int brief, FILE* file)
{
  ptrdiff_t m, n, nzmax, nz, p, j, *Ap, *Ai ;
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
    m->matrix0 = NULL;
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

  fillNumericsMatrix(m, storageType, (int)size0, (int)size1, data);

  return info;
}

void readInFileName(NumericsMatrix* const m, const char *filename)
{
  FILE* finput = fopen(filename, "r");
  printInFile(m, finput);
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
      printf("createNumericsMatrix :: storageType value %d not implemented yet !", mat->storageType);
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
      printf("createNumericsMatrix :: storageType value %d not implemented yet !", storageType);
      exit(EXIT_FAILURE);
  }

  fillNumericsMatrix(M, storageType, size0, size1, data);

  return M;
}


void fillNumericsMatrix(NumericsMatrix* M, int storageType, int size0, int size1, void* data)
{

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

void NM_clearSparseStorage(NumericsMatrix *A)
{
  NM_clearTriplet(A);
  NM_clearCSC(A);
  NM_clearCSCTranspose(A);
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

    if (A->matrix2->triplet)
    {
      A_ = A->matrix2->triplet;

      if (!B->matrix2->triplet)
      {
        B->matrix2->triplet = cs_spalloc(A_->m, A_->n, A_->nzmax, 0, 1);
      }

      B_ = B->matrix2->triplet;
    }
    else
    {
      assert (A->matrix2->csc);

      A_ = A->matrix2->csc;

      if (!B->matrix2->csc)
      {
        B->matrix2->csc = cs_spalloc(A_->m, A_->n, A_->nzmax, 0, 0);
      }

      B_ = B->matrix2->csc;
    }

    assert (A_);
    assert (B_);

    if (B_ ->nzmax < A_ ->nzmax)
    {
      B_->x = (double *) realloc(B_->x, A_->nzmax * sizeof(double));
      B_->i = (csi *) realloc(B_->i, A_->nzmax * sizeof(csi));
    }
    else if (!(B_->x))
    {
      B_->x = (double *) malloc(A_->nzmax * sizeof(double));
    }

    if (A_->nz >= 0)
    {
      /* triplet */
      B_->p = (csi *) realloc(B_->p, A_->nzmax * sizeof(csi));
    }
    else
    {
      if (B_->n < A_->n)
      {
        /* csc */
        B_-> p = (csi *) realloc(B_->p, (A_->n + 1) * sizeof(csi));
      }
    }

    B_->nzmax = A_->nzmax;
    B_->nz = A_->nz;
    B_->m = A_->m;
    B_->n = A_->n;

    memcpy(B_->x, A_->x, A_->nzmax * sizeof(double));
    memcpy(B_->i, A_->i, A_->nzmax * sizeof(csi));

    if (A_->nz >= 0)
    {
      memcpy(B_->p, A_->p, A_->nzmax * sizeof(csi));
    }
    else
    {
      memcpy(B_->p, A_->p, (A_->n + 1) * sizeof(csi));
    }


    /* invalidations */
    NM_clearDense(B);
    NM_clearSparseBlock(B);

    if (B_->nz >= 0)
    {
      NM_clearCSC(B);
      NM_clearCSCTranspose(B);
    }
    else
    {
      NM_clearTriplet(B);
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


/* NumericsMatrix : initialize triplet storage from sparse block storage */
CSparseMatrix* NM_triplet(NumericsMatrix* A)
{
  if (!NM_sparse(A)->triplet)
  {
    if (A->storageType == NM_SPARSE_BLOCK || A->storageType == NM_DENSE)
    {

      /* Invalidation of previously constructed csc storage. */
      /* If we want to avoid this -> rewrite cs_compress with reallocation. */
      NM_clearCSC(A);
      NM_clearCSCTranspose(A);

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
    }
  }
  assert (A->matrix2->triplet);

  return A->matrix2->triplet;
}

CSparseMatrix* NM_csc(NumericsMatrix *A)

{
  if(!NM_sparse(A)->csc)
  {
    assert(A->matrix2);
    A->matrix2->csc = cs_compress(NM_triplet(A)); /* triplet -> csc
                                                  * with allocation */

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

  default:
  {
    assert(A->storageType == NM_SPARSE);
    CHECK_RETURN(cs_aaxpy(alpha, NM_csc(A), x, beta, y));
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
        CHECK_RETURN(cs_aaxpy(alpha, NM_csc_trans(A), x, beta, y));
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
  switch(A->storageType)
  {
  case NM_DENSE:
  {
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, A->size0, B->size1, B->size1, alpha, A->matrix0, A->size0, B->matrix0, B->size0, beta, C->matrix0, A->size0);

    NM_clearSparseBlock(C);
    NM_clearSparseStorage(C);

    break;
  }
  case NM_SPARSE_BLOCK:
  {
    prodNumericsMatrixNumericsMatrix(alpha, A, B, beta, C);

    NM_clearDense(C);
    NM_clearSparseStorage(C);
    break;
  }
  case NM_SPARSE:
  {
    CSparseMatrix* result = cs_add(cs_multiply(NM_csc(A), NM_csc(B)),
                                   NM_csc(C), alpha, beta);

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

#ifdef WITH_MUMPS


MPI_Comm NM_MPI_com(NumericsMatrix* A)
{
  if (!A || (A && NM_linearSolverParams(A)->mpi_com == MPI_COMM_NULL))
  {
    int myid;
    int argc = 0;
    /* C99 requires that argv[argc] == NULL. With openmpi 1.8, we get a
     * segfault if this is not true */
    char *argv0 = NULL;
    char **argv = &argv0;
    CHECK_MPI(MPI_Init(&argc, &argv));
    CHECK_MPI(MPI_Comm_rank(MPI_COMM_WORLD, &myid));

    if (A)
    {
      NM_linearSolverParams(A)->mpi_com = MPI_COMM_WORLD;
      NM_linearSolverParams(A)->mpi_com_init = 1;
    }
  }

  if(A)
  {
    return NM_linearSolverParams(A)->mpi_com;
  }
  else
  {
    return MPI_COMM_WORLD;
  }
}

int* NM_MUMPS_irn(NumericsMatrix* A)
{

  if (NM_sparse(A)->triplet)
  {
    CSparseMatrix* triplet = NM_sparse(A)->triplet;
    csi nz = triplet->nz;

    int* iWork = NM_iWork(A, (int) (2*nz) + 1);

    for (int k=0; k<nz; ++k)
    {
      iWork [k + nz] = (int) (triplet->p [k]) + 1;
      iWork [k]      = (int) (triplet->i [k]) + 1;
    }

    iWork [2*nz] = (int) nz;
  }
  else
  {
    CSparseMatrix* csc = NM_sparse(A)->csc;
    csi nzmax = csc->nzmax ;

    int* iWork = NM_iWork(A, (int) (2*nzmax) + 1);

    csi n = csc->n ;
    csi nz = 0;
    csi* csci = csc->i ;
    csi* cscp = csc->p ;

    for (csi j=0; j<n; ++j)
    {
      for (csi p = cscp [j]; p < cscp [j+1]; ++p)
      {
        assert (csc->x [p] != 0.);
        nz++;
        iWork [p + nzmax] = (int) j;
        iWork [p]         = (int) csci [p];
      }
    }

    iWork [2*nzmax] = (int) nz;
  }

  return NM_iWork(A, 0);
}


int* NM_MUMPS_jcn(NumericsMatrix* A)
{
  if (NM_sparse(A)->triplet)
  {
    return NM_iWork(A, 0) + NM_sparse(A)->triplet->nz;
  }
  else
  {
    ptrdiff_t nzmax = NM_csc(A)->nzmax;
    return NM_iWork(A, 0) + nzmax;
  }
}


DMUMPS_STRUC_C* NM_MUMPS_id(NumericsMatrix* A)
{
  NumericsSparseLinearSolverParams* params = NM_linearSolverParams(A);

  if (!params->solver_data)
  {
    params->solver_data = malloc(sizeof(DMUMPS_STRUC_C));

    DMUMPS_STRUC_C* mumps_id = (DMUMPS_STRUC_C*) params->solver_data;

    // Initialize a MUMPS instance. Use MPI_COMM_WORLD.
    mumps_id->job = JOB_INIT;
    mumps_id->par = 1;
    mumps_id->sym = 0;

    if (NM_MPI_com(A) == MPI_COMM_WORLD)
    {
      mumps_id->comm_fortran = USE_COMM_WORLD;
    }
    else
    {
      mumps_id->comm_fortran = MPI_Comm_c2f(NM_MPI_com(A));
    }

    dmumps_c(mumps_id);

    if (verbose > 1)
    {
      mumps_id->ICNTL(1) = 6; // Error messages, standard output stream.
      mumps_id->ICNTL(2) = 6; // Diagnostics,    standard output stream.
      mumps_id->ICNTL(3) = 6; // Global infos,   standard output stream.

      mumps_id->ICNTL(4) = 4; // Errors, warnings and information on
                              // input, output parameters printed.

//      mumps_id->ICNTL(10) = 1; // One step of iterative refinment
      mumps_id->ICNTL(11) = 1; // Error analysis

    }
    else
    {
//      mumps_id->ICNTL(10) = 5;
//      mumps_id->CNTL(2) = 1e-14;

      mumps_id->ICNTL(1) = -1;
      mumps_id->ICNTL(2) = -1;
      mumps_id->ICNTL(3) = -1;
    }

    mumps_id->ICNTL(24) = 1; // Null pivot row detection see also CNTL(3) & CNTL(5)
    // ok for a cube on a plane & four contact points
    // computeAlartCurnierSTD != generated in this case...

    //mumps_id->CNTL(3) = ...;
    //mumps_id->CNTL(5) = ...;

  }
  DMUMPS_STRUC_C* mumps_id = (DMUMPS_STRUC_C*) params->solver_data;
  mumps_id->n = (int) NM_triplet(A)->n;
  mumps_id->irn = NM_MUMPS_irn(A);
  mumps_id->jcn = NM_MUMPS_jcn(A);

  int nz;
  if (NM_sparse(A)->triplet)
  {
    nz = (int) NM_sparse(A)->triplet->nz;
    mumps_id->nz = nz;
    mumps_id->a = NM_sparse(A)->triplet->x;
  }
  else
  {
    nz = NM_linearSolverParams(A)->iWork[2 * NM_csc(A)->nzmax];
    mumps_id->nz = nz;
    mumps_id->a = NM_sparse(A)->csc->x;
  }




  return (DMUMPS_STRUC_C*) params->solver_data;
}
#endif

int NM_gesv(NumericsMatrix* A, double *b)
{
  assert(A->size0 == A->size1);

  int info = 1;

  switch (A->storageType)
  {
  case NM_DENSE:
  {
    assert(A->matrix0);

    DGESV(A->size0, 1, A->matrix0, A->size0, NM_iWork(A, A->size0), b,
          A->size0, &info);
    break;
  }

  case NM_SPARSE_BLOCK: /* sparse block -> triplet -> csc */
  case NM_SPARSE:
  {
    switch (NM_linearSolverParams(A)->solver)
    {
    case NS_CS_LUSOL:
      info = !cs_lusol(1, NM_csc(A), b, DBL_EPSILON);
      break;

#ifdef WITH_MUMPS
    case NS_MUMPS:
    {
      /* the mumps instance is initialized (call with job=-1) */
      DMUMPS_STRUC_C* mumps_id = NM_MUMPS_id(A);

      mumps_id->rhs = b;
      mumps_id->job = 6;

      /* compute the solution */
      dmumps_c(mumps_id);

      /* clean the mumps instance */
      mumps_id->job = -2;
      dmumps_c(mumps_id);
      info = mumps_id->info[0];

      if (info > 0)
      {
        if (verbose > 0)
        {
          printf("NM_gesv: MUMPS fails : info(1)=%d, info(2)=%d\n", info, mumps_id->info[1]);
        }
      }
      if (verbose > 1)
      {
        printf("MUMPS : condition number %g\n", mumps_id->rinfog[9]);
        printf("MUMPS : component wise scaled residual %g\n", mumps_id->rinfog[6]);
        printf("MUMPS : \n");
      }

      /* Here we free mumps_id ...  */
      free(NM_linearSolverParams(A)->solver_data);
      NM_linearSolverParams(A)->solver_data = NULL;

      break;
    }
#endif
    default:
    {
      fprintf(stderr, "NM_gesv: unknown sparse linearsolver : %d\n", NM_linearSolverParams(A)->solver);
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
