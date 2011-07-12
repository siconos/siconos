/* Siconos-Numerics, Copyright INRIA 2005-2010.
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
#include "NumericsMatrix.h"
#include "LA.h"
#include "misc.h"

void prodNumericsMatrix(int sizeX, int sizeY, double alpha, const NumericsMatrix* const A, const double* const x, double beta, double* y)
{

  assert(A);
  assert(x);
  assert(y);
  assert(A->size0 == sizeY);
  assert(A->size1 == sizeX);

  int storage = A->storageType;

  /* double* storage */
  if (storage == 0)
  {
    DGEMV(LA_NOTRANS, sizeY, sizeX, alpha, A->matrix0, sizeY, x, 1, beta, y, 1);
  }
  /* SparseBlock storage */
  else if (storage == 1)
  {
    prodSBM(sizeX, sizeY, alpha, A->matrix1, x, beta, y);
  }
  else
  {
    fprintf(stderr, "Numerics, NumericsMatrix, product matrix - vector prod(A,x,y) failed, unknown storage type for A.\n");
    exit(EXIT_FAILURE);
  }
}

void prodNumericsMatrixNumericsMatrix(double alpha, const NumericsMatrix* const A, const NumericsMatrix* const B, double beta,  NumericsMatrix* C)
{

  int astorage = A->storageType;
  int bstorage = B->storageType;
  int cstorage = C->storageType;
  assert(A);
  assert(B);
  assert(C);
  assert(A->size1 == B->size0);
  assert(C->size0 == A->size0);
  assert(C->size1 == B->size1);


  /* double* storage */
  if ((astorage == 0) & (bstorage == 0) & (cstorage == 0))
  {
    /*      DGEMV(LA_NOTRANS, sizeY, sizeX, alpha, A->matrix0, sizeY, x, 1, beta, y, 1); */
    DGEMM(LA_NOTRANS, LA_NOTRANS, A->size0, B->size1, A->size1, alpha, A->matrix0, A->size0, B->matrix0, B->size0, beta, C->matrix0, C->size0);
  }
  /* SparseBlock storage */
  else if ((astorage == 1) & (bstorage == 1) & (cstorage == 1))
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
        y[row] += DDOT(sizeX, &mat[currentRowNumber + row], incx, x, incy);
    }
    else
    {
      for (int row = 0; row < sizeY; row++)
        y[row] = DDOT(sizeX, &mat[currentRowNumber + row], incx, x, incy);
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
      y[i] += DDOT(A->size0 , MM + i , incx , x , incy);
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
  if (m->storageType == 0)
  {
    if (m->matrix0)
      free(m->matrix0);
    m->matrix0 = NULL;
  }
  else
  {
    freeSBM(m->matrix1);
    free(m->matrix1);
  }
}


void display(const NumericsMatrix* const m)
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
    printf("[");
    for (int i = 0; i < m->size1 * m->size0; i++)
    {
      printf("%lf ", m->matrix0[i]);
      if ((i + 1) % m->size1 == 0)
        printf("\n");
    }
    printf("]");
    printf("\n (warning: column-major) \n");
  }
  else if (storageType == 1)
    printSBM(m->matrix1);
}
void displayRawbyRaw(const NumericsMatrix* const m)
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

void printInFile(const NumericsMatrix* const m, FILE* file)
{
  if (! m)
  {
    fprintf(stderr, "Numerics, NumericsMatrix printInFile failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  int storageType = m->storageType;
  fprintf(file, "%d\n", m->storageType);
  fprintf(file, "%d\n", m->size0);
  fprintf(file, "%d\n", m->size1);

  if (storageType == 0)
  {
    fprintf(file, "%i\t%i\n", m->size0, m->size1);
    for (int i = 0; i < m->size1 * m->size0; i++)
    {
      fprintf(file, "%32.24e ", m->matrix0[i]);
      if ((i + 1) % m->size1 == 0)
        fprintf(file, "\n");
    }
  }
  else if (storageType == 1)
    printInFileSBM(m->matrix1, file);
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
      CHECK_IO(fscanf(file, "%le", &(m->matrix0[i])));
      if ((i + 1) % m->size1 == 0)
        IGNORE_IO(fscanf(file, "\n"));
    }


  }
  else if (storageType == 1)
  {
    m->matrix0 = NULL;
    readInFileSBM(m->matrix1, file);
  }
}


void newFromFile(NumericsMatrix* const m, FILE *file)
{
  if (! m)
  {
    fprintf(stderr, "Numerics, NumericsMatrix newFromFile failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  CHECK_IO(fscanf(file, "%d", &(m->storageType)));
  CHECK_IO(fscanf(file, "%d", &(m->size0)));
  CHECK_IO(fscanf(file, "%d", &(m->size1)));
  int storageType = m->storageType;

  if (storageType == 0)
  {
    CHECK_IO(fscanf(file, "%d\t%d\n", &(m->size0), &(m->size1)));

    m->matrix0 = (double *)malloc(m->size1 * m->size0 * sizeof(double));

    for (int i = 0; i < m->size1 * m->size0; i++)
    {
      CHECK_IO(fscanf(file, "%le", &(m->matrix0[i])));
      if ((i + 1) % m->size1 == 0)
        IGNORE_IO(fscanf(file, "\n"));
    }


  }
  else if (storageType == 1)
  {
    m->matrix0 = NULL;
    m->matrix1 = (SparseBlockStructuredMatrix*)malloc(sizeof(SparseBlockStructuredMatrix));
    newFromFileSBM(m->matrix1, file);
  }
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
};

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
