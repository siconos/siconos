#include <stdio.h>
#include <stdlib.h>
#include "NumericsMatrix.h"
#include "LA.h"


void prod(int sizeX, int sizeY, const NumericsMatrix* const A, const double* const x, double* y, int init)
{
  if (A == NULL)
  {
    fprintf(stderr, "Numerics, NumericsMatrix, product matrix - vector prod(A,x,y) failed, A == NULL.\n");
    exit(EXIT_FAILURE);
  }
  /* check dim */
  if (A->size0 != sizeY || A->size1 != sizeX)
  {
    fprintf(stderr, "Numerics, NumericsMatrix, product matrix - vector prod(A,x,y) failed, unconsistent sizes.\n");
    exit(EXIT_FAILURE);
  }

  int storage = A->storageType;

  /* double* storage */
  if (storage == 0)
  {
    int incx = 1, incy = 1;
    double coef = 0.0; /* y = Ax */
    if (init == 0) /* y += Ax */
      coef = 1.0;
    DGEMV(LA_NOTRANS, sizeY, sizeX, 1.0, A->matrix0, sizeX, x, incx, coef, y, incy);
  }
  /* SparseBlock storage */
  else if (storage == 1)
    prodSBM(sizeY, A->matrix1, x, y, init);
  else
  {
    fprintf(stderr, "Numerics, NumericsMatrix, product matrix - vector prod(A,x,y) failed, unknown storage type for A.\n");
    exit(EXIT_FAILURE);
  }
}

void subRowProd(int sizeX, int sizeY, int currentRowNumber, const NumericsMatrix* A, const double* const x, double* y, int init)
{
  if (A == NULL)
  {
    fprintf(stderr, "Numerics, NumericsMatrix, product matrix - vector subRowProd(A,x,y) failed, A == NULL.\n");
    exit(EXIT_FAILURE);
  }
  /* check dim */
  if (A->size1 != sizeX || sizeY > A->size0)
  {
    fprintf(stderr, "Numerics, NumericsMatrix, product matrix - vector subRowProd(A,x,y) failed, unconsistent sizes.\n");
    exit(EXIT_FAILURE);
  }

  int storage = A->storageType;

  /* double* storage */
  if (storage == 0)
  {
    int incx = A->size0, incy = 1;
    double coef = 0.0; /* y = subAx */
    double* mat = A->matrix0;
    if (init == 0) /* y += subAx */
      coef = 1.0;
    for (int row = currentRowNumber; row < sizeY + currentRowNumber; row++)
      DDOT(sizeX, &mat[row], incx, x, incy);
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
