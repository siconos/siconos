/* Siconos-Numerics version 3.0.0, Copyright INRIA 2005-2008.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */
#include <stdio.h>
#include <stdlib.h>
#include "NumericsMatrix.h"
#include "LA.h"


void prod(int sizeX, int sizeY, double alpha, const NumericsMatrix* const A, const double* const x, double beta, double* y)
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
    DGEMV(LA_NOTRANS, sizeY, sizeX, alpha, A->matrix0, sizeX, x, incx, beta, y, incy);
  }
  /* SparseBlock storage */
  else if (storage == 1)
    prodSBM(sizeY, alpha, A->matrix1, x, beta, y);
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
  /* Checks inputs */
  if (A == NULL || x == NULL || y == NULL)
  {
    fprintf(stderr, "Numerics, NumericsMatrix, product matrix - vector rowProdNoDiag(A,x,y) failed, A or x or y == NULL.\n");
    exit(EXIT_FAILURE);
  }
  /* Check dim */
  if (A->size1 != sizeX || sizeY > A->size0)
  {
    fprintf(stderr, "Numerics, NumericsMatrix, product matrix - vector rowProdNoDiag(A,x,y) failed, unconsistent sizes.\n");
    exit(EXIT_FAILURE);
  }


  /* Checks storage type */
  int storage = A->storageType;

  /* double* storage */
  if (storage == 0)
  {
    fprintf(stderr, "Numerics, NumericsMatrix, product matrix - vector rowProdNoDiag(A,x,y) failed, not yet implemented for matrix not sparse.\n");
    exit(EXIT_FAILURE);
  }

  /*    exit(EXIT_FAILURE); */
  /*       if( currentRowNumber > A->size0 || (currentRowNumber+sizeY) > A->size0 ) */
  /*  { */
  /*    fprintf(stderr,"Numerics, NumericsMatrix, product matrix - vector rowProdNoDiag(A,x,y) failed, unconsistent sizes.\n"); */
  /*    exit(EXIT_FAILURE); */
  /*  } */
  /*       int incx = A->size0, incy = 1; */
  /*       double coef = 0.0; /\* y = subAx *\/ */
  /*       double* mat = A->matrix0; */
  /*       if(init==0) /\* y += subAx *\/ */
  /*  coef = 1.0; */
  /*       for(int row = currentRowNumber; row < sizeY+currentRowNumber; row++) */
  /*  DDOT(sizeX, &mat[row], incx, x, incy); */

  /* SparseBlock storage */
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
    if (m->matrix0 != NULL)
      free(m->matrix0);
    m->matrix0 = NULL;
  }
  else
    freeSBM(m->matrix1);
}


void display(const NumericsMatrix* const m)
{
  if (m == NULL)
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
  if (m == NULL)
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
