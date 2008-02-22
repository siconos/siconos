/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2007.
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

/*
  Tests functions for NumericsMatrix structure

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NumericsMatrix.h"
#include "LA.h"
#include <math.h>

void test_BuildNumericsMatrix(NumericsMatrix** MM1, NumericsMatrix** MM2)
{
  /* Build two equal Numerics Matrices, one with double* storage (MM1), the other with sparse storage (MM2)*/
  NumericsMatrix * M1 = *MM1;
  NumericsMatrix * M2 = *MM2;
  int n = 8;
  /* Double * storage (column-major) */
  double m0[] = {1, 2, 0, 5, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 2, 2, 4, 0, -1, 6, 0, 0, 1, 2, 3, 4, 0, 0, 1, 0, 0, 0, -1, 1, 0, 6, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 0, 0, 0, 0, 5, 2, 2, 2};
  M1->storageType = 0;
  M1->size0 = n;
  M1->size1 = n;
  M1->matrix0 = malloc(n * n * sizeof(double));
  /* Note: M1->matrix0 = &m0[0] results in strange behavior ... */
  int i;
  for (i = 0; i < n * n; i++)
    M1->matrix0[i] = m0[i];
  M1->matrix1 = NULL;
  /* Build a NumericsMatrix with sparse-block storage */
  M2->storageType = 1;
  M2->size0 = n;
  M2->size1 = n;
  M2->matrix0 = NULL;
  M2->matrix1 = malloc(1 * sizeof(*(M2->matrix1)));
  SparseBlockStructuredMatrix * SBM = M2->matrix1;
  SBM->nbblocks = 6;
  SBM->size = 3;
  int sizes[] = {4, 6, 8};
  SBM->blocksize = malloc(SBM->size * sizeof(int));
  for (i = 0; i < SBM->size; i++)
    SBM->blocksize[i] = sizes[i];
  int i1[] = {0, 0, 1, 1, 2, 2};
  SBM->RowIndex = malloc(SBM->nbblocks * sizeof(int));
  for (i = 0; i < SBM->nbblocks; i++)
    SBM->RowIndex[i] = i1[i];
  SBM->ColumnIndex = malloc(SBM->nbblocks * sizeof(int));
  int i2[] = {0, 1, 1, 2, 0, 2};
  for (i = 0; i < SBM->nbblocks; i++)
    SBM->ColumnIndex[i] = i2[i];
  SBM->block = malloc(SBM->nbblocks * sizeof(* (SBM->block)));
  double block0[] = {1, 2 , 0 , 5 , 2 , 1 , 0 , 0 , 0 , 0 , 1 , -1, 4, 0 , -1, 6};
  double block1[] = {3, 4, 0, 0, -1, 1, 0, 6};
  double block2[] = {1, 0, 0, 2};
  double block3[] = {0, 0, 5, 2};
  double block4[] = {0, 0, 0, 0, 2, 2, 1, 2};
  double block5[] = {2, -1, 2, 2};
  SBM->block[0] = malloc(16 * sizeof(double));
  SBM->block[1] = malloc(8 * sizeof(double));
  SBM->block[2] = malloc(4 * sizeof(double));
  SBM->block[3] = malloc(4 * sizeof(double));
  SBM->block[4] = malloc(8 * sizeof(double));
  SBM->block[5] = malloc(4 * sizeof(double));
  for (i = 0; i < 16; i++)
    SBM->block[0][i] = block0[i];
  for (i = 0; i < 8; i++)
    SBM->block[1][i] = block1[i];
  for (i = 0; i < 4; i++)
    SBM->block[2][i] = block2[i];
  for (i = 0; i < 4; i++)
    SBM->block[3][i] = block3[i];
  for (i = 0; i < 8; i++)
    SBM->block[4][i] = block4[i];
  for (i = 0; i < 4; i++)
    SBM->block[5][i] = block5[i];
}

int test_prod(NumericsMatrix* M1, NumericsMatrix* M2)
{

  printf("== Numerics tests: prod(NumericsMatrix,vector) == \n");
  int i , n = M1->size1;
  double * x = malloc(n * sizeof(double));
  double alpha = 2.3, beta = 1.9;
  double yref[n];
  double * y = malloc(n * sizeof(double));
  for (i = 0; i < n; i++)
  {
    x[i] = i + 1.0;
    yref[i] = 0.1 * i;
    y[i] = yref[i];
  }

  int incx = 1, incy = 1;
  DGEMV(LA_NOTRANS, n, n, alpha, M1->matrix0, n, x, incx, beta, yref, incy);

  prod(n, n, alpha, M1, x, beta, y);
  double tol = 1e-12;
  int info = 0;
  for (i = 0; i < n; i++)
  {
    if (fabs(y[i] - yref[i]) > tol) info = 1;
    //    printf("%lf\n", fabs(y[i]-yref[i]));
  }
  if (info == 0)
    printf("Step 0 ( y = alpha*A*x + beta*y, double* storage) ok ...\n");
  else
    printf("Step 0 ( y = alpha*A*x + beta*y, double* storage) failed ...\n");

  /* Sparse ... */
  for (i = 0; i < n; i++)
  {
    y[i] = 0.1 * i;
  }
  prod(n, n, alpha, M2, x, beta, y);
  for (i = 0; i < n; i++)
  {
    if (fabs(y[i] - yref[i]) > tol) info = 1;
    /*      printf("%lf\n", fabs(y[i]-yref[i])); */
  }
  if (info == 0)
    printf("Step 2 ( y = alpha*A*x + beta*y, sparse storage) ok ...\n");
  else
    printf("Step 2 ( y = alpha*A*x + beta*y,  sparsestorage) failed ...\n");


  free(x);
  free(y);
  printf("== End of test prod(NumericsMatrix,vector), result = %d\n", info);

  return info;
}

int test_subRowprod(NumericsMatrix* M1, NumericsMatrix* M2)
{

  printf("== Numerics tests: subRowProd(NumericsMatrix,vector) == \n");
  int i , n = M1->size1;
  double * x = malloc(n * sizeof(double));

  for (i = 0; i < n; i++)
  {
    x[i] = i + 1;
  }

  int min = 2;
  int max = 6;
  int sizeY = max - min;
  /* Computes yRef = subA*x, subA = A limited to row min to max*/
  double * y = malloc(sizeY * sizeof(double));
  double yref[sizeY];
  int incx = n, incy = 1;
  for (i = 0; i < sizeY; i++)
    yref[i] = DDOT(n, &(M1->matrix0[min + i]), incx, x, incy);

  subRowProd(n, sizeY, min, M1, x, y, 1);
  double tol = 1e-12;
  int info = 0;
  for (i = 0; i < sizeY; i++)
  {
    if (fabs(y[i] - yref[i]) > tol) info = 1;
    /*       printf("%lf\n", fabs(y[i]-yref[i])); */
  }
  if (info == 0)
    printf("Step 0 ( y = subA*x, double* storage) ok ...\n");
  else
    printf("Step 0 ( y = subA*x, double* storage) failed ...\n");

  /* += */
  subRowProd(n, sizeY, min, M1, x, y, 0);
  for (i = 0; i < sizeY; i++)
  {
    if (fabs(y[i] - 2 * yref[i]) > tol) info = 1;
    /*       printf("%lf\n", fabs(y[i]-2*yref[i]));  */
  }
  if (info == 0)
    printf("Step 1 ( y += subA*x, double* storage) ok ...\n");
  else
    printf("Step 1 ( y += subA*x, double* storage) failed ...\n");

  free(y);
  sizeY = 2;
  int pos = 1; // pos of the required row of blocks
  y = malloc(sizeY * sizeof(double));
  for (i = 0; i < sizeY; i++)
    yref[i] = DDOT(n, &(M1->matrix0[4 + i]), incx, x, incy);
  /* Sparse ... */
  subRowProd(n, sizeY, pos, M2, x, y, 1);
  for (i = 0; i < sizeY; i++)
  {
    if (fabs(y[i] - yref[i]) > tol) info = 1;
    //  printf("%lf\n", fabs(y[i]-yref[i]));
  }
  for (i = 0; i < sizeY; i++)
    yref[i] = DDOT(n, &(M1->matrix0[6 + i]), incx, x, incy);
  subRowProd(n, sizeY, pos + 1, M2, x, y, 1);
  for (i = 0; i < sizeY; i++)
  {
    if (fabs(y[i] - yref[i]) > tol) info = 1;
    //  printf("%lf\n", fabs(y[i]-yref[i]));
  }


  if (info == 0)
    printf("Step 2 ( y = subA*x, sparse storage) ok ...\n");
  else
    printf("Step 2 ( y = subA*x,  sparse storage) failed ...\n");

  /* Sparse, += ... */
  subRowProd(n, sizeY, pos + 1, M2, x, y, 0);
  for (i = 0; i < sizeY; i++)
  {
    if (fabs(y[i] - 2 * yref[i]) > tol) info = 1;
    /*       printf("%lf\n", fabs(y[i]-yref[i])); */
  }
  if (info == 0)
    printf("Step 3 ( y += subA*x, sparse storage) ok ...\n");
  else
    printf("Step 3 ( y += subA*x,  sparse storage) failed ...\n");


  free(x);
  free(y);
  printf("== End of test subRowProd(NumericsMatrix,vector), result = %d\n", info);

  return info;
}

int test_rowProdNoDiag(NumericsMatrix* M1, NumericsMatrix* M2)
{

  printf("== Numerics tests: rowProdNoDiag(NumericsMatrix,vector) == \n");
  int i , n = M1->size1;
  double * x = malloc(n * sizeof(double));

  for (i = 0; i < n; i++)
  {
    x[i] = i + 1;
  }

  int min = 2;
  int max = 6;
  int sizeY = max - min;
  /* Computes yRef = subA*x, subA = A limited to row min to max*/
  double * y = malloc(sizeY * sizeof(double));
  double yref[sizeY];
  //  int incx = n, incy =1;
  double tol = 1e-12;
  int info = 0;
  /*   for(i=0;i<sizeY;i++) */
  /*     yref[i]=DDOT(n, &(M1->matrix0[min+i]), incx, x, incy); */

  /*   rowProdNoDiag(n,sizeY,min,M1,x,y,1); */
  /*   for(i = 0; i< sizeY; i++) */
  /*     { */
  /*       if( fabs(y[i]-yref[i])>tol) info = 1;  */
  /* /\*       printf("%lf\n", fabs(y[i]-yref[i])); *\/ */
  /*     } */
  /*   if(info ==0) */
  /*     printf("Step 0 ( y = subA*x, double* storage) ok ...\n"); */
  /*   else */
  /*     printf("Step 0 ( y = subA*x, double* storage) failed ...\n"); */

  /*   /\* += *\/ */
  /*   rowProdNoDiag(n,sizeY,min,M1,x,y,0); */
  /*   for(i = 0; i< sizeY; i++) */
  /*     { */
  /*       if( fabs(y[i]-2*yref[i])>tol) info = 1;  */
  /*       /\*       printf("%lf\n", fabs(y[i]-2*yref[i]));  *\/ */
  /*     } */
  /*   if(info ==0) */
  /*     printf("Step 1 ( y += subA*x, double* storage) ok ...\n"); */
  /*   else */
  /*     printf("Step 1 ( y += subA*x, double* storage) failed ...\n"); */

  free(y);
  sizeY = 2;
  int pos = 1; // pos of the required row of blocks
  y = malloc(sizeY * sizeof(double));
  yref[0] = 40;
  yref[1] = 16;

  /* Sparse ... */
  rowProdNoDiag(n, sizeY, pos, M2, x, y, 1);
  for (i = 0; i < sizeY; i++)
  {
    if (fabs(y[i] - yref[i]) > tol) info = 1;
    //  printf("%lf\n", fabs(y[i]-yref[i]));
  }
  rowProdNoDiag(n, sizeY, pos + 1, M2, x, y, 1);
  yref[0] = 10;
  yref[1] = 14;
  for (i = 0; i < sizeY; i++)
  {
    if (fabs(y[i] - yref[i]) > tol) info = 1;
    //  printf("%lf\n", fabs(y[i]-yref[i]));
  }

  if (info == 0)
    printf("Step 2 ( y = subA*x, sparse storage) ok ...\n");
  else
    printf("Step 2 ( y = subA*x,  sparse storage) failed ...\n");

  /* Sparse, += ... */
  rowProdNoDiag(n, sizeY, pos + 1, M2, x, y, 0);
  for (i = 0; i < sizeY; i++)
  {
    if (fabs(y[i] - 2 * yref[i]) > tol) info = 1;
    /*       printf("%lf\n", fabs(y[i]-yref[i])); */
  }
  if (info == 0)
    printf("Step 3 ( y += subA*x, sparse storage) ok ...\n");
  else
    printf("Step 3 ( y += subA*x,  sparse storage) failed ...\n");


  free(x);
  free(y);
  printf("== End of test rowProdNoDiag(NumericsMatrix,vector), result = %d\n", info);

  return info;
}

int main(void)
{

  printf("========= Starts Numerics tests for NumericsMatrix ========= \n");
  NumericsMatrix * M1 = malloc(sizeof(*M1));
  NumericsMatrix * M2 = malloc(sizeof(*M2));

  test_BuildNumericsMatrix(&M1, &M2);
  printf("Construction ok ...\n");
  int info = test_prod(M1, M2);
  printf("End of Prod ...\n");
  info = test_subRowprod(M1, M2);
  printf("End of Sub-Prod ...\n");
  info = test_rowProdNoDiag(M1, M2);
  printf("End of Sub-Prod no diag ...\n");

  /* free memory */
  M1->matrix0 = NULL;
  free(M1->matrix0);
  freeSBM(M2->matrix1);
  free(M1);
  free(M2);
  printf("========= End Numerics tests for NumericsMatrix ========= \n");
  return info;
}

