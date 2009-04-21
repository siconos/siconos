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

/*
  Tests functions for NumericsMatrix structure

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NumericsMatrix.h"
#include "LA.h"
#include <math.h>

int test_BuildNumericsMatrix(NumericsMatrix** MM)
{

  NumericsMatrix * M1 = MM[0];
  NumericsMatrix * M2 = MM[1];
  NumericsMatrix * M3 = MM[2];
  NumericsMatrix * M4 = MM[3];


  int info = 0;
  /* Build two equal Numerics Matrices, one with double* storage (MM1), the other with sparse storage (MM2)*/
  int n = 8;
  /* Double * storage (column-major) */
  double m0[64] = {1, 2, 0, 5, 0, 0, 0, 0,
                   2, 1, 0, 0, 0, 0, 0, 0,
                   0, 0, 1, -1, 0, 0, 2, 2,
                   4, 0, -1, 6, 0, 0, 1, 2,
                   3, 4, 0, 0, 1, 0, 0, 0,
                   -1, 1, 0, 6, 0, 2, 0, 0,
                   0, 0, 0, 0, 0, 0, 2, -1,
                   0, 0, 0, 0, 5, 2, 2, 2
                  };
  M1->storageType = 0;
  M1->size0 = n;
  M1->size1 = n;
  M1->matrix0 = malloc(n * n * sizeof(double));
  /* Note: M1->matrix0 = &m0[0] results in strange behavior ... */
  int i;
  for (i = 0; i < n * n; i++)
    M1->matrix0[i] = m0[i];
  M1->matrix1 = NULL;

  int nn = 2;
  /* Double * storage (column-major) */
  double m00[] = {1, 2, 0, 5, 0, 0, 0, 0,
                  2, 1, 0, 0, 0, 0, 0, 0,
                  0, 0, 1, -1, 0, 0, 2, 2,
                  4, 0, -1, 6, 0, 0, 1, 2
                 };
  M3->storageType = 0;
  M3->size0 = n;
  M3->size1 = nn;
  M3->matrix0 = malloc(nn * n * sizeof(double));
  for (i = 0; i < n * n; i++)
    M3->matrix0[i] = m00[i];
  M3->matrix1 = NULL;


  /* Build a NumericsMatrix with sparse-block storage */
  M2->storageType = 1;
  M2->size0 = n;
  M2->size1 = n;
  M2->matrix0 = NULL;


  SparseBlockStructuredMatrix * SBM = (SparseBlockStructuredMatrix *)malloc(sizeof(SparseBlockStructuredMatrix));
  M2->matrix1 = SBM;
  SBM->nbblocks = 6;
  SBM->blocknumber0 = 3;

  SBM->blocksize0 = (int*)malloc(3 * sizeof(int));
  SBM->blocksize0[0] = 4;
  SBM->blocksize0[1] = 6;
  SBM->blocksize0[2] = 8;

  SBM->filled1 = 4;
  SBM->filled2 = SBM->nbblocks;

  SBM->index1_data = (size_t*)malloc((SBM->filled1) * sizeof(size_t));
  SBM->index1_data[0] = 0;
  SBM->index1_data[1] = 2;
  SBM->index1_data[2] = 4;
  SBM->index1_data[3] = 6;

  SBM->index2_data = (size_t*)malloc((SBM->filled2) * sizeof(size_t));
  SBM->index2_data[0] =  0;
  SBM->index2_data[1] =  1;
  SBM->index2_data[2] =  1;
  SBM->index2_data[3] =  2;
  SBM->index2_data[4] =  0;
  SBM->index2_data[5] =  2;

  SBM->block = malloc(SBM->nbblocks * sizeof(*(SBM->block)));
  double block0[16] = {1, 2 , 0 , 5 , 2 , 1 , 0 , 0 , 0 , 0 , 1 , -1, 4, 0 , -1, 6};
  double block1[8] = {3, 4, 0, 0, -1, 1, 0, 6};
  double block2[4] = {1, 0, 0, 2};
  double block3[4] = {0, 0, 5, 2};
  double block4[8] = {0, 0, 0, 0, 2, 2, 1, 2};
  double block5[4] = {2, -1, 2, 2};
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


  /* Build a NumericsMatrix with sparse-block storage */
  M4->storageType = 1;
  M4->size0 = n;
  M4->size1 = 4;
  M4->matrix0 = NULL;


  SparseBlockStructuredMatrix * SBM2 = (SparseBlockStructuredMatrix *)malloc(sizeof(SparseBlockStructuredMatrix));
  M4->matrix1 = SBM2;
  SBM2->nbblocks = 2;
  SBM2->blocknumber0 = 3;

  SBM2->blocksize0 = (int*)malloc(SBM2->blocknumber0 * sizeof(int));
  SBM2->blocksize0[0] = 4;
  SBM2->blocksize0[1] = 6;
  SBM2->blocksize0[2] = 8;

  SBM2->blocknumber1 = 1;
  SBM2->blocksize1 = (int*)malloc(SBM2->blocknumber1 * sizeof(int));
  SBM2->blocksize1[0] = 4;

  SBM2->filled1 = 4;
  SBM2->filled2 = SBM2->nbblocks;

  SBM2->index1_data = (size_t*)malloc((SBM2->filled1) * sizeof(size_t));
  SBM2->index1_data[0] = 0;
  SBM2->index1_data[1] = 1;
  SBM2->index1_data[2] = 1;
  SBM2->index1_data[3] = 2;

  SBM2->index2_data = (size_t*)malloc((SBM2->filled2) * sizeof(size_t));
  SBM2->index2_data[0] =  0;
  SBM2->index2_data[1] =  0;

  SBM2->block = malloc(SBM2->nbblocks * sizeof(*(SBM2->block)));
  double block00[16] = {1, 2 , 0 , 5 , 2 , 1 , 0 , 0 , 0 , 0 , 1 , -1, 4, 0 , -1, 6};
  double block40[8] = {0, 0, 0, 0, 2, 2, 1, 2};
  SBM2->block[0] = malloc(16 * sizeof(double));
  SBM2->block[1] = malloc(8 * sizeof(double));
  for (i = 0; i < 16; i++)
    SBM2->block[0][i] = block00[i];
  for (i = 0; i < 8; i++)
    SBM2->block[1][i] = block40[i];
  return info;
}

int test_prodNumericsMatrix(NumericsMatrix* M1, NumericsMatrix* M2)
{

  printf("== Numerics tests: prodNumericsMatrix(NumericsMatrix,vector) == \n");
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

  prodNumericsMatrix(n, n, alpha, M1, x, beta, y);
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
  prodNumericsMatrix(n, n, alpha, M2, x, beta, y);
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
  printf("== End of test prodNumericsMatrix(NumericsMatrix,vector), result = %d\n", info);

  return info;
}
int test_prodNumericsMatrixNumericsMatrix(NumericsMatrix** MM)
{


  NumericsMatrix * M1 = MM[0];
  NumericsMatrix * M2 = MM[0];
  NumericsMatrix * M3 = MM[0];
  NumericsMatrix * M4 = MM[0];


  int info = -1;
  printf("== Numerics tests: prodNumericsMatrixNumericsMatrix(NumericsMatrix,NumericsMatrix) == \n");
  int i, j, k , n = M1->size1;
  double alpha = 1.0, beta = 1.0;
  double tol = 1e-12;


  NumericsMatrix C;

  C.storageType = 0;
  C.size0 = M1->size0;
  C.size1 = M1->size1;
  C.matrix0 = (double *)malloc(C.size0 * C.size1 * sizeof(double));
  C.matrix1 = NULL;
  prodNumericsMatrixNumericsMatrix(alpha, M1, M1, beta,  &C);

  double * Cref = (double *)malloc(C.size0 * C.size1 * sizeof(double));
  double sum;
  for (i = 0; i < C.size0; i++)
  {
    for (j = 0; j < C.size1; j++)
    {
      sum = 0.0;
      for (k = 0; k < M1->size1; k++)
      {
        sum = sum + M1->matrix0[i + k * M1->size1] * M1->matrix0[k + j * M1->size1];
      }
      Cref[i + j * C.size1] = sum;

    }
  }
  double err = 0.0;
  for (i = 0; i < C.size0; i++)
  {
    for (j = 0; j < C.size1; j++)
    {
      err += (Cref[i + j * C.size1] - C.matrix0[i + j * C.size1]) * (Cref[i + j * C.size1] - C.matrix0[i + j * C.size1]);
    }
  }
  if (err < tol)
  {
    info = 0;
  }

  if (info == 0)
    printf("Step 0 ( C = alpha*A*B + beta*C, double* storage, square matrix ) ok ...\n");
  else
  {
    printf("Step 0 ( C = alpha*A*B + beta*C, double* storage, square matrix) failed ...\n");
    return info;
  }


  NumericsMatrix C2;

  C2.storageType = 0;
  C2.size0 = M1->size0;
  C2.size1 = M3->size1;
  C2.matrix0 = (double *)malloc(C2.size0 * C2.size1 * sizeof(double));
  C2.matrix1 = NULL;
  prodNumericsMatrixNumericsMatrix(alpha, M1, M3, beta,  &C2);

  double * C2ref = (double *)malloc(C2.size0 * C2.size1 * sizeof(double));
  for (i = 0; i < C2.size0; i++)
  {
    for (j = 0; j < C2.size1; j++)
    {
      sum = 0.0;
      for (k = 0; k < M1->size1; k++)
      {
        sum = sum + M1->matrix0[i + k * M1->size1] * M3->matrix0[k + j * M3->size1];
      }
      C2ref[i + j * C2.size1] = sum;

    }
  }
  err = 0.0;
  for (i = 0; i < C2.size0; i++)
  {
    for (j = 0; j < C2.size1; j++)
    {
      err += (C2ref[i + j * C2.size1] - C2.matrix0[i + j * C2.size1]) * (C2ref[i + j * C2.size1] - C2.matrix0[i + j * C2.size1]);
    }
  }
  if (err < tol)
  {
    info = 0;
  }

  if (info == 0)
    printf("Step 1 ( C = alpha*A*B + beta*C, double* storage, non square) ok ...\n");
  else
  {
    printf("Step 1 ( C = alpha*A*B + beta*C, double* storage, non square) failed ...\n");
    return info;
  }







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

  int i, nmm = 4 ;
  NumericsMatrix ** NMM = malloc(nmm * sizeof(NumericsMatrix *)) ;


  for (i = 0 ; i < nmm; i++)
  {
    NMM[i] = malloc(sizeof(NumericsMatrix));
  }


  test_BuildNumericsMatrix(NMM);
  printf("Construction ok ...\n");
  int info = test_prodNumericsMatrix(NMM[0], NMM[1]);
  printf("End of ProdNumericsMatrix ...\n");
  info = test_prodNumericsMatrixNumericsMatrix(NMM);
  printf("End of ProdNumericsMatrixNumericsMatrix ...\n");
  info = test_subRowprod(NMM[0], NMM[1]);
  printf("End of Sub-Prod ...\n");
  info = test_rowProdNoDiag(NMM[0], NMM[1]);
  printf("End of Sub-Prod no diag ...\n");

  /* free memory */

  for (i = 0 ; i < nmm; i++)
  {
    free(NMM[i]->matrix0);
    /*    free(NMM[i]->matrix1->blocksize0); */
    /*    free(NMM[i]->matrix1->blocksize1); */
    /*    free(NMM[i]->matrix1->index1_data); */
    /*    free(NMM[i]->matrix1->index2_data); */
    free(NMM[i]->matrix1);
    free(NMM[i]);
  }
  free(NMM);



  printf("========= End Numerics tests for NumericsMatrix ========= \n");
  return info;
}

