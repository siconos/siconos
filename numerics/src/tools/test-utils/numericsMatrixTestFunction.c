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

/*
  Tests functions for NumericsMatrix structure

 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "SiconosConfig.h"
#include "NumericsMatrix.h"
#include "SiconosLapack.h"
#include "numericsMatrixTestFunction.h"
#include "sanitizer.h"
#include "SparseBlockMatrix.h"
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
  M1->matrix0 = (double *)malloc(n * n * sizeof(double));

  int i;
  for (i = 0; i < n * n; i++)
    M1->matrix0[i] = m0[i];

  int nn = 4;
  /* Double * storage (column-major) */
  double m00[] = {1, 2, 0, 5, 0, 0, 0, 0,
                  2, 1, 0, 0, 0, 0, 0, 0,
                  0, 0, 1, -1, 0, 0, 2, 2,
                  4, 0, -1, 6, 0, 0, 1, 2
                 };
  M3->storageType = 0;
  M3->size0 = n;
  M3->size1 = nn;
  M3->matrix0 = (double *)malloc(nn * n * sizeof(double));
  for (i = 0; i < nn * n; i++)
    M3->matrix0[i] = m00[i];

  /* Build a NumericsMatrix with sparse-block storage */
  M2->storageType = 1;
  M2->size0 = n;
  M2->size1 = n;

  SparseBlockStructuredMatrix * SBM = (SparseBlockStructuredMatrix *)malloc(sizeof(SparseBlockStructuredMatrix));
  M2->matrix1 = SBM;
  SBM->nbblocks = 6;
  SBM->blocknumber0 = 3;
  SBM->blocknumber1 = SBM->blocknumber0 ;

  SBM->blocksize0 = (unsigned int*)malloc(3 * sizeof(unsigned int));
  SBM->blocksize0[0] = 4;
  SBM->blocksize0[1] = 6;
  SBM->blocksize0[2] = 8;
  SBM->blocksize1 = SBM->blocksize0;


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

  SBM->block = (double **)malloc(SBM->nbblocks * sizeof(*(SBM->block)));
  double block0[16] = {1, 2 , 0 , 5 , 2 , 1 , 0 , 0 , 0 , 0 , 1 , -1, 4, 0 , -1, 6};
  double block1[8] = {3, 4, 0, 0, -1, 1, 0, 6};
  double block2[4] = {1, 0, 0, 2};
  double block3[4] = {0, 0, 5, 2};
  double block4[8] = {0, 0, 0, 0, 2, 2, 1, 2};
  double block5[4] = {2, -1, 2, 2};
  SBM->block[0] = (double *)malloc(16 * sizeof(double));
  SBM->block[1] = (double *)malloc(8 * sizeof(double));
  SBM->block[2] = (double *)malloc(4 * sizeof(double));
  SBM->block[3] = (double *)malloc(4 * sizeof(double));
  SBM->block[4] = (double *)malloc(8 * sizeof(double));
  SBM->block[5] = (double *)malloc(4 * sizeof(double));
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
  /* SBM_print(SBM); */

  /* M1 and M2 must have the same values.*/
  double tol = 1e-12;
  int j;
  for (i = 0; i < M1->size0; i++)
  {
    for (j = 0; j < M1->size1; j++)
    {
      if (fabs(M1->matrix0[i + j * M1->size0] - SBM_get_value(M2->matrix1, i, j)) > tol) info = 1;

      /*    printf("%i\t%i\n",i,j); */
      /*    printf("%lf\n",M1->matrix0[i+j*M1->size0]-SBM_get_value(M2->matrix1,i,j));  */
      /*    printf("%lf\n",M1->matrix0[i+j*M1->size0]); */
      /*    printf("%lf\n",SBM_get_value(M2->matrix1,i,j));   */

      if (info == 1) break;
    }
    if (info == 1) break ;
  }


  /* Build a NumericsMatrix with sparse-block storage */
  M4->storageType = 1;
  M4->size0 = n;
  M4->size1 = 4;

  M4->matrix1 = (SparseBlockStructuredMatrix *)malloc(sizeof(SparseBlockStructuredMatrix));
  SparseBlockStructuredMatrix * SBM2 = M4->matrix1;

  SBM2->nbblocks = 2;
  SBM2->blocknumber0 = 3;

  SBM2->blocksize0 = (unsigned int*)malloc(SBM2->blocknumber0 * sizeof(unsigned int));
  SBM2->blocksize0[0] = 4;
  SBM2->blocksize0[1] = 6;
  SBM2->blocksize0[2] = 8;
  SBM2->blocksize1 = SBM2->blocksize0;

  SBM2->blocknumber1 = 1;
  SBM2->blocksize1 = (unsigned int*)malloc(SBM2->blocknumber1 * sizeof(unsigned int));
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

  SBM2->block = (double **)malloc(SBM2->nbblocks * sizeof(*(SBM2->block)));
  double block00[16] = {1, 2 , 0 , 5 , 2 , 1 , 0 , 0 , 0 , 0 , 1 , -1, 4, 0 , -1, 6};
  double block40[8] = {0, 0, 0, 0, 2, 2, 1, 2};
  SBM2->block[0] = (double *)malloc(16 * sizeof(double));
  SBM2->block[1] = (double *)malloc(8 * sizeof(double));
  for (i = 0; i < 16; i++)
    SBM2->block[0][i] = block00[i];
  for (i = 0; i < 8; i++)
    SBM2->block[1][i] = block40[i];


  /*   SBM_print(SBM2); */
  /* M3 and M4 must have the same values.*/

  for (i = 0; i < M3->size0; i++)
  {
    for (j = 0; j < M3->size1; j++)
    {
      if (fabs(M3->matrix0[i + j * M3->size0] - SBM_get_value(M4->matrix1, i, j)) > tol) info = 1;

      /*    printf("%i\t%i\n",i,j); */
      /*    printf("%lf\n",M3->matrix0[i+j*M3->size0]-SBM_get_value(M4->matrix1,i,j)); */
      /*    printf("%lf\n",M3->matrix0[i+j*M3->size0]); */
      /*    printf("%lf\n",SBM_get_value(M4->matrix1,i,j)); */

      if (info == 1) break;
    }
    if (info == 1) break ;
  }
  return info;
}
/* ============================================================================================================================== */

int test_prodNumericsMatrix(NumericsMatrix** MM)
{
  NumericsMatrix* M1 =  MM[0];
  NumericsMatrix* M2 =  MM[1];
  NumericsMatrix* M3 =  MM[2];
  NumericsMatrix* M4 =  MM[3];




  printf("== Numerics tests: prodNumericsMatrix(NumericsMatrix,vector) == \n");
  int i , n = M1->size1, m = 4;

  double * x = (double *)malloc(n * sizeof(double));
  double * x2 = (double *)malloc(m * sizeof(double));
  double alpha = 2.3, beta = 1.9;
  double * yref = (double *)malloc(n * sizeof(double));
  double * yref2 = (double *)malloc(n * sizeof(double));;
  double * y = (double *)malloc(n * sizeof(double));
  double * y2 = (double *)malloc(n * sizeof(double));
  for (i = 0; i < n; i++)
  {
    x[i] = i + 1.0;
    yref[i] = 0.1 * i;
    yref2[i] = 0.1 * i;
    y[i] = yref[i];
    y2[i] = yref2[i];
  }
  x2[0] = 0;
  x2[1] = 0;
  x2[2] = 0;
  x2[3] = 0;
  int incx = 1, incy = 1;
  cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, alpha, M1->matrix0, n, x, incx, beta, yref, incy);

  NM_gemv(alpha, M1, x, beta, y);
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


  cblas_dgemv(CblasColMajor, CblasNoTrans, n, m, alpha, M3->matrix0, n, x2, incx, beta, yref2, incy);

  NM_gemv(alpha, M3, x2, beta, y2);
  for (i = 0; i < n; i++)
  {
    if (fabs(y2[i] - yref2[i]) > tol) info = 1;
    /*           printf("%lf\n", fabs(y2[i]-yref2[i])); */
    /*            printf("%lf\n",y2[i]); */
    /*            printf("%lf\n",yref2[i]); */
  }
  if (info == 0)
    printf("Step 1 ( y = alpha*A*x + beta*y, double* storage, non square) ok ...\n");
  else
    printf("Step 1 ( y = alpha*A*x + beta*y, double* storage, non square) failed ...\n");






  /* Sparse ... */
  for (i = 0; i < n; i++)
  {
    y[i] = 0.1 * i;
    y2[i] = 0.1 * i;
  }
  NM_gemv(alpha, M2, x, beta, y);
  for (i = 0; i < n; i++)
  {
    if (fabs(y[i] - yref[i]) > tol) info = 1;
    /*       printf("%lf\n", fabs(y[i]-yref[i]));  */
    /*       printf("%lf\n", y[i]);  */
  }
  if (info == 0)
    printf("Step 2 ( y = alpha*A*x + beta*y, sparse storage) ok ...\n");
  else
    printf("Step 2 ( y = alpha*A*x + beta*y,  sparsestorage) failed ...\n");



  NM_gemv(alpha, M4, x2, beta, y2);
  for (i = 0; i < n; i++)
  {

    if (fabs(y2[i] - yref2[i]) > tol) info = 1;
    /*      printf("%lf\n", fabs(y2[i]-yref2[i]));  */
    /*       printf("%lf\n",y2[i]); */
    /*       printf("%lf\n",yref2[i]); */
  }

  if (info == 0)
    printf("Step 3 ( y = alpha*A*x + beta*y, sparse storage, non square) ok ...\n");
  else
    printf("Step 3 ( y = alpha*A*x + beta*y,  sparsestorage, non square) failed ...\n");



  free(x);
  free(x2);
  free(y);
  free(y2);
  free(yref);
  free(yref2);

  printf("== End of test NM_gemv(result = %d\n", info);

  return info;
}

/* ============================================================================================================================== */

int test_prodNumericsMatrixNumericsMatrix(NumericsMatrix** MM)
{


  NumericsMatrix * M1 = MM[0];
  NumericsMatrix * M2 = MM[1];
  NumericsMatrix * M3 = MM[2];
  NumericsMatrix * M4 = MM[3];


  int info = -1;
  printf("== Numerics tests: NM_gemm(NumericsMatrix,NumericsMatrix) == \n");
  int i, j, k;
  double alpha = 1.0, beta = 0.0;
  double tol = 1e-12;

  double * C2ref = 0;

  SparseBlockStructuredMatrix * SBM3 = 0;
  SparseBlockStructuredMatrix * SBM4 = 0;

  NumericsMatrix C;
  NM_null(&C);

  C.storageType = 0;
  C.size0 = M1->size0;
  C.size1 = M1->size1;
  C.matrix0 = (double *)malloc(C.size0 * C.size1 * sizeof(double));
  MSAN_INIT_VAR(C.matrix0, C.size0 * C.size1);
  NM_gemm(alpha, M1, M1, beta,  &C);

  double * Cref = (double *)malloc(C.size0 * C.size1 * sizeof(double));
  double sum;
  for (i = 0; i < C.size0; i++)
  {
    for (j = 0; j < C.size1; j++)
    {
      sum = 0.0;
      for (k = 0; k < M1->size1; k++)
      {
        sum = sum + M1->matrix0[i + k * M1->size0] * M1->matrix0[k + j * M1->size0];
      }
      Cref[i + j * C.size0] = sum;

    }
  }
  double err = 0.0;
  for (i = 0; i < C.size0; i++)
  {
    for (j = 0; j < C.size1; j++)
    {
      printf("Cref[%i+%i*%i]= %lf\t\t", i, j, C.size0, Cref[i + j * C.size0]);
      printf("Cmatrix0[%i+%i*%i]= %lf\t", i, j, C.size0, C.matrix0[i + j * C.size0]);
      err += (Cref[i + j * C.size0] - C.matrix0[i + j * C.size0]) * (Cref[i + j * C.size0] - C.matrix0[i + j * C.size0]);
      printf("err = %lf\n", err);
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
    goto exit_1;
  }


  NumericsMatrix C2;
  NM_null(&C2);

  C2.storageType = 0;
  C2.size0 = M1->size0;
  C2.size1 = M3->size1;
  C2.matrix0 = (double *)malloc(C2.size0 * C2.size1 * sizeof(double));
  MSAN_INIT_VAR(C2.matrix0, C2.size0 * C2.size1);
  NM_gemm(alpha, M1, M3, beta,  &C2);

  C2ref = (double *)malloc(C2.size0 * C2.size1 * sizeof(double));
  for (i = 0; i < C2.size0; i++)
  {
    for (j = 0; j < C2.size1; j++)
    {
      sum = 0.0;
      for (k = 0; k < M1->size1; k++)
      {
        sum = sum + M1->matrix0[i + k * M1->size0] * M3->matrix0[k + j * M3->size0];
      }
      C2ref[i + j * C2.size0] = sum;
      /*    printf("C2ref(%i,%i)=%f\n", i,j,C2ref[i+j*C2.size0] ); */
    }
  }
  err = 0.0;
  for (i = 0; i < C2.size0; i++)
  {
    for (j = 0; j < C2.size1; j++)
    {
      err += (C2ref[i + j * C2.size0] - C2.matrix0[i + j * C2.size0]) * (C2ref[i + j * C2.size0] - C2.matrix0[i + j * C2.size0]);
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
    goto exit_2;
  }

  NumericsMatrix C3;
  NM_null(&C3);

  C3.storageType = 1;
  C3.size0 = M2->size0;
  C3.size1 = M2->size1;
  SBM3 = (SparseBlockStructuredMatrix *)malloc(sizeof(SparseBlockStructuredMatrix));
  C3.matrix1 = SBM3;
  beta = 1.0;
  i = 1;
  // while (i > 0)

  SBM_alloc_for_gemm(M2->matrix1, M2->matrix1, SBM3);
  /* SBM_print(SBM3); */


  NM_gemm(alpha, M2, M2, beta,  &C3);
  //SBM_free(SBM3);
  printf("i= %i\n", i++);


  /*     Check if it is correct */

  /* C3 and CRef must have the same values.*/

  for (i = 0; i < C3.size0; i++)
  {
    for (j = 0; j < C3.size1; j++)
    {
      if (fabs(Cref[i + j * C3.size0] - SBM_get_value(C3.matrix1, i, j)) > tol) info = 1;

      /*    printf("%i\t%i\n",i,j); */
      /*    printf("%lf\n",fabs(Cref[i+j*C3.size0]-SBM_get_value(C3.matrix1,i,j) )); */
      /*     printf("%lf\n",Cref[i+j*C3.size0]);  */
      /*     printf("%lf\n",SBM_get_value(C3.matrix1,i,j));  */

      if (info == 1) break;
    }
    if (info == 1) break ;
  }


  if (info == 0)
    printf("Step 2 ( C = alpha*A*B + beta*C, sparse storage) ok ...\n");
  else
  {
    printf("Step 2 ( C = alpha*A*B + beta*C, sparse storage) failed ...\n");
    goto exit_3;
  }

  NumericsMatrix C4;

  NM_null(&C4);
  C4.storageType = 1;
  C4.size0 = M2->size0;
  C4.size1 = M4->size1;
  SBM4 = (SparseBlockStructuredMatrix *)malloc(sizeof(SparseBlockStructuredMatrix));
  C4.matrix1 = SBM4;

  SBM_alloc_for_gemm(M2->matrix1, M4->matrix1, SBM4);
  /* SBM_print(SBM4); */

  NM_gemm(alpha, M2, M4, beta,  &C4);


  /*     Check if it is correct */

  /* C4 and C2Ref must have the same values.*/

  for (i = 0; i < C4.size0; i++)
  {
    for (j = 0; j < C4.size1; j++)
    {
      if (fabs(C2ref[i + j * C4.size0] - SBM_get_value(C4.matrix1, i, j)) > tol) info = 1;

      /*    printf("%i\t%i\n",i,j); */
      /*    printf("%lf\n",fabs(C2ref[i+j*C4.size0]-SBM_get_value(C4.matrix1,i,j) )); */
      /*     printf("%lf\n",C2ref[i+j*C4.size0]); */
      /*     printf("%lf\n",SBM_get_value(C4.matrix1,i,j)); */

      if (info == 1) break;
    }
    if (info == 1) break ;
  }


  if (info == 0)
    printf("Step 3 ( C = alpha*A*B + beta*C, sparse storage) ok ...\n");
  else
  {
    printf("Step 3 ( C = alpha*A*B + beta*C, sparse storage) failed ...\n");
  }

  NM_free(&C4);
exit_3:
  NM_free(&C3);
exit_2:
  free(C2.matrix0);
  free(C2ref);
exit_1:
  free(Cref);
  free(C.matrix0);
  return info;
}

int test_NM_row_prod(NumericsMatrix* M1, NumericsMatrix* M2)
{

  printf("== Numerics tests: NM_row_prod(NumericsMatrix,vector) == \n");
  int i , n = M1->size1;
  double * x = (double *)malloc(n * sizeof(double));

  for (i = 0; i < n; i++)
  {
    x[i] = i + 1;
  }

  int min = 2;
  int max = 6;
  int sizeY = max - min;
  /* Computes yRef = subA*x, subA = A limited to row min to max*/
  double * y = (double *)malloc(sizeY * sizeof(double));
  double yref[4];
  int incx = n, incy = 1;
  for (i = 0; i < sizeY; i++)
    yref[i] = cblas_ddot(n, &(M1->matrix0[min + i]), incx, x, incy);

  NM_row_prod(n, sizeY, min, M1, x, y, 1);
  double tol = 1e-12;
  int info = 0;
  for (i = 0; i < sizeY; i++)
  {
    if (fabs(y[i] - yref[i]) > tol) info = 1;
    /*        printf("%lf\n", fabs(y[i]-yref[i]));  */
    /*           printf("%lf\n", y[i]); */
    /*           printf("%lf\n", yref[i]); */
  }
  if (info == 0)
    printf("Step 0 ( y = subA*x, double* storage) ok ...\n");
  else
    printf("Step 0 ( y = subA*x, double* storage) failed ...\n");

  /* += */
  NM_row_prod(n, sizeY, min, M1, x, y, 0);
  for (i = 0; i < sizeY; i++)
  {
    if (fabs(y[i] - 2 * yref[i]) > tol) info = 1;
    /*        printf("%lf\n", fabs(y[i]-2*yref[i]));  */
    /*           printf("%lf\n", y[i]); */
    /*           printf("%lf\n", 2*yref[i]); */
  }
  if (info == 0)
    printf("Step 1 ( y += subA*x, double* storage) ok ...\n");
  else
    printf("Step 1 ( y += subA*x, double* storage) failed ...\n");





  free(y);
  sizeY = 2;
  int pos = 1; // pos of the required row of blocks
  y = (double *)malloc(sizeY * sizeof(double));

  for (i = 0; i < sizeY; i++)
  {
    y[i] = 0.0;
    yref[i] = cblas_ddot(n, &(M1->matrix0[4 + i]), incx, x, incy);
  }
  /* Sparse ... */
  NM_row_prod(n, sizeY, pos, M2, x, y, 1);
  for (i = 0; i < sizeY; i++)
  {
    if (fabs(y[i] - yref[i]) > tol) info = 1;
    //  printf("%lf\n", fabs(y[i]-yref[i]));
  }
  for (i = 0; i < sizeY; i++)
    yref[i] = cblas_ddot(n, &(M1->matrix0[6 + i]), incx, x, incy);
  NM_row_prod(n, sizeY, pos + 1, M2, x, y, 1);
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
  NM_row_prod(n, sizeY, pos + 1, M2, x, y, 0);
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
  printf("== End of test NM_row_prod(NumericsMatrix,vector), result = %d\n", info);

  return info;
}
int test_NM_row_prod_non_square(NumericsMatrix* M3, NumericsMatrix* M4)
{

  printf("== Numerics tests: subRowProd_non_square(NumericsMatrix,vector) == \n");
  int i , n = M3->size0, m = M3->size1;
  double * x = (double *)malloc(m * sizeof(double));

  for (i = 0; i < m; i++)
  {
    x[i] = i + 1;
  }

  int min = 1;
  int max = 3;
  int sizeY = max - min;
  int sizeX = m;
  /* Computes yRef = subA*x, subA = A limited to row min to max*/
  double * y = (double *)malloc(sizeY * sizeof(double));
  double yref[2];
  int incx = n, incy = 1;
  for (i = 0; i < sizeY; i++)
    yref[i] = cblas_ddot(m, &(M3->matrix0[min + i]), incx, x, incy);

  NM_row_prod(sizeX, sizeY, min, M3, x, y, 1);
  double tol = 1e-12;
  int info = 0;
  for (i = 0; i < sizeY; i++)
  {
    if (fabs(y[i] - yref[i]) > tol) info = 1;
    /*       printf("%lf\n", fabs(y[i]-yref[i]));  */
    /*           printf("%lf\n", y[i]); */
    /*           printf("%lf\n", yref[i]); */
  }
  if (info == 0)
    printf("Step 0 ( y = subA*x, double* storage _non_square) ok ...\n");
  else
    printf("Step 0 ( y = subA*x, double* storage _non_square) failed ...\n");

  /* += */
  NM_row_prod(sizeX, sizeY, min, M3, x, y, 0);
  for (i = 0; i < sizeY; i++)
  {
    if (fabs(y[i] - 2 * yref[i]) > tol) info = 1;
    /*         printf("%lf\n", fabs(y[i]-2*yref[i]));  */
    /*           printf("%lf\n", y[i]); */
    /*           printf("%lf\n", 2*yref[i]); */
  }
  if (info == 0)
    printf("Step 1 ( y += subA*x, double* storage _non_square) ok ...\n");
  else
    printf("Step 1 ( y += subA*x, double* storage _non_square) failed ...\n");


  free(y);


  sizeY = 2;
  int pos = 1; // pos of the required row of blocks
  y = (double *)malloc(sizeY * sizeof(double));

  for (i = 0; i < sizeY; i++)
  {
    y[i] = 0.0;
    yref[i] = cblas_ddot(m, &(M3->matrix0[4 + i]), incx, x, incy);
  }
  /* Sparse ... */
  NM_row_prod(sizeX, sizeY, pos, M4, x, y, 1);
  for (i = 0; i < sizeY; i++)
  {
    if (fabs(y[i] - yref[i]) > tol) info = 1;
    //  printf("%lf\n", fabs(y[i]-yref[i]));
  }
  for (i = 0; i < sizeY; i++)
    yref[i] = cblas_ddot(m, &(M3->matrix0[6 + i]), incx, x, incy);
  NM_row_prod(sizeX, sizeY, pos + 1, M4, x, y, 1);
  for (i = 0; i < sizeY; i++)
  {
    if (fabs(y[i] - yref[i]) > tol) info = 1;
    //  printf("%lf\n", fabs(y[i]-yref[i]));
  }


  if (info == 0)
    printf("Step 2 ( y = subA*x, sparse storage _non_square) ok ...\n");
  else
    printf("Step 2 ( y = subA*x,  sparse storage _non_square) failed ...\n");

  /* Sparse, += ... */
  NM_row_prod(sizeX, sizeY, pos + 1, M4, x, y, 0);
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
  printf("== End of test NM_row_prod(NumericsMatrix,vector), result = %d\n", info);

  return info;
}

int test_NM_row_prod_no_diag(NumericsMatrix* M1, NumericsMatrix* M2)
{

  printf("== Numerics tests: NM_row_prod_no_diag(NumericsMatrix,vector) == \n");
  int i , n = M1->size1;
  double * x = (double *)malloc(n * sizeof(double));

  for (i = 0; i < n; i++)
  {
    x[i] = i + 1;
  }

  int min = 2;
  int max = 6;
  int sizeY = max - min;
  /* Computes yRef = subA*x, subA = A limited to row min to max*/
  double * y = (double *)malloc(sizeY * sizeof(double));
  double yref[4];
  //  int incx = n, incy =1;
  double tol = 1e-12;
  int info = 0;
  /*   int incx=1,incy=1; */

  /*   for(i=0;i<sizeY;i++) */
  /*     yref[i]=cblas_ddot(n, &(M1->matrix0[min+i]), incx, x, incy); */

  /*   NM_row_prod_no_diag(n,sizeY,min,M1,x,y,1); */

  /*   for(i = 0; i< sizeY; i++) */
  /*     { */
  /*       if( fabs(y[i]-yref[i])>tol) info = 1; */
  /* /\*       printf("%lf\n", fabs(y[i]-yref[i])); *\/ */
  /*     } */
  /*   if(info ==0) */
  /*     printf("Step 0 ( y = subA*x, double* storage) ok ...\n"); */
  /*   else */
  /*     printf("Step 0 ( y = subA*x, double* storage) failed ...\n"); */

  /*   /\* += *\/ */
  /*   NM_row_prod_no_diag(n,sizeY,min,M1,x,y,0); */
  /*   for(i = 0; i< sizeY; i++) */
  /*     { */
  /*       if( fabs(y[i]-2*yref[i])>tol) info = 1; */
  /*       /\*       printf("%lf\n", fabs(y[i]-2*yref[i]));  *\/ */
  /*     } */
  /*   if(info ==0) */
  /*     printf("Step 1 ( y += subA*x, double* storage) ok ...\n"); */
  /*   else */
  /*     printf("Step 1 ( y += subA*x, double* storage) failed ...\n"); */

  free(y);
  sizeY = 2;
  int pos = 1; // pos of the required row of blocks
  y = (double *)calloc(sizeY, sizeof(double));


  yref[0] = 40;
  yref[1] = 16;

  /* Sparse ... */
  NM_row_prod_no_diag(n, sizeY, pos, SIZE_MAX, M2, x, y, NULL, 1);
  for (i = 0; i < sizeY; i++)
  {
    if (fabs(y[i] - yref[i]) > tol) info = 1;
    //  printf("%lf\n", fabs(y[i]-yref[i]));
  }
  NM_row_prod_no_diag(n, sizeY, pos + 1, SIZE_MAX, M2, x, y, NULL, 1);
  yref[0] = 10;
  yref[1] = 14;
  for (i = 0; i < sizeY; i++)
  {
    if (fabs(y[i] - yref[i]) > tol) info = 1;
  }

  if (info == 0)
    printf("Step 2 ( y = subA*x, sparse storage) ok ...\n");
  else
    printf("Step 2 ( y = subA*x,  sparse storage) failed ...\n");

  /* Sparse, += ... */
  NM_row_prod_no_diag(n, sizeY, pos + 1, SIZE_MAX, M2, x, y, NULL, 0);
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
  printf("== End of test NM_row_prod_no_diag(NumericsMatrix,vector), result = %d\n", info);

  return info;
}
int test_NM_row_prod_no_diag_non_square(NumericsMatrix* M3, NumericsMatrix* M4)
{

  printf("== Numerics tests: NM_row_prod_no_diag_non_square(NumericsMatrix,vector) == \n");
  int i ,  m = M3->size1;
  double * x = (double *)malloc(m * sizeof(double));

  for (i = 0; i < m; i++)
  {
    x[i] = i + 1;
  }

  int min = 2;
  int max = 6;
  int sizeY = max - min;
  int sizeX = m;
  /* Computes yRef = subA*x, subA = A limited to row min to max*/
  double * y = (double *)malloc(sizeY * sizeof(double));
  double yref[4];
  //  int incx = n, incy =1;
  double tol = 1e-12;
  int info = 0;
  /*   for(i=0;i<sizeY;i++) */
  /*     yref[i]=cblas_ddot(n, &(M3->matrix0[min+i]), incx, x, incy); */

  /*   NM_row_prod_no_diag(n,sizeY,min,M3,x,y,1); */
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
  /*   NM_row_prod_no_diag(n,sizeY,min,M3,x,y,0); */
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
  y = (double *)malloc(sizeY * sizeof(double));
  y[0] = 0;
  y[1] = 0;
  yref[0] = 0;
  yref[1] = 0;

  /* Sparse ... */
  NM_row_prod_no_diag(sizeX, sizeY, pos, SIZE_MAX, M4, x, y, NULL, 1);
  for (i = 0; i < sizeY; i++)
  {
    if (fabs(y[i] - yref[i]) > tol) info = 1;
    /*       printf("%lf\n", fabs(y[i]-yref[i])); */
    /*       printf("%lf\n", y[i]); */
    /*       printf("%lf\n", yref[i]); */
  }

  printf("\n");
  yref[0] = 10;
  yref[1] = 14;

  NM_row_prod_no_diag(sizeX, sizeY, pos + 1, SIZE_MAX, M4, x, y, NULL, 1);

  for (i = 0; i < sizeY; i++)
  {
    if (fabs(y[i] - yref[i]) > tol) info = 1;
    /*      printf("%lf\n", fabs(y[i]-yref[i])); */
    /*             printf("%lf\n", y[i]); */
    /*             printf("%lf\n", yref[i]); */
  }

  if (info == 0)
    printf("Step 2 ( y = subA*x, sparse storage) ok ...\n");
  else
    printf("Step 2 ( y = subA*x,  sparse storage) failed ...\n");

  /* Sparse, += ... */
  NM_row_prod_no_diag(sizeX, sizeY, pos + 1, SIZE_MAX, M4, x, y, NULL, 0);
  for (i = 0; i < sizeY; i++)
  {
    if (fabs(y[i] - 2 * yref[i]) > tol) info = 1;
    /*           printf("%lf\n", fabs(y[i]-yref[i])); */
    /*             printf("%lf\n", y[i]); */
    /*             printf("%lf\n", yref[i]); */
  }
  if (info == 0)
    printf("Step 3 ( y += subA*x, sparse storage) ok ...\n");
  else
    printf("Step 3 ( y += subA*x,  sparse storage) failed ...\n");


  free(x);
  free(y);
  printf("== End of test NM_row_prod_no_diag_non_square(NumericsMatrix,vector), result = %d\n", info);

  return info;
}
int test_SBM_row_to_dense(SparseBlockStructuredMatrix *M)
{
  double * denseRes = (double*) malloc(M->blocksize0[M->blocknumber0 - 1] * M->blocksize1[M->blocknumber1 - 1] * sizeof(double));
  unsigned int curRow = 0;
  unsigned int nbCol = M->blocksize1[M->blocknumber1 - 1];
  for (unsigned int i = 0; i < M->blocknumber0; i++)
  {
    unsigned int lLin = 0;
    unsigned int nbBlockRow = M->blocksize0[i] - curRow;
    SBM_row_to_dense(M, i, denseRes, 0, nbBlockRow);
    for (unsigned int lin = curRow; lin < M->blocksize0[i]; lin++)
    {
      unsigned int lCol = 0;
      for (unsigned int col = 0; col < nbCol; col++)
      {
        if (fabs(SBM_get_value(M, lin, col) - denseRes[lLin + lCol * (nbBlockRow)]) > 10e-12)
        {
          free(denseRes);
          return 1;
        }
        lCol++;
      }
      lLin++;
    }
    curRow = M->blocksize0[i];
  }
  curRow = 0;
  for (unsigned int i = 0; i < M->blocknumber0; i++)
  {

    //    int lLin=0;
    //    int nbBlockRow=M->blocksize0[i]-curRow;
    SBM_row_to_dense(M, i, denseRes, curRow, M->blocksize0[M->blocknumber0 - 1]);
    curRow = M->blocksize0[i];
  }

  double * denseRes2 = (double*) malloc(M->blocksize0[M->blocknumber0 - 1] * M->blocksize1[M->blocknumber1 - 1] * sizeof(double));

  SBM_to_dense(M, denseRes2);
  for (unsigned int n = 0; n < M->blocksize0[M->blocknumber0 - 1]*M->blocksize1[M->blocknumber1 - 1]; n++)
    if (fabs(denseRes2[n] - denseRes[n]) > 10e-12)
    {
      free(denseRes);
      free(denseRes2);
      return 1;
    }

  free(denseRes);
  free(denseRes2);
  return 0;
}
int test_SBM_row_permutation(SparseBlockStructuredMatrix *M)
{
  SparseBlockStructuredMatrix MRes;
  unsigned int nbRow = M->blocknumber0;
  unsigned int * rowBlockIndex = (unsigned int*) malloc(nbRow * sizeof(unsigned int));
  unsigned int * mark = (unsigned int*) malloc(nbRow * sizeof(unsigned int));
  for (unsigned int i = 0; i < nbRow; i++)
  {
    mark[i] = 0;
  }
  for (unsigned int i = 0; i < nbRow; i++)
  {
    int candidate = rand() % nbRow;
    while (mark[candidate])
      candidate = rand() % nbRow;
    rowBlockIndex[i] = candidate;
    mark[candidate] = 1;
  }
  SBM_row_permutation(rowBlockIndex, M, &MRes);
  double * denseMRes = (double*) malloc(M->blocksize0[M->blocknumber0 - 1] * M->blocksize1[M->blocknumber1 - 1] * sizeof(double));
  SBM_to_dense(&MRes, denseMRes);
  double * denseM = (double*) malloc(M->blocksize0[M->blocknumber0 - 1] * M->blocksize1[M->blocknumber1 - 1] * sizeof(double));
  unsigned int curRow = 0;
  unsigned int nbRowInM = M->blocksize0[M->blocknumber0 - 1];
  for (unsigned int i = 0; i < nbRow; i++)
  {
    unsigned int rowInM = rowBlockIndex[i];
    unsigned int nbRow = 0;
    if (rowInM)
      nbRow = M->blocksize0[rowInM] - M->blocksize0[rowInM - 1];
    else
      nbRow = M->blocksize0[rowInM];
    SBM_row_to_dense(M, rowInM, denseM, curRow, nbRowInM);
    curRow += nbRow;
  }
  for (unsigned int n = 0; n < M->blocksize0[M->blocknumber0 - 1]*M->blocksize1[M->blocknumber1 - 1]; n++)
    if (fabs(denseMRes[n] - denseM[n]) > 10e-12)
    {
      free(denseM);
      free(denseMRes);
      free(rowBlockIndex);
      free(mark);
      SBMfree(&MRes, 0);
      return 1;
    }
  free(denseM);
  free(denseMRes);
  free(rowBlockIndex);
  free(mark);
  SBMfree(&MRes, 0);
  return 0;
}
int test_SBM_column_permutation(SparseBlockStructuredMatrix *M)
{
  SparseBlockStructuredMatrix MRes;
  int nbCol = M->blocknumber1;
  unsigned int * colBlockIndex = (unsigned int*) malloc(nbCol * sizeof(unsigned int));
  int * mark = (int*) malloc(nbCol * sizeof(int));
  for (int i = 0; i < nbCol; i++)
    mark[i] = 0;
  for (int i = 0; i < nbCol; i++)
  {
    int candidate = rand() % nbCol;
    while (mark[candidate])
      candidate = rand() % nbCol;
    colBlockIndex[i] = candidate;
    mark[candidate] = 1;
  }
  SBM_column_permutation(colBlockIndex, M, &MRes);
  free(colBlockIndex);
  free(mark);
  SBMfree(&MRes, 0);
  return 0;
}
