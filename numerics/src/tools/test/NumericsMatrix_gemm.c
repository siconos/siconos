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

/*
  Tests functions for NumericsMatrix structure

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NumericsMatrix.h"
#include <math.h>
#include "numericsMatrixTestFunction.h"
#include "SiconosLapack.h"
#include "sanitizer.h"
#include "SparseBlockMatrix.h"
#include "NumericsSparseMatrix.h"
/* #define DEBUG_NOCOLOR */
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "debug.h"


/* ============================================================================================================================== */

static void add_initial_value_square_1(NumericsMatrix * M)
{
  int i=0, j=0;
  for (i=0; i < 4 ; i++)
  {
    for (j=0; j < 4 ; j++)
      NM_zentry(M,i,j,1.0+i+j);
  }
  for (i=0; i < 4 ; i++)
  {
    for (j=4; j < 6 ; j++)
      NM_zentry(M,i,j,2.0+i+j);
  }
  for (i=4; i < 6 ; i++)
  {
    for (j=4; j < 6 ; j++)
      NM_zentry(M,i,j,3.0+i+j);
  }
  for (i=4; i < 6 ; i++)
  {
    for (j=6; j < 8 ; j++)
      NM_zentry(M,i,j,4.0+i+j);
  }
  for (i=6; i < 8 ; i++)
  {
    for (j=0; j < 4 ; j++)
      NM_zentry(M,i,j,5.0+i+j);
  }
  for (i=6; i < 8 ; i++)
  {
    for (j=6; j < 8 ; j++)
      NM_zentry(M,i,j,6.0+i+j);
  }
  
}
static void add_initial_value_rectangle_1(NumericsMatrix * M)
{
  int i=0, j=0;
  for (i=0; i < 4 ; i++)
  {
    for (j=0; j < 4 ; j++)
      NM_zentry(M,i,j,1.0+i+j);
  }
  for (i=6; i < 8 ; i++)
  {
    for (j=0; j < 4 ; j++)
      NM_zentry(M,i,j,5.0+i+j);
  }
}

static int NM_gemm_test(NumericsMatrix** MM, double alpha, double beta)
{
  printf("\n == Numerics tests: NM_gemm(NumericsMatrix,NumericsMatrix) == \n");
  printf("Starts NM_gemm_test for alpha = %e and beta=%e\n",alpha,beta);

  NumericsMatrix * M1 = MM[0];
  NumericsMatrix * M2 = MM[1];
  NumericsMatrix * M3 = MM[2];
  NumericsMatrix * M4 = MM[3];


  int info = -1;
  int i, j, k;

  double tol = 1e-14;

  NumericsMatrix C;
  NM_null(&C);

  C.storageType = NM_DENSE;
  C.size0 = M1->size0;
  C.size1 = M1->size1;
  C.matrix0 = (double *)calloc(C.size0 * C.size1 , sizeof(double));
  MSAN_INIT_VAR(C.matrix0, C.size0 * C.size1);
  add_initial_value_square_1(&C);
  DEBUG_EXPR(NM_display(&C));

  NM_gemm(alpha, M1, M1, beta,  &C);

  NumericsMatrix * Cref= NM_create(NM_DENSE,C.size0, C.size1);
  add_initial_value_square_1(Cref);

  /* gemm by hand */
  double sum;
  for (i = 0; i < C.size0; i++)
  {
    for (j = 0; j < C.size1; j++)
    {
      sum = beta* Cref->matrix0[i + j * C.size0];
      for (k = 0; k < M1->size1; k++)
      {
        sum = sum + alpha* M1->matrix0[i + k * M1->size0] * M1->matrix0[k + j * M1->size0];
      }
      Cref->matrix0[i + j * C.size0] = sum;
    }
  }

  
  double err = 0.0;
  for (i = 0; i < C.size0; i++)
  {
    for (j = 0; j < C.size1; j++)
    {
      DEBUG_PRINTF("Cref[%i+%i*%i]= %lf\t\t", i, j, C.size0, Cref->matrix0[i + j * C.size0]);
      DEBUG_PRINTF("Cmatrix0[%i+%i*%i]= %lf\t", i, j, C.size0, C.matrix0[i + j * C.size0]);
      err += (Cref->matrix0[i + j * C.size0] - C.matrix0[i + j * C.size0]) * (Cref->matrix0[i + j * C.size0] - C.matrix0[i + j * C.size0]);
      DEBUG_PRINTF("err = %lf\n", err);
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
  C2.matrix0 = (double *)calloc(C2.size0 * C2.size1, sizeof(double));
  MSAN_INIT_VAR(C2.matrix0, C2.size0 * C2.size1);
  add_initial_value_rectangle_1(&C2);

  NM_gemm(alpha, M1, M3, beta,  &C2);

  NumericsMatrix * C2ref = NM_create(NM_DENSE,C2.size0,C2.size1);
  add_initial_value_rectangle_1(C2ref);
  

  for (i = 0; i < C2.size0; i++)
  {
    for (j = 0; j < C2.size1; j++)
    {
      sum = beta  * C2ref->matrix0[i + j * C2.size0] ;
      for (k = 0; k < M1->size1; k++)
      {
        sum = sum + alpha *  M1->matrix0[i + k * M1->size0] * M3->matrix0[k + j * M3->size0];
      }
      C2ref->matrix0[i + j * C2.size0] = sum;
    }
  }
  err = 0.0;
  for (i = 0; i < C2.size0; i++)
  {
    for (j = 0; j < C2.size1; j++)
    {
      err += (C2ref->matrix0[i + j * C2.size0] - C2.matrix0[i + j * C2.size0]) * (C2ref->matrix0[i + j * C2.size0] - C2.matrix0[i + j * C2.size0]);
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
  C3.storageType = NM_SPARSE_BLOCK;
  C3.size0 = M2->size0;
  C3.size1 = M2->size1;
  /* This step is not necessary for NM_gemm but conveniently creates a zero matrix with the right structure */
  C3.matrix1 = SBM_zero_matrix_for_multiply(M2->matrix1, M2->matrix1);
  add_initial_value_square_1(&C3);
  DEBUG_EXPR(SBM_print(C3.matrix1););

  NM_gemm(alpha, M2, M2, beta,  &C3);
  DEBUG_EXPR(NM_display(&C3));
  DEBUG_EXPR(NM_dense_display(Cref->matrix0,M2->size0,M2->size1,M2->size0));
  /*     Check if it is correct */

  /* C3 and CRef must have the same values.*/

  info = SBM_dense_equal(C3.matrix1, Cref->matrix0, tol);

  if (info == 0)
    printf("Step 2 ( C = alpha*A*B + beta*C, SBM storage) ok ...\n");
  else
  {
    printf("Step 2 ( C = alpha*A*B + beta*C, SBM storage) failed ...\n");
    goto exit_3;
  }

  NumericsMatrix C4;
  NM_null(&C4);
  C4.storageType = NM_SPARSE_BLOCK;
  C4.size0 = M2->size0;
  C4.size1 = M4->size1;
 /* This step is not necessary for NM_gemm but conveniently creates a zero matrix with the right structure */
  C4.matrix1 = SBM_zero_matrix_for_multiply(M2->matrix1, M4->matrix1);
  add_initial_value_rectangle_1(&C4);
  
  NM_gemm(alpha, M2, M4, beta,  &C4);

  DEBUG_EXPR(NM_display(&C4));
  DEBUG_EXPR(NM_dense_display(C2ref->matrix0,M2->size0,M4->size1,M2->size0));
  /*     Check if it is correct */
  /* C4 and C2Ref must have the same values.*/

  info = SBM_dense_equal(C4.matrix1, C2ref->matrix0, tol);

  if (info == 0)
    printf("Step 3 ( C = alpha*A*B + beta*C, SBM storage, non square) ok ...\n");
  else
  {
    printf("Step 3 ( C = alpha*A*B + beta*C, SBM storage, non square) failed ...\n");
    goto exit_4;
  }

  
  NumericsMatrix * M5 = test_matrix_5();
  DEBUG_EXPR(NM_display(M5););
  assert(NM_equal(M5,M2));

  NumericsMatrix * C5 = NM_create(NM_SPARSE,M5->size0, M5->size1);
  NM_triplet_alloc(C5,0);
  C5->matrix2->origin= NSM_TRIPLET;
  add_initial_value_square_1(C5);
  
  NM_gemm(alpha, M5, M5, beta, C5);

  info = NM_dense_equal(C5,Cref->matrix0,tol);
  
  if (info == 0)
    printf("Step 4 ( C = alpha*A*B + beta*C, NM_SPARSE storage, square) ok ...\n");
  else
  {
    printf("Step 4 ( C = alpha*A*B + beta*C, NM_SPARSE storage, square) failed ...\n");
    goto exit_5;
  }
  
  NumericsMatrix * M6 = test_matrix_6();
  DEBUG_EXPR(NM_display(M6););
  assert(NM_equal(M6,M4));

  NumericsMatrix * C6 = NM_create(NM_SPARSE,M2->size0, M4->size1);
  NM_triplet_alloc(C6,0);
  C6->matrix2->origin= NSM_TRIPLET;
  add_initial_value_rectangle_1(C6);
  
  NM_gemm(alpha, M2, M4, beta, C6);

  info = NM_dense_equal(C6,C2ref->matrix0,tol);
  
  if (info == 0)
    printf("Step 5 ( C = alpha*A*B + beta*C, NM_SPARSE storage, non square) ok ...\n");
  else
  {
    printf("Step 5 ( C = alpha*A*B + beta*C, NM_SPARSE storage, non square) failed ...\n");
  }

  
exit_5:
  NM_free(M5);
  NM_free(C5);
exit_4:
  NM_free(&C4);
exit_3:
  NM_free(&C3);
exit_2:
  free(C2.matrix0);
  NM_free(C2ref);
exit_1:
  NM_free(Cref);
  free(C.matrix0);
  return info;
}

int main(void)
{

  printf("========= Starts Numerics tests for NumericsMatrix ========= \n");

  int i, nmm = 4 ;
  NumericsMatrix ** NMM = (NumericsMatrix **)malloc(nmm * sizeof(NumericsMatrix *)) ;

  int info = test_build_first_4_NM(NMM);
  if (info != 0)
  {
    printf("Construction failed ...\n");
    return info;
  }
  printf("Construction ok ...\n");

  info = NM_gemm_test(NMM,1.0,0.0);
  if (info != 0)
  {
    printf("End of ProdNumericsMatrix : unsucessfull\n");
    return info;
  }
  info = NM_gemm_test(NMM,1.0,1.0);
  if (info != 0)
  {
    printf("End of ProdNumericsMatrix : unsucessfull\n");
    return info;
  }
  info = NM_gemm_test(NMM,0.0,1.0);
  if (info != 0)
  {
    printf("End of ProdNumericsMatrix : unsucessfull\n");
    return info;
  }
  info = NM_gemm_test(NMM,0.5,0.5);
  if (info != 0)
  {
    printf("End of ProdNumericsMatrix : unsucessfull\n");
    return info;
  }
  
  printf("End of ProdNumericsMatrix ...\n");
  if (info != 0) return info;
  /* free memory */

  for (i = 0 ; i < nmm; i++)
  {
    NM_free(NMM[i]);
    free(NMM[i]);
  }

  free(NMM);



  printf("========= End Numerics tests for NumericsMatrix ========= \n");
  return info;
}
