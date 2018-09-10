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
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
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
static void add_initial_value_square_2(NumericsMatrix * M)
{
  int i=0, j=0;
  for (i=0; i < 4 ; i++)
  {
    for (j=0; j < 4 ; j++)
      NM_zentry(M,i,j,1.0+i+j);
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

/* hand made gemm of nxp A matrix and pxm B matrix */
static void dense_gemm_by_hand(double alpha, double * A, double * B, int n, int m, int p, double beta, double *C)
{
  double sum =0.0;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
    {
      sum = beta  * C[i + j * n] ;
      for (int k = 0; k < p ; k++)
      {
        sum = sum + alpha *  A[i + k * n] * B[k + j * p];
      }
      C[i + j * n] = sum;
    }
  }
}

static double dense_comparison(double * C, int n, int m, double *Cref)
{
  double err = 0.0;
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m ; j++)
    {
      DEBUG_PRINTF("Cref[%i+%i*%i]= %lf\t\t", i, j, n, Cref[i + j * n]);
      DEBUG_PRINTF("C[%i+%i*%i]= %lf\t", i, j, n, C[i + j * n]);
      err += (Cref[i + j * n] - C[i + j * n]) * (Cref[i + j * n] - C[i + j * n]);
      DEBUG_PRINTF("err = %lf\n", err);
    }
  }
  return err;
}


static int SBM_gemm_without_allocation_test(NumericsMatrix** MM, double alpha, double beta)
{
  printf("\n == Numerics tests: SBM_gemm_without_allocation_test == \n");
  printf("Starts SBM_gemm_without_allocation_test for alpha = %e and beta=%e\n",alpha,beta);

  NumericsMatrix * M1 = MM[0];
  NumericsMatrix * M2 = MM[1];
  NumericsMatrix * M3 = MM[2];
  NumericsMatrix * M4 = MM[3];


  int info = -1;


  double tol = 1e-14;

  /***********************************************************/
  /* C = alpha*A*B + beta*C, double* storage, square matrix  */
  /***********************************************************/

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
  dense_gemm_by_hand(alpha, M1->matrix0, M1->matrix0, M1->size0, M1->size1, M1->size0, beta,  Cref->matrix0);
  double err = dense_comparison(C.matrix0, C.size0, C.size1, Cref->matrix0);

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


  /***********************************************************/
  /* C = alpha*A*B + beta*C, double* storage, non square     */
  /***********************************************************/

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

  dense_gemm_by_hand(alpha, M1->matrix0, M3->matrix0, M1->size0, M3->size1, M1->size1, beta,  C2ref->matrix0);
  err = dense_comparison(C2.matrix0, C2.size0, C2.size1, C2ref->matrix0);
  
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
  
  /***********************************************************/
  /* C = alpha*A*B + beta*C, SBM storage,  square            */
  /***********************************************************/
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

  dense_gemm_by_hand(alpha, M1->matrix0, M1->matrix0, M1->size0, M1->size1, M1->size0, beta,  Cref->matrix0);

  SBM_gemm_without_allocation(alpha,M2->matrix1,M2->matrix1,beta, C3.matrix1);
  DEBUG_EXPR(NM_display(&C3));
  DEBUG_EXPR(NM_dense_display(Cref->matrix0,M2->size0,M2->size1,M2->size0));
  info = SBM_dense_equal(C3.matrix1, Cref->matrix0, tol);


  
  /*  Check if it is correct */

  /* C3 and CRef must have the same values.*/

  info = SBM_dense_equal(C3.matrix1, Cref->matrix0, tol);

  if (info == 0)
    printf("Step 2 ( C = alpha*A*B + beta*C, SBM storage) ok ...\n");
  else
  {
    printf("Step 2 ( C = alpha*A*B + beta*C, SBM storage) failed ...\n");
    goto exit_3;
  }
  
  /***********************************************************/
  /* C = alpha*A*B + beta*C, SBM storage,  non square        */
  /***********************************************************/
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

  SBM_gemm_without_allocation(alpha,M2->matrix1,M4->matrix1,beta, C4.matrix1);
  dense_gemm_by_hand(alpha, M1->matrix0, M3->matrix0, M1->size0, M3->size1, M1->size1, beta,  C2ref->matrix0);
  
  DEBUG_EXPR(NM_display(&C4));
  DEBUG_EXPR(NM_dense_display(C2ref->matrix0,M2->size0,M4->size1,M2->size0));
  
  info = SBM_dense_equal(C4.matrix1, C2ref->matrix0, tol);

  if (info == 0)
    printf("Step 3 ( C = alpha*A*B + beta*C, SBM storage, non square) ok ...\n");
  else
  {
    printf("Step 3 ( C = alpha*A*B + beta*C, SBM storage, non square) failed ...\n");
    goto exit_4;
  }

  /***********************************************************/
  /* C = alpha*A*B + beta*C, double* storage, square matrix, empty column of blocks  */
  /***********************************************************/

  NumericsMatrix * M9 = test_matrix_9();


  NumericsMatrix * C7 = NM_create(NM_DENSE, M9->size0, M9->size1);
  MSAN_INIT_VAR(C7->matrix0, C7->size0 * C7->size1);

  add_initial_value_square_2(C7);

  NM_gemm(alpha, M9, M9, beta, C7);

  NumericsMatrix * C3ref= NM_create(NM_DENSE,C7->size0, C7->size1);
  add_initial_value_square_2(C3ref);
  /* gemm by hand */
  dense_gemm_by_hand(alpha, M9->matrix0, M9->matrix0, M9->size0, M9->size1, M9->size0, beta,  C3ref->matrix0);
  DEBUG_EXPR(NM_display(C7));
  DEBUG_EXPR(NM_dense_display(C3ref->matrix0,M9->size0,M9->size1,M9->size0));
  err = dense_comparison(C7->matrix0, C7->size0, C7->size1, C3ref->matrix0);

  if (err < tol)
  {
    info = 0;
  }

  if (info == 0)
    printf("Step 6 ( C = alpha*A*B + beta*C, double* storage, square matrix, empty column of blocks ) ok ...\n");
  else
  {
    printf("Step 6 ( C = alpha*A*B + beta*C, double* storage, square matrix, empty column of blocks) failed ...\n");
    goto exit_7;
  }


  /* /\**********************************************************************\/ */
  /* /\* C = alpha*A*B + beta*C, SBM storage, empty column of blocks        *\/ */
  /* /\**********************************************************************\/ */

  NumericsMatrix * M10 = test_matrix_10();
  DEBUG_EXPR(NM_display(M10););

  NumericsMatrix * C8 = NM_create(NM_SPARSE_BLOCK,M10->size0, M10->size1);
  /* This step is not necessary for NM_gemm but conveniently creates a zero matrix with the right structure */
  C8->matrix1 = SBM_zero_matrix_for_multiply(M10->matrix1, M10->matrix1);

  add_initial_value_square_2(C8);

  NM_gemm(alpha, M10, M10, beta, C8);

 
  SBM_gemm_without_allocation(alpha, M10->matrix1, M10->matrix1, beta, C8->matrix1);
  dense_gemm_by_hand(alpha, M9->matrix0, M9->matrix0, M9->size0, M9->size1, M9->size0, beta,  C3ref->matrix0);

  DEBUG_EXPR(NM_display(C8));
  DEBUG_EXPR(NM_dense_display(C3ref->matrix0,M10->size0,M10->size1,M10->size0));
 


  info = NM_dense_equal(C8,C3ref->matrix0,tol);

  if (info == 0)
    printf("Step 7 ( C = alpha*A*B + beta*C, NM_SPARSE_BLOCK storage,  empty column of blocks) ok ...\n");
  else
  {
    printf("Step 7 ( C = alpha*A*B + beta*C, NM_SPARSE_BLOCK storage,  empty column of blocks) failed ...\n");
    goto exit_8;
  }



  /* /\************************************************************************************\/ */
  /* /\* C = alpha*A*B + beta*C, SBM storage, empty column of blocks, extra blocks       *\/ */
  /* /\************************************************************************************\/ */

  NumericsMatrix * C20 = test_matrix_20();
  DEBUG_EXPR(NM_display(C20););

  add_initial_value_square_2(C20);

  NM_gemm(alpha, M10, M10, beta, C20);

  SBM_gemm_without_allocation(alpha, M10->matrix1, M10->matrix1, beta, C20->matrix1);
 
  DEBUG_EXPR(NM_display(C20));
  DEBUG_EXPR(NM_dense_display(C3ref->matrix0,M10->size0,M10->size1,M10->size0));

  info = NM_dense_equal(C20,C3ref->matrix0,tol);

  if (info == 0)
    printf("Step 8 ( C = alpha*A*B + beta*C, NM_SPARSE_BLOCK storage,  empty column of blocks, extra blocks) ok ...\n");
  else
  {
    printf("Step 8 ( C = alpha*A*B + beta*C, NM_SPARSE_BLOCK storage,  empty column of blocks, extra blocks) failed ...\n");
    goto exit_9;
  }


exit_9:
  NM_free(C20);
exit_8:
  NM_free(M10);
  NM_free(C8);
exit_7:
  NM_free(M9);
  NM_free(C7);
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

  printf("========= Starts Numerics tests for SBM_gemm_without_allocation ========= \n");

  int i, nmm = 4 ;
  NumericsMatrix ** NMM = (NumericsMatrix **)malloc(nmm * sizeof(NumericsMatrix *)) ;

  int info = test_build_first_4_NM(NMM);
  if (info != 0)
  {
    printf("Construction failed ...\n");
    return info;
  }
  printf("Construction ok ...\n");

  info = SBM_gemm_without_allocation_test(NMM,1.0,0.0);
  if (info != 0)
  {
    printf("End of Numerics tests for SBM_gemm_without_allocation : unsucessfull\n");
    return info;
  }
  info = SBM_gemm_without_allocation_test(NMM,1.0,1.0);
  if (info != 0)
  {
    printf("End of Numerics tests for SBM_gemm_without_allocation : unsucessfull\n");
    return info;
  }
  info = SBM_gemm_without_allocation_test(NMM,0.0,1.0);
  if (info != 0)
  {
    printf("End of Numerics tests for SBM_gemm_without_allocation : unsucessfull\n");
    return info;
  }
  info = SBM_gemm_without_allocation_test(NMM,0.5,0.5);
  if (info != 0)
  {
    printf("End of Numerics tests for SBM_gemm_without_allocation : unsucessfull\n");
    return info;
  }

  printf("End of Numerics tests for SBM_gemm_without_allocation : Sucessfull ...\n");
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
