
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
#include "SparseBlockMatrix.h"
#include "SiconosLapack.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sanitizer.h"

#include "NumericsSparseMatrix.h"

#include "CSparseMatrix.h"
#include "CSparseMatrix.h"
#include "numerics_verbose.h"

/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"


#include "SBM_test.h"


static int SBM_add_test1(double tol, double alpha, double beta)
{

  int info=0;
  FILE *file = fopen("data/SBM1.dat", "r");
  SparseBlockStructuredMatrix * M = SBM_new_from_file(file);
  fclose(file);
  /* SBM_print(M); */



  SparseBlockStructuredMatrix * C = SBM_add(M,M,alpha,beta);
  /* SBM_print(C); */

  int nm = M->blocksize0[M->blocknumber0-1] * M->blocksize1[M->blocknumber1-1];
  double * M_dense = (double *) malloc(nm*sizeof(double));

  SBM_to_dense(M, M_dense);

  double * C_dense = (double *) malloc(nm*sizeof(double));

  cblas_dscal(nm, 0.0, C_dense, 1);
  cblas_daxpy(nm, alpha, M_dense, 1, C_dense, 1);
  cblas_daxpy(nm, beta, M_dense, 1, C_dense, 1);

  info = SBM_dense_equal(C, C_dense, tol);

  free(M_dense);
  free(C_dense);

  SBM_clear(M);
  SBM_clear(C);
  return info;

}

static int SBM_add_test2(double tol, double alpha, double beta)
{
  printf("========= Starts SBM tests SBM_add_test2 for alpha = %e and beta = %e ========= \n", alpha, beta);
  int info = 0;
  NumericsMatrix *M2 = test_matrix_2();
  SparseBlockStructuredMatrix * SBM2= M2->matrix1;
  DEBUG_EXPR(SBM_print(SBM2););

  NumericsMatrix *M10 = test_matrix_10();
  SparseBlockStructuredMatrix * SBM10= M10->matrix1;
  DEBUG_EXPR(SBM_print(SBM10););


  SparseBlockStructuredMatrix * C2 = SBM_add(SBM2,SBM10,alpha,beta);
  DEBUG_EXPR(SBM_print(C2););

  SparseBlockStructuredMatrix * C3 = SBM_add(SBM10,SBM2,alpha,beta);
  DEBUG_EXPR(SBM_print(C3););

  int n =  M2->size0 ;
  int m =  M2->size1; 
  int nm = n*m;
  double * M2_dense = (double *) malloc(nm*sizeof(double));
  double * M10_dense = (double *) malloc(nm*sizeof(double));
  SBM_to_dense(SBM2, M2_dense);
  SBM_to_dense(SBM10, M10_dense);

  double * C2_dense = (double *) malloc(nm*sizeof(double));

  cblas_dscal(nm, 0.0, C2_dense, 1);
  
  cblas_daxpy(nm, alpha, M2_dense, 1, C2_dense, 1);

  cblas_daxpy(nm, beta, M10_dense, 1, C2_dense, 1);

  DEBUG_EXPR(NM_dense_display(C2_dense,n,m,n ));
  info = SBM_dense_equal(C2, C2_dense, tol);
  if (info == 1)
    return info;

  double * C3_dense = (double *) malloc(nm*sizeof(double));

  cblas_dscal(nm, 0.0, C3_dense, 1);
  cblas_daxpy(nm, alpha, M10_dense, 1, C3_dense, 1);
  cblas_daxpy(nm, beta, M2_dense, 1, C3_dense, 1);
  DEBUG_EXPR(NM_dense_display(C3_dense,n,m,n ));
  
  info = SBM_dense_equal(C3, C3_dense, tol);
  if (info == 1)
    return info;

  SBM_add_without_allocation(SBM2,SBM10,alpha,beta,C2,0.0);
  DEBUG_EXPR(SBM_print(C2););
  info = SBM_dense_equal(C2, C2_dense, tol);
  if (info == 1)
    return info;
  
  SBM_add_without_allocation(SBM10,SBM2,alpha,beta,C3,0.0);
  DEBUG_EXPR(SBM_print(C3););
  info = SBM_dense_equal(C3, C3_dense, tol);

  if (info == 1)
    return info;
  
  NM_clear(M2);
  NM_clear(M10);
  SBM_clear(C2);
  SBM_clear(C3);
  free(C2_dense);
  free(C3_dense);
  free(M2_dense);
  free(M10_dense);

  return info;


}


int SBM_add_test_all(void)
{
  int info =0;
  double tol = 1e-14;

  printf("========= Starts SBM tests SBM_add  ========= \n");

  info = SBM_add_test1(tol, 1.0, 1.0);
  if (info == 1)
  {
    printf("========= Ends SBM tests SBM_add  :  Unsuccessfull ========= \n");
    return info;
  }
  info = SBM_add_test1(tol, 0.0, 1.0);
  if (info == 1)
  {
    printf("========= Ends SBM tests SBM_add  :  Unsuccessfull ========= \n");
    return info;
  }
  info = SBM_add_test1(tol, 1.0, 0.0);
  if (info == 1)
  {
    printf("========= Ends SBM tests SBM_add  :  Unsuccessfull ========= \n");
    return info;
  }
  info = SBM_add_test1(tol, 0.5, 0.5);
  if (info == 1)
  {
    printf("========= Ends SBM tests SBM_add  :  Unsuccessfull ========= \n");
    return info;
  }

  info = SBM_add_test2(tol, 1.0, 1.0);
  if (info == 1)
  {
    printf("========= Ends SBM tests SBM_add  :  Unsuccessfull ========= \n");
    return info;
  }
  info = SBM_add_test2(tol, 0.0, 1.0);
  if (info == 1)
  {
    printf("========= Ends SBM tests SBM_add  :  Unsuccessfull ========= \n");
    return info;
  }
  info = SBM_add_test2(tol, 1.0, 0.0);
  if (info == 1)
  {
    printf("========= Ends SBM tests SBM_add  :  Unsuccessfull ========= \n");
    return info;
  }
  info = SBM_add_test2(tol, 0.5, 0.5);
  if (info == 1)
  {
    printf("========= Ends SBM tests SBM_add  :  Unsuccessfull ========= \n");
    return info;
  }
  printf("========= Ends SBM tests SBM_add  :  successfull ========= \n");


return info;
}

int test_SBM_column_permutation_all(void)
{

  printf("========= Starts SBM tests 3 for SBM ========= \n");
  
  FILE *file = fopen("data/SBM1.dat", "r");
  SparseBlockStructuredMatrix * M = SBM_new_from_file(file);
  fclose(file);
  /*alloc enough memory */
  int res = test_SBM_column_permutation(M);
  if (res)
  {
    printf("========= Failed SBM tests 3 for SBM  ========= \n");
    return 1;
  }
  SBM_clear(M);
  file = fopen("data/SBM2.dat", "r");
  M = SBM_new_from_file(file);
  fclose(file);
  res = test_SBM_column_permutation(M);
  if (res)
  {
    printf("========= Failed SBM tests 3 for SBM  ========= \n");
    return 1;
  }
  SBM_clear(M);
  printf("\n========= Succed SBM tests 3 for SBM  ========= \n");
  return 0;

}


int SBM_extract_component_3x3_all(void)
{

  printf("========= Starts SBM tests 7 for SBM ========= \n");
  FILE *file = fopen("data/SBM1.dat", "r");
  SparseBlockStructuredMatrix * M = SBM_new_from_file(file);
  fclose(file);
  SBM_print(M);
  
  SparseBlockStructuredMatrix * N = SBM_new();
  unsigned int row_components[1] = {0};
  unsigned int row_components_size =1;
  unsigned int col_components[1] = {0};
  unsigned int col_components_size =1;
  SBM_extract_component_3x3(M, N, row_components, row_components_size, col_components, col_components_size   );
  SBM_print(N);

  SparseBlockStructuredMatrix * T = SBM_new();
  unsigned int row_components_T[2] = {1,2};
  unsigned int row_components_size_T =2;
  unsigned int col_components_T[2] = {1,2};
  unsigned int col_components_size_T =2;
  SBM_extract_component_3x3(M, T, row_components_T, row_components_size_T, col_components_T, col_components_size_T   );
  SBM_print(T);

  SparseBlockStructuredMatrix * NT = SBM_new();
  unsigned int row_components_NT[2] = {0};
  unsigned int row_components_size_NT =1;
  
  unsigned int col_components_NT[2] = {1,2};
  unsigned int col_components_size_NT =2;
  SBM_extract_component_3x3(M, NT, row_components_NT, row_components_size_NT, col_components_NT, col_components_size_NT   );
  SBM_print(NT);
  
  SparseBlockStructuredMatrix * TN = SBM_new();
  unsigned int row_components_TN[2] = {1,2};
  unsigned int row_components_size_TN =2;
  
  unsigned int col_components_TN[2] = {0};
  unsigned int col_components_size_TN =1;
  SBM_extract_component_3x3(M, TN, row_components_TN, row_components_size_TN, col_components_TN, col_components_size_TN   );
  SBM_print(TN);
  
  
  
  int res = test_SBM_row_to_dense(M);
  if (res)
  {
    printf("========= Failed SBM tests 7 for SBM  ========= \n");
    return 1;
  }

  SBM_clear(M);
  
  return 0;

}


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
  NM_clear(C20);
exit_8:
  NM_clear(M10);
  NM_clear(C8);
exit_7:
  NM_clear(M9);
  NM_clear(C7);
exit_4:
  NM_clear(&C4);
exit_3:
  NM_clear(&C3);
exit_2:
  free(C2.matrix0);
  NM_clear(C2ref);
exit_1:
  NM_clear(Cref);
  free(C.matrix0);
  return info;
}

int SBM_gemm_without_allocation_all(void)
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
    NM_clear(NMM[i]);
    free(NMM[i]);
  }

  free(NMM);



  printf("========= End Numerics tests for NumericsMatrix ========= \n");
  return info;
}


static int SBM_multiply_test1(double tol)
{

  int info=0;

  FILE *file = fopen("data/SBM1.dat", "r");
  SparseBlockStructuredMatrix * M = SBM_new_from_file(file);
  fclose(file);
  DEBUG_EXPR(SBM_print(M););

  SparseBlockStructuredMatrix * C = SBM_multiply(M,M);
  DEBUG_EXPR(SBM_print(C););
  
  int n = M->blocksize0[M->blocknumber0-1];
  int m = M->blocksize1[M->blocknumber1-1];
  int nm = n*m;
  
  double * M_dense = (double *) malloc(nm*sizeof(double));
  SBM_to_dense(M, M_dense);

  double * C_dense = (double *) malloc(nm*sizeof(double));

  double beta = 0.0;
  double alpha=1.0;
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, m, n,
              alpha, M_dense, n, M_dense, n, beta, C_dense, n);

  DEBUG_EXPR(NM_dense_display(C_dense,n,m,n););

  info = SBM_dense_equal(C, C_dense, tol);
  
  free(M_dense);
  free(C_dense);

  SBM_clear(M);
  SBM_clear(C);
  return info;
}

static int SBM_multiply_test2(double tol)
{
  printf("========= Starts SBM tests SBM_multiply_test2 ========= \n");
  int info = 0;
  NumericsMatrix *M2 = test_matrix_2();
  SparseBlockStructuredMatrix * SBM2= M2->matrix1;
  DEBUG_EXPR(SBM_print(SBM2););

  NumericsMatrix *M10 = test_matrix_10();
  SparseBlockStructuredMatrix * SBM10= M10->matrix1;
  DEBUG_EXPR(SBM_print(SBM10););


  SparseBlockStructuredMatrix * C2 = SBM_multiply(SBM2,SBM10);
  DEBUG_EXPR(SBM_print(C2););

  SparseBlockStructuredMatrix * C3 = SBM_multiply(SBM10,SBM2);
  DEBUG_EXPR(SBM_print(C3););

  int n =  M2->size0 ;
  int m =  M2->size1; 
  int nm = n*m;
  double * M2_dense = (double *) malloc(nm*sizeof(double));
  double * M10_dense = (double *) malloc(nm*sizeof(double));
  SBM_to_dense(SBM2, M2_dense);
  SBM_to_dense(SBM10, M10_dense);

  double * C2_dense = (double *) malloc(nm*sizeof(double));

  double beta  =0.0;
  double alpha =1.0;
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, m, n,
              alpha, M2_dense, n, M10_dense, n, beta, C2_dense, n);

  DEBUG_EXPR(NM_dense_display(C2_dense,n,m,n););

  info = SBM_dense_equal(C2, C2_dense, tol);
  if (info == 1)
    return info;

  double * C3_dense = (double *) malloc(nm*sizeof(double));
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, m, n,
              alpha, M10_dense, n, M2_dense, n, beta, C3_dense, n);

  DEBUG_EXPR(NM_dense_display(C3_dense,n,m,n););


  DEBUG_EXPR(NM_dense_display(C3_dense,n,m,n ));
  
  info = SBM_dense_equal(C3, C3_dense, tol);
  if (info == 1)
    return info;

  NM_clear(M2);
  NM_clear(M10);
  SBM_clear(C2);
  SBM_clear(C3);
  free(C2_dense);
  free(C3_dense);
  free(M2_dense);
  free(M10_dense);

  return info;


}
static int SBM_multiply_test3(double tol)
{
  printf("========= Starts SBM tests SBM_multiply_test3  ========= \n");
  int info = 0;
  NumericsMatrix *M2 = test_matrix_2();
  SparseBlockStructuredMatrix * SBM2= M2->matrix1;
  DEBUG_EXPR(SBM_print(SBM2););

  NumericsMatrix *M4 = test_matrix_4();
  SparseBlockStructuredMatrix * SBM4= M4->matrix1;
  DEBUG_EXPR(SBM_print(SBM4););


  SparseBlockStructuredMatrix * C2 = SBM_multiply(SBM2,SBM4);
  DEBUG_EXPR(SBM_print(C2););


  int n =  M2->size0 ;

  int m =  M4->size1 ;
  
  
  double * M2_dense = (double *) malloc(n*n*sizeof(double));
  double * M4_dense = (double *) malloc(n*m*sizeof(double));
  SBM_to_dense(SBM2, M2_dense);
  SBM_to_dense(SBM4, M4_dense);


  double * C2_dense = (double *) malloc(n*m*sizeof(double));

  double beta  =0.0;
  double alpha =1.0;
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, m, n,
              alpha, M2_dense, n, M4_dense, n, beta, C2_dense, n);

  DEBUG_EXPR(NM_dense_display(C2_dense,n,m,n););

  info = SBM_dense_equal(C2, C2_dense, tol);
  

  NM_clear(M2);
  NM_clear(M4);
  SBM_clear(C2);
  free(C2_dense);
  free(M2_dense);
  free(M4_dense);

  return info;


}
int SBM_multiply_test_all(void)
{

  int info =0;
  double tol = 1e-14;

  printf("========= Starts SBM tests SBM_multiply ========= \n");

  info = SBM_multiply_test1(tol);
  if (info == 1)
  {
    printf("========= Ends SBM tests SBM_multiply  :  Unsuccessfull ========= \n");
    return info;
  }
  info = SBM_multiply_test2(tol);
  if (info == 1)
  {
    printf("========= Ends SBM tests SBM_multiply  :  Unsuccessfull ========= \n");
    return info;
  }
  info = SBM_multiply_test3(tol);
  if (info == 1)
  {
    printf("========= Ends SBM tests SBM_multiply  :  Unsuccessfull ========= \n");
    return info;
  }
 
  printf("========= Ends SBM tests SBM_multiply  :  successfull ========= \n");
  return info;

}


int test_SBM_row_permutation_all(void)
{

  printf("========= Starts SBM tests 2 for SBM ========= \n");
  FILE *file = fopen("data/SBM1.dat", "r");
  SparseBlockStructuredMatrix * M = SBM_new_from_file(file);
  fclose(file);
  /*alloc enough memory */
  int res = test_SBM_row_permutation(M);
  if (res)
  {
    printf("========= Failed SBM tests 2 for SBM  ========= \n");
    return 1;
  }

  SBM_clear(M);

  
  file = fopen("data/SBM2.dat", "r");
  SparseBlockStructuredMatrix * M2 = SBM_new_from_file(file);
  fclose(file);
  res = test_SBM_row_permutation(M2);
  if (res)
  {
    printf("========= Failed SBM tests 2 for SBM  ========= \n");
    return 1;
  }
  SBM_clear(M2);
  printf("\n========= Succed SBM tests 2 for SBM  ========= \n");
  return 0;




}


int test_SBM_row_to_dense_all(void)
{

  printf("========= Starts SBM tests 1 for SBM ========= \n");
  
  FILE *file = fopen("data/SBM1.dat", "r");
  SparseBlockStructuredMatrix * M =SBM_new_from_file(file);
  fclose(file);
  /*alloc enough memory */
  int res = test_SBM_row_to_dense(M);
  if (res)
  {
    printf("========= Failed SBM tests 1 for SBM  ========= \n");
    return 1;
  }

  SBM_clear(M);



  file = fopen("data/SBM2.dat", "r");
  SparseBlockStructuredMatrix * M2 = SBM_new_from_file(file);
  fclose(file);
  res = test_SBM_row_to_dense(M2);
  if (res)
  {
    printf("========= Failed SBM tests 1 for SBM  ========= \n");
    return 1;
  }
  SBM_clear(M2);
  printf("\n========= Succeeded SBM tests 1 for SBM  ========= \n");
  return 0;

}



int SBM_to_dense_all(void)
{
  int res;
  printf("========= Starts SBM tests 4 for SBM ========= \n");
  
  FILE *file = fopen("data/SBM2.dat", "r");
  SparseBlockStructuredMatrix * M =SBM_new_from_file(file);
  SBM_print(M);
  fclose(file);
  /*alloc enough memory */
  CSparseMatrix sparseMat;
  res = SBM_to_sparse_init_memory(M, &sparseMat);
  if (res)
  {
    printf("========= Failed SBM tests 4 for SBM  ========= \n");
    return 1;
  }

  res = SBM_to_sparse(M, &sparseMat);
  if (res)
  {
    printf("========= Failed SBM tests 4 for SBM  ========= \n");
    return 1;
  }
  cs_print(&sparseMat, 1);
  CSparseMatrix_spfree_on_stack(&sparseMat);

  int n = M->blocksize0[M->blocknumber0 - 1];
  int m = M->blocksize1[M->blocknumber1 - 1];
  double * denseMat = (double *)malloc(n * m * sizeof(double));
  SBM_to_dense(M, denseMat);
  if (res)
  {
    printf("========= Failed SBM tests 4 for SBM  ========= \n");
    return 1;
  }
  printf("[");
  for (int i = 0; i < n * m; i++)
  {
    printf("%lf ", denseMat[i]);
    if ((i + 1) % m == 0)
      printf("\n");
  }
  printf("]");
  printf("\n (warning: column-major) \n");

  free(denseMat);
  SBM_clear(M);
  printf("\n========= Succed SBM tests 4 for SBM  ========= \n");
  return 0;

}



int SBM_to_sparse_all(void)
{
  int res;
  printf("========= Starts SBM tests 4 for SBM ========= \n");
  FILE *file = fopen("data/SBM2.dat", "r");
  SparseBlockStructuredMatrix * M = SBM_new_from_file(file);
  SBM_print(M);
  fclose(file);
  /*alloc enough memory */
  CSparseMatrix sparseMat;
  res = SBM_to_sparse_init_memory(M, &sparseMat);
  if (res)
  {
    printf("========= Failed SBM tests 4 for SBM  ========= \n");
    return 1;
  }

  res = SBM_to_sparse(M, &sparseMat);
  if (res)
  {
    printf("========= Failed SBM tests 4 for SBM  ========= \n");
    return 1;
  }
  cs_print(&sparseMat, 1);
  CSparseMatrix_spfree_on_stack(&sparseMat);

  int n = M->blocksize0[M->blocknumber0 - 1];
  int m = M->blocksize1[M->blocknumber1 - 1];
  double * denseMat = (double *)malloc(n * m * sizeof(double));
  SBM_to_dense(M, denseMat);
  if (res)
  {
    printf("========= Failed SBM tests 4 for SBM  ========= \n");
    return 1;
  }
  printf("[");
  for (int i = 0; i < n * m; i++)
  {
    printf("%lf ", denseMat[i]);
    if ((i + 1) % m == 0)
      printf("\n");
  }
  printf("]");
  printf("\n (warning: column-major) \n");

  free(denseMat);
  printf("NUMERICS_SBM_FREE_BLOCK value %d", NUMERICS_SBM_FREE_BLOCK);
  SBM_clear(M);
  printf("\n========= Succed SBM tests 4 for SBM  ========= \n");
  return 0;

}

static void add_initial_value_square_1a(NumericsMatrix * M)
{
  int i=0, j=0;
  for (i=0; i < 4 ; i++)
  {
    for (j=0; j < 4 ; j++)
      NM_zentry(M,i,j,1.0);
  }
  for (i=0; i < 4 ; i++)
  {
    for (j=4; j < 6 ; j++)
      NM_zentry(M,i,j,2.0);
  }
  for (i=4; i < 6 ; i++)
  {
    for (j=4; j < 6 ; j++)
      NM_zentry(M,i,j,3.0);
  }
  for (i=4; i < 6 ; i++)
  {
    for (j=6; j < 8 ; j++)
      NM_zentry(M,i,j,4.0);
  }
  for (i=6; i < 8 ; i++)
  {
    for (j=0; j < 4 ; j++)
      NM_zentry(M,i,j,5.0);
  }
  for (i=6; i < 8 ; i++)
  {
    for (j=6; j < 8 ; j++)
      NM_zentry(M,i,j,6.0);
  }
}
static void add_initial_value_square_SBM_1(SparseBlockStructuredMatrix * M)
{
  int i=0, j=0;
  for (i=0; i < 4 ; i++)
  {
    for (j=0; j < 4 ; j++)
      CHECK_RETURN(SBM_zentry(M,i,j,1.0));
  }
  for (i=0; i < 4 ; i++)
  {
    for (j=4; j < 6 ; j++)
      CHECK_RETURN(SBM_zentry(M,i,j,2.0));
  }
  for (i=4; i < 6 ; i++)
  {
    for (j=4; j < 6 ; j++)
      SBM_zentry(M,i,j,3.0);
  }
  for (i=4; i < 6 ; i++)
  {
    for (j=6; j < 8 ; j++)
      SBM_zentry(M,i,j,4.0);
  }
  for (i=6; i < 8 ; i++)
  {
    for (j=0; j < 4 ; j++)
      SBM_zentry(M,i,j,5.0);
  }
  for (i=6; i < 8 ; i++)
  {
    for (j=6; j < 8 ; j++)
      SBM_zentry(M,i,j,6.0);
  }
}


static int SBM_zentry_test1(double tol)
{

  int info=0;

  NumericsMatrix * M2 = test_matrix_2();

  NumericsMatrix * C= NM_create(NM_DENSE, M2->size0, M2->size1);
  add_initial_value_square_1a(C);


  NumericsMatrix * C2= NM_create(NM_SPARSE_BLOCK, M2->size0, M2->size1);
  C2->matrix1 = SBM_zero_matrix_for_multiply(M2->matrix1, M2->matrix1);
  add_initial_value_square_SBM_1(C2->matrix1);
  
  DEBUG_EXPR(NM_display(C));
  DEBUG_EXPR(NM_display(C2));
  info = NM_dense_equal(C2,C->matrix0,tol);
  NM_clear(M2);
  NM_clear(C2);
  NM_clear(C);
  return info;
}

static int SBM_zentry_test2(double tol)
{

  int info=0;

  NumericsMatrix * M2 = test_matrix_2();
  
  /* CHECK_RETURN(SBM_zentry(M2->matrix1,0,8,1.0)); */
  /* CHECK_RETURN(SBM_zentry(M2->matrix1,8,0,1.0)); */
  /* CHECK_RETURN(SBM_zentry(M2->matrix1,0,7,1.0)); */
  
  info = SBM_zentry(M2->matrix1,0,8,1.0);
  info = SBM_zentry(M2->matrix1,8,0,1.0);
  info = SBM_zentry(M2->matrix1,0,7,1.0);

  NM_clear(M2);

  return info;
}


int SBM_zentry_all(void)
{

  int info =0;
  double tol = 1e-14;

  printf("========= Starts SBM tests SBM_zentry  ========= \n");

  info = SBM_zentry_test1(tol);
  if (info == 1)
  {
    printf("========= Ends SBM tests SBM_zentry  :  Unsuccessfull ========= \n");
    return info;
  }
  info = SBM_zentry_test2(tol);
  if (info == 1)
  {
    printf("========= Ends SBM tests SBM_zentry  :  Unsuccessfull ========= \n");
    return info;
  }
  printf("========= Ends SBM tests SBM_zentry  :  successfull ========= \n");

  return info;

}


int main()
{

  int info =  SBM_add_test_all();

  info += SBM_zentry_all();
  
  info += SBM_to_dense_all();

  info += SBM_to_sparse_all();

  info += SBM_multiply_test_all();
  
  info += SBM_gemm_without_allocation_all();

  info += test_SBM_row_to_dense_all();

  info += test_SBM_column_permutation_all();

  info += test_SBM_row_permutation_all();

  info += SBM_extract_component_3x3_all();

  return info;
}
