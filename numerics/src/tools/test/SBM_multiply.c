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
/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"


static int SBM_multiply_test1(double tol)
{

  int info=0;
  SparseBlockStructuredMatrix M;
  FILE *file = fopen("data/SBM1.dat", "r");
  SBM_new_from_file(&M, file);
  fclose(file);
  DEBUG_EXPR(SBM_print(&M););

  SparseBlockStructuredMatrix * C = SBM_multiply(&M,&M);
  DEBUG_EXPR(SBM_print(C););
  
  int n = M.blocksize0[M.blocknumber0-1];
  int m = M.blocksize1[M.blocknumber1-1];
  int nm = n*m;
  
  double * M_dense = (double *) malloc(nm*sizeof(double));
  SBM_to_dense(&M, M_dense);

  double * C_dense = (double *) malloc(nm*sizeof(double));

  double beta = 0.0;
  double alpha=1.0;
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, m, n,
              alpha, M_dense, n, M_dense, n, beta, C_dense, n);

  DEBUG_EXPR(NM_dense_display(C_dense,n,m,n););

  info = SBM_dense_equal(C, C_dense, tol);
  
  free(M_dense);
  free(C_dense);

  SBMfree(&M, NUMERICS_SBM_FREE_BLOCK);
  SBM_free(C);
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

  NM_free(M2);
  NM_free(M10);
  SBM_free(C2);
  SBM_free(C3);
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
  

  NM_free(M2);
  NM_free(M4);
  SBM_free(C2);
  free(C2_dense);
  free(M2_dense);
  free(M4_dense);

  return info;


}
int main(void)
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
