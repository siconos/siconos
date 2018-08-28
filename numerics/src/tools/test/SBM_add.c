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


static int SBM_add_test1(double tol, double alpha, double beta)
{

  int info=0;
  SparseBlockStructuredMatrix M;
  FILE *file = fopen("data/SBM1.dat", "r");
  SBM_new_from_file(&M, file);
  fclose(file);
  /* SBM_print(&M); */



  SparseBlockStructuredMatrix * C = SBM_add(&M,&M,alpha,beta);
  /* SBM_print(C); */

  int nm = M.blocksize0[M.blocknumber0-1] * M.blocksize1[M.blocknumber1-1];
  double * M_dense = (double *) malloc(nm*sizeof(double));

  SBM_to_dense(&M, M_dense);

  double * C_dense = (double *) malloc(nm*sizeof(double));

  cblas_dscal(nm, 0.0, C_dense, 1);
  cblas_daxpy(nm, alpha, M_dense, 1, C_dense, 1);
  cblas_daxpy(nm, beta, M_dense, 1, C_dense, 1);

  info = SBM_dense_equal(C, C_dense, tol);

  free(M_dense);
  free(C_dense);

  SBMfree(&M, NUMERICS_SBM_FREE_BLOCK);
  SBM_free(C);
  return info;

}

static int SBM_add_test2(double tol, double alpha, double beta)
{
  printf("========= Starts SBM tests SBM_add_test2 for alpha = %e and beta = %e ========= \n", alpha, beta);
  int info = 0;
  NumericsMatrix *M2 = test_matrix_2();
  SparseBlockStructuredMatrix * SBM2= M2->matrix1;
  DEBUG_EXPR(SBM_print(SBM2););

  NumericsMatrix *M5 = test_matrix_5();
  SparseBlockStructuredMatrix * SBM5= M5->matrix1;
  DEBUG_EXPR(SBM_print(SBM5););


  SparseBlockStructuredMatrix * C2 = SBM_add(SBM2,SBM5,alpha,beta);
  DEBUG_EXPR(SBM_print(C2););

  SparseBlockStructuredMatrix * C3 = SBM_add(SBM5,SBM2,alpha,beta);
  DEBUG_EXPR(SBM_print(C3););

  int n =  M2->size0 ;
  int m =  M2->size1; 
  int nm = n*m;
  double * M2_dense = (double *) malloc(nm*sizeof(double));
  double * M5_dense = (double *) malloc(nm*sizeof(double));
  SBM_to_dense(SBM2, M2_dense);
  SBM_to_dense(SBM5, M5_dense);

  double * C2_dense = (double *) malloc(nm*sizeof(double));

  cblas_dscal(nm, 0.0, C2_dense, 1);
  
  cblas_daxpy(nm, alpha, M2_dense, 1, C2_dense, 1);

  cblas_daxpy(nm, beta, M5_dense, 1, C2_dense, 1);

  DEBUG_EXPR(NM_dense_display(C2_dense,n,m,n ));
  info = SBM_dense_equal(C2, C2_dense, tol);
  if (info == 1)
    return info;

  double * C3_dense = (double *) malloc(nm*sizeof(double));

  cblas_dscal(nm, 0.0, C3_dense, 1);
  cblas_daxpy(nm, alpha, M5_dense, 1, C3_dense, 1);
  cblas_daxpy(nm, beta, M2_dense, 1, C3_dense, 1);
  DEBUG_EXPR(NM_dense_display(C3_dense,n,m,n ));
  
  info = SBM_dense_equal(C3, C3_dense, tol);
  if (info == 1)
    return info;

  SBM_add_without_allocation(SBM2,SBM5,alpha,beta,C2);
  info = SBM_dense_equal(C2, C2_dense, tol);
  DEBUG_EXPR(SBM_print(C2););
  if (info == 1)
    return info;
  SBM_add_without_allocation(SBM5,SBM2,alpha,beta,C3);
  info = SBM_dense_equal(C3, C3_dense, tol);
  DEBUG_EXPR(SBM_print(C3););
  if (info == 1)
    return info;
  
  NM_free(M2);
  NM_free(M5);
  SBM_free(C2);
  SBM_free(C3);
  free(C2_dense);
  free(C3_dense);
  free(M2_dense);
  free(M5_dense);

  return info;


}





int main(void)
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
