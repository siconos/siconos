/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

#include <assert.h>                      // for assert
#include <math.h>                        // for fabs
#include <stdint.h>                      // for SIZE_MAX
#include <stdio.h>                       // for printf, fclose, fopen, NULL
#include <stdlib.h>                      // for free, malloc, calloc
#include <float.h>                       // for DBL_EPSILON
#include "SiconosBlas.h"                 // for cblas_ddot, cblas_dgemv, cbl...
#include "CSparseMatrix_internal.h"               // for CS_INT, cs_print, cs
#include "NumericsFwd.h"                 // for NumericsMatrix, SparseBlockS...
#include "NumericsMatrix.h"              // for NumericsMatrix, NM_clear, NM_...
#include "NumericsSparseMatrix.h"        // for NumericsSparseMatrix, NSM_TR...
#include "NumericsVector.h"              // for NV_equal
#include "SparseBlockMatrix.h"           // for SBM_zero_matrix_for_multiply
#include "debug.h"                       // for DEBUG_EXPR, DEBUG_PRINTF
#include "numericsMatrixTestFunction.h"  // for test_build_first_4_NM, NM_de...
#include "numerics_verbose.h"            // for numerics_error
#include "sanitizer.h"                   // for MSAN_INIT_VAR


#ifdef WITH_MUMPS
#include "NM_MUMPS.h"
#endif

#ifndef SIZE_MAX
# ifdef __SIZE_MAX__
#  define SIZE_MAX __SIZE_MAX__
# else
#  define SIZE_MAX std::numeric_limits<size_t>::max()
# endif
#endif
static int NM_read_write_test(void)
{

  printf("========= Starts Numerics tests for NumericsMatrix ========= \n");

  int i, nmm = 4 ;
  NumericsMatrix ** NMM = (NumericsMatrix **)malloc(nmm * sizeof(NumericsMatrix *)) ;
  NumericsMatrix ** Mread = (NumericsMatrix **)malloc(nmm * sizeof(NumericsMatrix *)) ;

  int info = test_build_first_4_NM(NMM);

  if(info != 0)
  {
    printf("Construction failed ...\n");
    goto free;
  }
  printf("Construction ok ...\n");

  /* Test of various I/O functions */

  for(i = 0 ; i < nmm; i++)
  {

    printf("test on NMM[%i]\n", i);

    NM_display(NMM[i]);
    NM_display_row_by_row(NMM[i]);
    FILE * foutput = fopen("testprintInfile.dat", "w");
    NM_write_in_file(NMM[i], foutput);
    fclose(foutput);
    FILE * finput = fopen("testprintInfile.dat", "r");
    NM_read_in_file(NMM[i], finput);
    fclose(finput);

    FILE * finput2 = fopen("testprintInfile.dat", "r");
    Mread[i] = NM_new_from_file(finput2);
    fclose(finput2);

    char  filename[50] = "testprintInfileName.dat";
    NM_write_in_filename(NMM[i], filename);
    NM_read_in_filename(NMM[i], filename);
    printf("end of test on NMM[%i]\n", i);

  }
  for(i = 0 ; i < nmm; i++, i++)
  {
    FILE * foutput2 = fopen("testprintInfileForScilab.dat", "w");
    NM_write_in_file_scilab(NMM[i], foutput2);
    fclose(foutput2);
  }

  /* free memory */
free:
  for(i = 0 ; i < nmm; i++)
  {
    NM_clear(NMM[i]);
    free(NMM[i]);
    NM_clear(Mread[i]);
    free(Mread[i]);
  }

  free(NMM);
  free(Mread);



  printf("========= End Numerics tests for NumericsMatrix ========= \n");
  return info;
}



static int NM_add_to_diag3_test(NumericsMatrix* M, double alpha)
{
  printf("\n == Numerics tests: NM_add_to_diag3(...) == \n");
  printf("Starts NM_add_to_diag3_test for alpha = %e\n",alpha);

  int info =-1;

  /***********************************************************/
  /* C = C + alpha +I NM_SPARSE_BLOCK storage, square matrix  */
  /***********************************************************/
  int n = M->size0;
  int m = M->size1;
  NumericsMatrix * C1= NM_new();
  NM_copy(M, C1);

  NumericsMatrix * Cref = NM_create(NM_DENSE, n, m);
  SBM_to_dense(M->matrix1,Cref->matrix0);


  NM_add_to_diag3(C1, alpha);
  DEBUG_EXPR(NM_display(C1););


  double * Id = (double *) calloc(n*m, sizeof(double));
  for(int i = 0; i < n; i++)
  {
    Id[i + i  *n  ] =1.0;
  }

  cblas_daxpy(n*m, alpha, Id, 1, Cref->matrix0, 1);
  DEBUG_EXPR(NM_display(Cref););
  info = NM_dense_equal(C1, Cref->matrix0, 1e-14);

  if(info == 0)
    printf("Step 0 ( C = C + alpha*I,  NM_SPARSE_BLOCK storage, square matrix) ok ...\n");
  else
  {
    printf("Step 0 (C = C + alpha*I, NM_SPARSE_BLOCK storage, square matrix) failed ...\n");
    free(Id);
    NM_clear(Cref);
  }

  /***********************************************************/
  /* C = C + alpha +I NM_SPARSE storage, square matrix  */
  /***********************************************************/
  NumericsMatrix * C2= NM_create(NM_SPARSE, n, m);
  NM_copy_to_sparse(M,C2);
  NM_add_to_diag3(C2, alpha);

  info = NM_dense_equal(C2, Cref->matrix0, 1e-14);

  if(info == 0)
    printf("Step 1 ( C = C + alpha*I,  NM_SPARSE storage, square matrix) ok ...\n");
  else
  {
    printf("Step 1 (C = C + alpha*I, NM_SPARSE storage, square matrix) failed ...\n");
    NM_clear(C2);
  }

  NM_clear(C2);
  free(Id);
  NM_clear(Cref);

  return info;
}


static int  NM_add_to_diag3_test_all(void)
{

  printf("========= Starts Numerics tests for NM_add_to_diag3 ========= \n");

  FILE *file = fopen("data/SBM1.dat", "r");
  SparseBlockStructuredMatrix * SBM = SBM_new_from_file(file);
  fclose(file);

  NumericsMatrix * M = NM_create(NM_SPARSE_BLOCK, SBM->blocksize0[SBM->blocknumber0-1],SBM->blocksize1[SBM->blocknumber1-1]);
  M->matrix1=SBM;


  int info = NM_add_to_diag3_test(M,1.0);
  if(info != 0)
  {
    printf("End of  : Numerics tests for NM_add_to_diag3unsucessfull\n");
    return info;
  }

  printf("End of Numerics tests for NM_add_to_diag3 ...\n");
  if(info != 0) return info;
  /* free memory */


  NM_clear(M);

  printf("========= End Numerics tests for NM_add_to_diag3 ========= \n");
  return info;
}

int triplet_to_dense(void);
int csc_to_dense(void);

int triplet_to_dense(void)
{
  int info =1;
  NumericsMatrix *A;
  const char * filename =  "./data/NSM_triplet_162x162.dat";
  A = NM_create_from_filename(filename);
  /* NM_display(A); */

  NumericsMatrix *B = NM_create(NM_DENSE,A->size0,A->size1);

  info =  NM_to_dense(A, B);
  if(info != 0)
    numerics_error("triplet_to_dense", "conversion error.");

  /* NM_display(B); */

  return (int) NM_equal(A,B) -1 ;
}

int csc_to_dense(void)
{
  int info =1;
  NumericsMatrix *A;
  /* char * filename =  "./data/NSM_triplet_162x162.dat"; */
  /* A = NM_create_from_filename(filename);   */
  /* /\* NM_display(A); *\/ */
  /* CSparseMatrix* A_csc = NM_csc(A); */
  /* NumericsSparseMatrix * C_NSM =  numericsSparseMatrix_new(); */
  /* C_NSM->csc = A_csc ; */
  /* C_NSM->origin= NS_CSC; */
  /* NumericsMatrix*  A_CSC  = NM_create_from_data(NM_SPARSE, A->size0, A->size1, (void*)C_NSM); */
  /* NM_display(A_CSC); */
  /* NM_write_in_filename(A_CSC, "./data/NSM_csc_162x162.dat"); */


  const char * filename =  "./data/NSM_csc_162x162.dat";
  A = NM_new_from_filename(filename);
  /* NM_display(A); */


  NumericsMatrix *B = NM_create(NM_DENSE,A->size0,A->size1);
  info =  NM_to_dense(A, B);
  if(info != 0)
    numerics_error("csc_to_dense", "conversion error.\n");

  /* NM_display(B);   */

  return (int) NM_equal(A,B) -1 ;

}

static int to_dense_test(void)
{

  int info = triplet_to_dense();
  printf("triplet_to_dense() :  info = %i\n", info);
  info += csc_to_dense();
  printf("csc_to_dense() :  info = %i\n", info);

  return info;
}


/* ============================================================================================================================== */

static void add_initial_value_square_1(NumericsMatrix * M)
{
  int i=0, j=0;
  for(i=0; i < 4 ; i++)
  {
    for(j=0; j < 4 ; j++)
      NM_zentry(M,i,j,1.0+i+j);
  }
  for(i=0; i < 4 ; i++)
  {
    for(j=4; j < 6 ; j++)
      NM_zentry(M,i,j,2.0+i+j);
  }
  for(i=4; i < 6 ; i++)
  {
    for(j=4; j < 6 ; j++)
      NM_zentry(M,i,j,3.0+i+j);
  }
  for(i=4; i < 6 ; i++)
  {
    for(j=6; j < 8 ; j++)
      NM_zentry(M,i,j,4.0+i+j);
  }
  for(i=6; i < 8 ; i++)
  {
    for(j=0; j < 4 ; j++)
      NM_zentry(M,i,j,5.0+i+j);
  }
  for(i=6; i < 8 ; i++)
  {
    for(j=6; j < 8 ; j++)
      NM_zentry(M,i,j,6.0+i+j);
  }
}

static void add_initial_value_square_2(NumericsMatrix * M)
{
  int i=0, j=0;
  for(i=0; i < 4 ; i++)
  {
    for(j=0; j < 4 ; j++)
      NM_zentry(M,i,j,1.0+i+j);
  }
  for(i=4; i < 6 ; i++)
  {
    for(j=6; j < 8 ; j++)
      NM_zentry(M,i,j,4.0+i+j);
  }
  for(i=6; i < 8 ; i++)
  {
    for(j=0; j < 4 ; j++)
      NM_zentry(M,i,j,5.0+i+j);
  }
  for(i=6; i < 8 ; i++)
  {
    for(j=6; j < 8 ; j++)
      NM_zentry(M,i,j,6.0+i+j);
  }
}



static void add_initial_value_rectangle_1(NumericsMatrix * M)
{
  for(int i=0; i < 4 ; i++)
  {
    for(int j=0; j < 4 ; j++)
      NM_zentry(M,i,j,1.0+i+j);
  }
  for(int i=6; i < 8 ; i++)
  {
    for(int j=0; j < 4 ; j++)
      NM_zentry(M,i,j,5.0+i+j);
  }
}

/* hand made gemm of nxp A matrix and pxm B matrix */
static void dense_gemm_by_hand(double alpha, double * A, double * B, int n, int m, int p, double beta, double *C)
{
  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j < m; j++)
    {
      double sum = beta  * C[i + j * n] ;
      for(int k = 0; k < p ; k++)
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
  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j < m ; j++)
    {
      DEBUG_PRINTF("Cref[%i+%i*%i]= %lf\t\t", i, j, n, Cref[i + j * n]);
      DEBUG_PRINTF("C[%i+%i*%i]= %lf\t", i, j, n, C[i + j * n]);
      err += (Cref[i + j * n] - C[i + j * n]) * (Cref[i + j * n] - C[i + j * n]);
      DEBUG_PRINTF("err = %lf\n", err);
    }
  }
  return err;
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
  double tol = 1e-14;

  /***********************************************************/
  /* C = alpha*A*B + beta*C, double* storage, square matrix  */
  /***********************************************************/

  NumericsMatrix C;
  NM_null(&C);
  C.storageType = NM_DENSE;
  C.size0 = M1->size0;
  C.size1 = M1->size1;
  C.matrix0 = (double *)calloc(C.size0 * C.size1, sizeof(double));

  MSAN_INIT_VAR(C.matrix0, C.size0 * C.size1);
  add_initial_value_square_1(&C);
  DEBUG_EXPR(NM_display(&C));

  NM_gemm(alpha, M1, M1, beta,  &C);

  NumericsMatrix * Cref= NM_create(NM_DENSE,C.size0, C.size1);
  add_initial_value_square_1(Cref);
  /* gemm by hand */
  dense_gemm_by_hand(alpha, M1->matrix0, M1->matrix0, M1->size0, M1->size1, M1->size0, beta,  Cref->matrix0);
  double err = dense_comparison(C.matrix0, C.size0, C.size1, Cref->matrix0);
  if(err < tol)
  {
    info = 0;
  }

  if(info == 0)
    printf("Step 0 ( C = alpha*A*B + beta*C, double* storage, square matrix ) ok ...\n");
  else
  {
    printf("Step 0 ( C = alpha*A*B + beta*C, double* storage, square matrix) failed ...\n");
    NM_clear(Cref);
    free(C.matrix0);
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

  if(err < tol)
  {
    info = 0;
  }

  if(info == 0)
    printf("Step 1 ( C = alpha*A*B + beta*C, double* storage, non square) ok ...\n");
  else
  {
    printf("Step 1 ( C = alpha*A*B + beta*C, double* storage, non square) failed ...\n");
    free(C2.matrix0);
    NM_clear(C2ref);
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
  /*     Check if it is correct */

  /* C3 and CRef must have the same values.*/

  info = SBM_dense_equal(C3.matrix1, Cref->matrix0, tol);

  if(info == 0)
    printf("Step 2 ( C = alpha*A*B + beta*C, SBM storage) ok ...\n");
  else
  {
    printf("Step 2 ( C = alpha*A*B + beta*C, SBM storage) failed ...\n");
    NM_clear(&C3);
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

  if(info == 0)
    printf("Step 3 ( C = alpha*A*B + beta*C, SBM storage, non square) ok ...\n");
  else
  {
    printf("Step 3 ( C = alpha*A*B + beta*C, SBM storage, non square) failed ...\n");
    NM_clear(&C4);
  }

  /***********************************************************/
  /* C = alpha*A*B + beta*C, NM_SPARSE storage,   square     */
  /***********************************************************/

  NumericsMatrix * M5 = test_matrix_5();
  DEBUG_EXPR(NM_display(M5););
  assert(NM_equal(M5,M2));

  NumericsMatrix * C5 = NM_create(NM_SPARSE,M5->size0, M5->size1);
  NM_triplet_alloc(C5,0);
  C5->matrix2->origin= NSM_TRIPLET;
  add_initial_value_square_1(C5);

  NM_gemm(alpha, M5, M5, beta, C5);

  info = NM_dense_equal(C5,Cref->matrix0,tol);

  if(info == 0)
    printf("Step 4 ( C = alpha*A*B + beta*C, NM_SPARSE storage, square) ok ...\n");
  else
  {
    printf("Step 4 ( C = alpha*A*B + beta*C, NM_SPARSE storage, square) failed ...\n");
    NM_clear(M5);
    NM_clear(C5);
  }

  /***********************************************************/
  /* C = alpha*A*B + beta*C, NM_SPARSE storage, non square  */
  /***********************************************************/

  NumericsMatrix * M6 = test_matrix_6();
  DEBUG_EXPR(NM_display(M6););
  assert(NM_equal(M6,M4));

  NumericsMatrix * C6 = NM_create(NM_SPARSE,M2->size0, M4->size1);
  NM_triplet_alloc(C6,0);
  C6->matrix2->origin= NSM_TRIPLET;
  add_initial_value_rectangle_1(C6);

  NM_gemm(alpha, M2, M4, beta, C6);

  info = NM_dense_equal(C6,C2ref->matrix0,tol);

  if(info == 0)
    printf("Step 5 ( C = alpha*A*B + beta*C, NM_SPARSE storage, non square) ok ...\n");
  else
  {
    printf("Step 5 ( C = alpha*A*B + beta*C, NM_SPARSE storage, non square) failed ...\n");
    NM_clear(M6);
    NM_clear(C6);
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

  if(err < tol)
  {
    info = 0;
  }

  if(info == 0)
    printf("Step 6 ( C = alpha*A*B + beta*C, double* storage, square matrix, empty column of blocks ) ok ...\n");
  else
  {
    printf("Step 6 ( C = alpha*A*B + beta*C, double* storage, square matrix, empty column of blocks) failed ...\n");
    NM_clear(M9);
    NM_clear(C7);
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

  DEBUG_EXPR(NM_display(C8));
  DEBUG_EXPR(NM_dense_display(C3ref->matrix0,M10->size0,M10->size1,M10->size0));

  info = NM_dense_equal(C8,C3ref->matrix0,tol);

  if(info == 0)
    printf("Step 7 ( C = alpha*A*B + beta*C, NM_SPARSE_BLOCK storage,  empty column of blocks) ok ...\n");
  else
  {
    printf("Step 7 ( C = alpha*A*B + beta*C, NM_SPARSE_BLOCK storage,  empty column of blocks) failed ...\n");
    NM_clear(M10);
    NM_clear(C8);
  }



  /* /\************************************************************************************\/ */
  /* /\* C = alpha*A*B + beta*C, SBM storage, empty column of blocks, extra blocks       *\/ */
  /* /\************************************************************************************\/ */

  NumericsMatrix * C20 = test_matrix_20();
  DEBUG_EXPR(NM_display(C20););

  add_initial_value_square_2(C20);

  NM_gemm(alpha, M10, M10, beta, C20);

  DEBUG_EXPR(NM_display(C20));
  DEBUG_EXPR(NM_dense_display(C3ref->matrix0,M10->size0,M10->size1,M10->size0));

  info = NM_dense_equal(C20,C3ref->matrix0,tol);

  if(info == 0)
    printf("Step 8 ( C = alpha*A*B + beta*C, NM_SPARSE_BLOCK storage,  empty column of blocks, extra blocks) ok ...\n");
  else
  {
    printf("Step 8 ( C = alpha*A*B + beta*C, NM_SPARSE_BLOCK storage,  empty column of blocks, extra blocks) failed ...\n");
    NM_clear(C20);
  }

  NM_clear(C20);
  NM_clear(M10);
  NM_clear(C8);
  NM_clear(M9);
  NM_clear(C7);
  NM_clear(M6);
  NM_clear(C6);
  NM_clear(M5);
  NM_clear(C5);
  NM_clear(&C4);
  NM_clear(&C3);
  free(C2.matrix0);
  NM_clear(C2ref);
  NM_clear(Cref);
  free(C.matrix0);
  return info;
}




static int gemm_square_triplet()
{


  int size0 =3;
  int size1 =3;

  // product of triplet matrices into triplet matrix
  NumericsMatrix * A  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(A,0);
  A->matrix2->origin= NSM_TRIPLET;
  NM_zentry(A, 0, 0, 1);
  NM_zentry(A, 0, 1, 2);
  NM_zentry(A, 0, 2, 3);
  NM_zentry(A, 1, 1, 2);
  NM_zentry(A, 1, 2, 3);
  NM_display(A);


  NumericsMatrix * B  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(B,0);
  B->matrix2->origin= NSM_TRIPLET;
  NM_zentry(B, 0, 0, 1);
  NM_zentry(B, 1, 1, 2);
  NM_zentry(B, 2, 2, 3);
  NM_display(B);

  NumericsMatrix * C  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(C,0);
  C->matrix2->origin= NSM_TRIPLET;

  double alpha = 1.0;
  double beta =0.0;
  NM_gemm(alpha, A,B, beta, C);
  NM_display(C);

  NumericsMatrix * Cref  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(Cref,0);
  Cref->matrix2->origin= NSM_TRIPLET;
  NM_zentry(Cref, 0, 0, 1);
  NM_zentry(Cref, 0, 1, 4);
  NM_zentry(Cref, 0, 2, 9);
  NM_zentry(Cref, 1, 1, 4);
  NM_zentry(Cref, 1, 2, 9);
  NM_display(Cref);

  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(C,%i,%i) =%e \n", i,j,NM_get_value(C,i,j) ); */
  /*   } */
  /* } */
  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(Cref,%i,%i) =%e \n", i,j,NM_get_value(Cref,i,j) ); */
  /*   } */
  /* } */


  printf("NM_equal(C,Cref) =%i \n", NM_equal(C,Cref));
  return (int)!NM_equal(C,Cref);


}

static int gemm_square_csc()
{


  int size0 =3;
  int size1 =3;

  // product of csc matrices into csc matrix
  NumericsMatrix * A  = NM_create(NM_SPARSE, size0, size1);
  NM_csc_empty_alloc(A,0);
  A->matrix2->origin= NSM_CSC;
  NM_zentry(A, 0, 0, 1);
  NM_zentry(A, 0, 1, 2);
  NM_zentry(A, 0, 2, 3);
  NM_zentry(A, 1, 1, 2);
  NM_zentry(A, 1, 2, 3);
  /* NM_display(A); */


  NumericsMatrix * B  = NM_create(NM_SPARSE, size0, size1);
  NM_csc_empty_alloc(B,0);
  B->matrix2->origin= NSM_CSC;
  NM_zentry(B, 0, 0, 1);
  NM_zentry(B, 1, 1, 2);
  NM_zentry(B, 2, 2, 3);
  /* NM_display(B); */

  NumericsMatrix * C  = NM_create(NM_SPARSE, size0, size1);
  int nzmax = 10;
  NM_csc_empty_alloc(C,nzmax);
  C->matrix2->origin= NSM_CSC;

  NM_display(C);
  double alpha = 1.0;
  double beta = 0.0;
  NM_gemm(alpha, A,B, beta, C);
  NM_display(C);

  NumericsMatrix * Cref  = NM_create(NM_SPARSE, size0, size1);
  NM_csc_empty_alloc(Cref,0);
  Cref->matrix2->origin= NSM_CSC;
  NM_zentry(Cref, 0, 0, 1);
  NM_zentry(Cref, 0, 1, 4);
  NM_zentry(Cref, 0, 2, 9);
  NM_zentry(Cref, 1, 1, 4);
  NM_zentry(Cref, 1, 2, 9);
  NM_display(Cref);

  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(C,%i,%i) =%e \n", i,j,NM_get_value(C,i,j) ); */
  /*   } */
  /* } */
  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(Cref,%i,%i) =%e \n", i,j,NM_get_value(Cref,i,j) ); */
  /*   } */
  /* } */
  return (int)!NM_equal(C,Cref);;

}
static int gemm_square_triplet_into_csc()
{
  int size0 =3;
  int size1 =3;

  // product of triplet matrices into triplet matrix
  NumericsMatrix * A  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(A,0);
  A->matrix2->origin= NSM_TRIPLET;
  NM_zentry(A, 0, 0, 1);
  NM_zentry(A, 0, 1, 2);
  NM_zentry(A, 0, 2, 3);
  NM_zentry(A, 1, 1, 2);
  NM_zentry(A, 1, 2, 3);
  NM_display(A);




  NumericsMatrix * B  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(B,0);
  B->matrix2->origin= NSM_TRIPLET;
  NM_zentry(B, 0, 0, 1);
  NM_zentry(B, 1, 1, 2);
  NM_zentry(B, 2, 2, 3);
  NM_display(B);

  NumericsMatrix * C  = NM_create(NM_SPARSE, size0, size1);
  int nzmax = 10;
  NM_csc_empty_alloc(C,nzmax);
  C->matrix2->origin= NSM_CSC;

  double alpha = 1.0;
  double beta =0.0;
  NM_gemm(alpha, A,B, beta, C);
  NM_display(C);
  NumericsMatrix * Cref  = NM_create(NM_SPARSE, size0, size1);
  NM_csc_empty_alloc(Cref,0);
  Cref->matrix2->origin= NSM_CSC;
  NM_zentry(Cref, 0, 0, 1);
  NM_zentry(Cref, 0, 1, 4);
  NM_zentry(Cref, 0, 2, 9);
  NM_zentry(Cref, 1, 1, 4);
  NM_zentry(Cref, 1, 2, 9);
  NM_display(Cref);

  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(C,%i,%i) =%e \n", i,j,NM_get_value(C,i,j) ); */
  /*   } */
  /* } */
  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(Cref,%i,%i) =%e \n", i,j,NM_get_value(Cref,i,j) ); */
  /*   } */
  /* } */
  return (int)!NM_equal(C,Cref);;
}

static int gemm_rectangle_triplet()
{


  int size0 =3;
  int size1 =9;

  // product of triplet matrices into triplet matrix
  NumericsMatrix * A  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(A,0);
  A->matrix2->origin= NSM_TRIPLET;
  NM_zentry(A, 0, 0, 1);
  NM_zentry(A, 0, 1, 2);
  NM_zentry(A, 0, 2, 3);
  NM_zentry(A, 1, 1, 2);
  NM_zentry(A, 1, 2, 3);
  NM_zentry(A, 2, 6, 2);
  NM_zentry(A, 2, 5, 22);
  NM_display(A);




  NumericsMatrix * B  = NM_create(NM_SPARSE, size1, size0);
  NM_triplet_alloc(B,0);
  B->matrix2->origin= NSM_TRIPLET;
  NM_zentry(B, 0, 0, 1);
  NM_zentry(B, 1, 1, 2);
  NM_zentry(B, 2, 2, 3);
  NM_zentry(B, 3, 0, 1);
  NM_zentry(B, 4, 1, 2);
  NM_zentry(B, 5, 2, 3);
  NM_display(B);

  NumericsMatrix * C  = NM_create(NM_SPARSE, size0, size0);
  NM_triplet_alloc(C,0);
  C->matrix2->origin= NSM_TRIPLET;

  double alpha = 1.0;
  double beta =0.0;
  NM_gemm(alpha, A,B, beta, C);
  NM_display(C);

  NumericsMatrix * Cref  = NM_create(NM_SPARSE, size0, size0);
  NM_triplet_alloc(Cref,0);
  Cref->matrix2->origin= NSM_TRIPLET;
  NM_zentry(Cref, 0, 0, 1);
  NM_zentry(Cref, 0, 1, 4);
  NM_zentry(Cref, 0, 2, 9);
  NM_zentry(Cref, 1, 1, 4);
  NM_zentry(Cref, 1, 2, 9);
  NM_zentry(Cref, 2, 2, 66);
  NM_display(Cref);

  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(C,%i,%i) =%e \n", i,j,NM_get_value(C,i,j) ); */
  /*   } */
  /* } */
  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(Cref,%i,%i) =%e \n", i,j,NM_get_value(Cref,i,j) ); */
  /*   } */
  /* } */


  printf("NM_equal(C,Cref) =%i \n", NM_equal(C,Cref));
  return (int)!NM_equal(C,Cref);


}

static int gemm_square_triplet_1()
{


  int size0 =3;
  int size1 =3;

  // product of triplet matrices into triplet matrix
  NumericsMatrix * A  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(A,0);
  A->matrix2->origin= NSM_TRIPLET;
  NM_zentry(A, 0, 0, 1);
  NM_zentry(A, 0, 1, 2);
  NM_zentry(A, 0, 2, 3);
  NM_zentry(A, 1, 1, 2);
  NM_zentry(A, 1, 2, 3);
  NM_display(A);


  NumericsMatrix * B  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(B,0);
  B->matrix2->origin= NSM_TRIPLET;
  NM_zentry(B, 0, 0, 1);
  NM_zentry(B, 1, 1, 2);
  NM_zentry(B, 2, 2, 3);
  NM_display(B);

  NumericsMatrix * C  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(C,0);
  C->matrix2->origin= NSM_TRIPLET;

  double alpha = 1.0;
  double beta =0.0;
  NM_gemm(alpha, A,B, beta, C);
  NM_display(C);

  NumericsMatrix * Cref  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(Cref,0);
  Cref->matrix2->origin= NSM_TRIPLET;
  NM_zentry(Cref, 0, 0, 1);
  NM_zentry(Cref, 0, 1, 4);
  NM_zentry(Cref, 0, 2, 9);
  NM_zentry(Cref, 1, 1, 4);
  NM_zentry(Cref, 1, 2, 9);
  NM_display(Cref);

  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(C,%i,%i) =%e \n", i,j,NM_get_value(C,i,j) ); */
  /*   } */
  /* } */
  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(Cref,%i,%i) =%e \n", i,j,NM_get_value(Cref,i,j) ); */
  /*   } */
  /* } */


  printf("NM_equal(C,Cref) =%i \n", NM_equal(C,Cref));
  return (int)!NM_equal(C,Cref);


}

static int gemm_square_csc_1()
{


  int size0 =3;
  int size1 =3;

  // product of csc matrices into csc matrix
  NumericsMatrix * A  = NM_create(NM_SPARSE, size0, size1);
  NM_csc_empty_alloc(A,0);
  A->matrix2->origin= NSM_CSC;
  NM_zentry(A, 0, 0, 1);
  NM_zentry(A, 0, 1, 2);
  NM_zentry(A, 0, 2, 3);
  NM_zentry(A, 1, 1, 2);
  NM_zentry(A, 1, 2, 3);
  /* NM_display(A); */


  NumericsMatrix * B  = NM_create(NM_SPARSE, size0, size1);
  NM_csc_empty_alloc(B,0);
  B->matrix2->origin= NSM_CSC;
  NM_zentry(B, 0, 0, 1);
  NM_zentry(B, 1, 1, 2);
  NM_zentry(B, 2, 2, 3);
  /* NM_display(B); */

  NumericsMatrix * C  = NM_create(NM_SPARSE, size0, size1);
  int nzmax = 10;
  NM_csc_empty_alloc(C,nzmax);
  C->matrix2->origin= NSM_CSC;

  NM_display(C);
  double alpha = 1.0;
  double beta = 0.0;
  NM_gemm(alpha, A,B, beta, C);
  NM_display(C);

  NumericsMatrix * Cref  = NM_create(NM_SPARSE, size0, size1);
  NM_csc_empty_alloc(Cref,0);
  Cref->matrix2->origin= NSM_CSC;
  NM_zentry(Cref, 0, 0, 1);
  NM_zentry(Cref, 0, 1, 4);
  NM_zentry(Cref, 0, 2, 9);
  NM_zentry(Cref, 1, 1, 4);
  NM_zentry(Cref, 1, 2, 9);
  NM_display(Cref);

  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(C,%i,%i) =%e \n", i,j,NM_get_value(C,i,j) ); */
  /*   } */
  /* } */
  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(Cref,%i,%i) =%e \n", i,j,NM_get_value(Cref,i,j) ); */
  /*   } */
  /* } */
  return (int)!NM_equal(C,Cref);;

}
static int gemm_square_triplet_into_csc_1()
{
  int size0 =3;
  int size1 =3;

  // product of triplet matrices into triplet matrix
  NumericsMatrix * A  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(A,0);
  A->matrix2->origin= NSM_TRIPLET;
  NM_zentry(A, 0, 0, 1);
  NM_zentry(A, 0, 1, 2);
  NM_zentry(A, 0, 2, 3);
  NM_zentry(A, 1, 1, 2);
  NM_zentry(A, 1, 2, 3);
  NM_display(A);




  NumericsMatrix * B  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(B,0);
  B->matrix2->origin= NSM_TRIPLET;
  NM_zentry(B, 0, 0, 1);
  NM_zentry(B, 1, 1, 2);
  NM_zentry(B, 2, 2, 3);
  NM_display(B);

  NumericsMatrix * C  = NM_create(NM_SPARSE, size0, size1);
  int nzmax = 10;
  NM_csc_empty_alloc(C,nzmax);
  C->matrix2->origin= NSM_CSC;

  double alpha = 1.0;
  double beta =0.0;
  NM_gemm(alpha, A,B, beta, C);
  NM_display(C);
  NumericsMatrix * Cref  = NM_create(NM_SPARSE, size0, size1);
  NM_csc_empty_alloc(Cref,0);
  Cref->matrix2->origin= NSM_CSC;
  NM_zentry(Cref, 0, 0, 1);
  NM_zentry(Cref, 0, 1, 4);
  NM_zentry(Cref, 0, 2, 9);
  NM_zentry(Cref, 1, 1, 4);
  NM_zentry(Cref, 1, 2, 9);
  NM_display(Cref);

  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(C,%i,%i) =%e \n", i,j,NM_get_value(C,i,j) ); */
  /*   } */
  /* } */
  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(Cref,%i,%i) =%e \n", i,j,NM_get_value(Cref,i,j) ); */
  /*   } */
  /* } */
  return (int)!NM_equal(C,Cref);;
}

static int gemm_rectangle_triplet_1()
{


  int size0 =3;
  int size1 =9;

  // product of triplet matrices into triplet matrix
  NumericsMatrix * A  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(A,0);
  A->matrix2->origin= NSM_TRIPLET;
  NM_zentry(A, 0, 0, 1);
  NM_zentry(A, 0, 1, 2);
  NM_zentry(A, 0, 2, 3);
  NM_zentry(A, 1, 1, 2);
  NM_zentry(A, 1, 2, 3);
  NM_zentry(A, 2, 6, 2);
  NM_zentry(A, 2, 5, 22);
  NM_display(A);




  NumericsMatrix * B  = NM_create(NM_SPARSE, size1, size0);
  NM_triplet_alloc(B,0);
  B->matrix2->origin= NSM_TRIPLET;
  NM_zentry(B, 0, 0, 1);
  NM_zentry(B, 1, 1, 2);
  NM_zentry(B, 2, 2, 3);
  NM_zentry(B, 3, 0, 1);
  NM_zentry(B, 4, 1, 2);
  NM_zentry(B, 5, 2, 3);
  NM_display(B);

  NumericsMatrix * C  = NM_create(NM_SPARSE, size0, size0);
  NM_triplet_alloc(C,0);
  C->matrix2->origin= NSM_TRIPLET;

  double alpha = 1.0;
  double beta =0.0;
  NM_gemm(alpha, A,B, beta, C);
  NM_display(C);

  NumericsMatrix * Cref  = NM_create(NM_SPARSE, size0, size0);
  NM_triplet_alloc(Cref,0);
  Cref->matrix2->origin= NSM_TRIPLET;
  NM_zentry(Cref, 0, 0, 1);
  NM_zentry(Cref, 0, 1, 4);
  NM_zentry(Cref, 0, 2, 9);
  NM_zentry(Cref, 1, 1, 4);
  NM_zentry(Cref, 1, 2, 9);
  NM_zentry(Cref, 2, 2, 66);
  NM_display(Cref);

  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(C,%i,%i) =%e \n", i,j,NM_get_value(C,i,j) ); */
  /*   } */
  /* } */
  /* for(int i=0; i < C->size0; i++) */
  /* { */
  /*   for(int j=0; j < C->size1; j++) */
  /*   { */
  /*     printf("NM_get_value(Cref,%i,%i) =%e \n", i,j,NM_get_value(Cref,i,j) ); */
  /*   } */
  /* } */


  printf("NM_equal(C,Cref) =%i \n", NM_equal(C,Cref));
  return (int)!NM_equal(C,Cref);


}





static int gemm_test_sparse(void)
{

  int info = gemm_square_triplet();
  info += gemm_square_csc();
  info +=  gemm_square_triplet_into_csc();
  info +=  gemm_rectangle_triplet();
  info += gemm_square_triplet_1();
  info += gemm_square_csc_1();
  info +=  gemm_square_triplet_into_csc_1();
  info +=  gemm_rectangle_triplet_1();

  return info;
}



static int NM_gemm_test_all(void)
{

  printf("========= Starts Numerics tests for NumericsMatrix ========= \n");

  int i, nmm = 4 ;
  NumericsMatrix ** NMM = (NumericsMatrix **)malloc(nmm * sizeof(NumericsMatrix *)) ;

  int info = test_build_first_4_NM(NMM);
  if(info != 0)
  {
    printf("Construction failed ...\n");
    return info;
  }
  printf("Construction ok ...\n");

  info = NM_gemm_test(NMM,1.0,0.0);
  if(info != 0)
  {
    printf("End of ProdNumericsMatrix : unsucessfull\n");
    return info;
  }
  info = NM_gemm_test(NMM,1.0,1.0);
  if(info != 0)
  {
    printf("End of ProdNumericsMatrix : unsucessfull\n");
    return info;
  }
  info = NM_gemm_test(NMM,0.0,1.0);
  if(info != 0)
  {
    printf("End of ProdNumericsMatrix : unsucessfull\n");
    return info;
  }
  info = NM_gemm_test(NMM,0.5,0.5);
  if(info != 0)
  {
    printf("End of ProdNumericsMatrix : unsucessfull\n");
    return info;
  }

  info = gemm_test_sparse();
  if(info != 0)
  {
    printf("End of ProdNumericsMatrix : unsucessfull\n");
    return info;
  }


  printf("End of ProdNumericsMatrix ...\n");



  /* free memory */
  for(i = 0 ; i < nmm; i++)
  {
    NM_clear(NMM[i]);
    free(NMM[i]);
  }

  free(NMM);



  printf("========= End Numerics tests for NumericsMatrix ========= \n");
  return info;
}


static int NM_insert_dense_test()
{
  int info = 0;
  size_t Asize0 = 5, Asize1 = 6;
  size_t Bsize0 = 3, Bsize1 = 3;
  size_t start_i = 1, start_j = 1;
  size_t end_i = start_i + Bsize0;
  size_t end_j = start_j + Bsize1;

  /* create and fill the dense matrix A */
  NumericsMatrix * A_dense = NM_create(NM_DENSE, Asize0, Asize1);
  for(size_t i = 0; i < Asize0; ++i)
    for(size_t j = 0; j < Asize1; ++j)
      NM_zentry(A_dense, i, j, 10.0);

  /* create and fill the dense matrix B */
  NumericsMatrix * B_dense = NM_create(NM_DENSE, Bsize0, Bsize1);
  for(size_t i = 0; i < Bsize0; ++i)
    for(size_t j = 0; j < Bsize1; ++j)
      NM_zentry(B_dense, i, j, 999.0);

  /* create an expected result */
  NumericsMatrix * AB_dense = NM_create(NM_DENSE, Asize0, Asize1);
  NM_copy(A_dense, AB_dense);
  for(size_t i = start_i; i < end_i; ++i)
    for(size_t j = start_j; j < end_j; ++j)
      NM_zentry(AB_dense, i, j, 999.0);

  NM_insert(A_dense, B_dense, 1, 1);

  //NM_display(A_dense);
  if(!NM_equal(A_dense, AB_dense))
    info = 1;

  NM_clear(A_dense);
  NM_clear(B_dense);
  NM_clear(AB_dense);
  free(A_dense);
  free(B_dense);
  free(AB_dense);
  printf("== End of test NM_insert_dense_test(result = %d)\n", info);
  return info;
}

static int NM_insert_sparse_test()
{
  int info = 0;
  size_t Asize0 = 5, Asize1 = 6;
  size_t Bsize0 = 3, Bsize1 = 3;
  size_t start_i = 1, start_j = 1;
  size_t end_i = start_i + Bsize0;
  size_t end_j = start_j + Bsize1;

  /* create and fill the dense matrix A */
  NumericsMatrix * A_sparse = NM_create(NM_SPARSE, Asize0, Asize1);
  NM_triplet_alloc(A_sparse, 12);
  A_sparse->matrix2->origin = NSM_TRIPLET;

  for(size_t i = 0; i < Asize0; i += 2)
    for(size_t j = 0; j < Asize1; j += 2)
      NM_zentry(A_sparse, i, j, 10.0);

  /* create and fill the dense matrix B */
  NumericsMatrix * B_sparse = NM_create(NM_SPARSE, Bsize0, Bsize1);
  NM_triplet_alloc(B_sparse, 4);
  B_sparse->matrix2->origin = NSM_TRIPLET;

  for(size_t i = 0; i < Bsize0; i += 2)
    for(size_t j = 0; j < Bsize1; j += 2)
      NM_zentry(B_sparse, i, j, 999.0);

  /* create an expected result */
  NumericsMatrix * AB_sparse = NM_create(NM_SPARSE, Asize0, Asize1);
  NM_copy(A_sparse, AB_sparse);
  for(size_t i = start_i; i < end_i; i += 2)
    for(size_t j = start_j; j < end_j; j += 2)
      NM_zentry(AB_sparse, i, j, 999.0);

  NM_insert(A_sparse, B_sparse, 1, 1);

  //NM_display(A_sparse);
  if(!NM_equal(A_sparse, AB_sparse))
    info = 1;

  NM_clear(A_sparse);
  NM_clear(B_sparse);
  NM_clear(AB_sparse);
  free(A_sparse);
  free(B_sparse);
  free(AB_sparse);

  printf("== End of test NM_insert_sparse_test(result = %d)\n", info);
  return info;
}


CS_INT cs_print(const cs *A, CS_INT brief);


static int NM_gemv_test(NumericsMatrix** MM)
{
  NumericsMatrix* M1 =  MM[0];
  NumericsMatrix* M2 =  MM[1];
  NumericsMatrix* M3 =  MM[2];
  NumericsMatrix* M4 =  MM[3];

  printf("== Numerics tests: NM_gemv(NumericsMatrix,vector) == \n");
  int i, n = M1->size1, m = 4;

  double * x = (double *)malloc(n * sizeof(double));
  double * x2 = (double *)malloc(m * sizeof(double));
  double alpha = 2.3, beta = 1.9;
  double * yref = (double *)malloc(n * sizeof(double));
  double * yref2 = (double *)malloc(n * sizeof(double));;
  double * y = (double *)malloc(n * sizeof(double));
  double * y2 = (double *)malloc(n * sizeof(double));

  for(i = 0; i < n; i++)
  {
    x[i] = i + 1.0;
    yref[i] = 0.1 * i;
    yref2[i] = 0.1 * i;
    y[i] = yref[i];
    y2[i] = yref2[i];
  }
  x2[0] = 0.1;
  x2[1] = 0.2;
  x2[2] = 0.3;
  x2[3] = 0.4;
  int incx = 1, incy = 1;
  cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, alpha, M1->matrix0, n, x, incx, beta, yref, incy);

  NM_gemv(alpha, M1, x, beta, y);

  double tol = 1e-12;
  int info = 0;

  if(NV_equal(y, yref, n, tol))
    printf("Step 0 ( y = alpha*A*x + beta*y, double* storage) ok ...\n");
  else
  {
    printf("Step 0 ( y = alpha*A*x + beta*y, double* storage) failed ...\n");
    info=1;
    goto free_and_return;
  }

  /* sparse storage test for M1 */
  for(i=0; i<n; i++) y[i]=0.1*i;
  NM_csc(M1);
  M1->storageType = NM_SPARSE;

  NM_gemv(alpha, M1, x, beta, y);

  if(NV_equal(y, yref, n, tol))
    printf("Step 0 ( y = alpha*A*x + beta*y, csc storage) ok ...\n");
  else
  {
    printf("Step 0 ( y = alpha*A*x + beta*y, csc storage) failed ...\n");
    info=1;
    goto free_and_return;
  }
  /* end of sparse storage test for M1 */

  cblas_dgemv(CblasColMajor, CblasNoTrans, n, m, alpha, M3->matrix0, n, x2, incx, beta, yref2, incy);

  NM_gemv(alpha, M3, x2, beta, y2);

  if(NV_equal(y2, yref2, n, tol))
    printf("Step 1 ( y = alpha*A*x + beta*y, double* storage, non square) ok ...\n");
  else
  {
    printf("Step 1 ( y = alpha*A*x + beta*y, double* storage, non square) failed ...\n");
    info=1;
    goto free_and_return;
  }

  /* sparse storage test for M3 */
  for(i=0; i<n; i++) y2[i]=0.1*i;
  NM_csc(M3);
  M3->storageType = NM_SPARSE;

  cs_print(M3->matrix2->csc, 0);

  NM_gemv(alpha, M3, x2, beta, y2);

  if(NV_equal(y2, yref2, n, tol))
    printf("Step 1 ( y = alpha*A*x + beta*y, csc storage, non square) ok ...\n");
  else
  {
    printf("Step 1 ( y = alpha*A*x + beta*y, csc storage, non square) failed ...\n");
    info=1;
    goto free_and_return;
  }
  /* end of sparse storage test for M3 */




  /* Sparse Block... */
  for(i = 0; i < n; i++)
  {
    y[i] = 0.1 * i;
    y2[i] = 0.1 * i;
  }
  NM_gemv(alpha, M2, x, beta, y);

  if(NV_equal(y, yref, n, tol))
    printf("Step 2 ( y = alpha*A*x + beta*y, SBM storage) ok ...\n");
  else
  {
    printf("Step 2 ( y = alpha*A*x + beta*y,  SBM  storage) failed ...\n");
    info=1;
    goto free_and_return;
  }

  /* sparse storage test for M2 */
  for(i=0; i<n; i++) y[i]=0.1*i;
  NM_csc(M2);
  M2->storageType = NM_SPARSE;

  NM_gemv(alpha, M2, x, beta, y);

  if(NV_equal(y, yref, n, tol))
    printf("Step 2 ( y = alpha*A*x + beta*y, csc storage) ok ...\n");
  else
  {
    printf("Step 2 ( y = alpha*A*x + beta*y, csc storage) failed ...\n");
    info=1;
    goto free_and_return;
  }
  /* end of sparse storage test for M2 */

  NM_gemv(alpha, M4, x2, beta, y2);

  if(NV_equal(y2, yref2, n, tol))
    printf("Step 3 ( y = alpha*A*x + beta*y, SBM storage, non square) ok ...\n");
  else
  {
    printf("Step 3 ( y = alpha*A*x + beta*y,  SBM storage, non square) failed ...\n");
    info=1;
    goto free_and_return;
  }

  /* sparse storage test for M4 */
  for(i=0; i<n; i++) y2[i]=0.1*i;
  NM_csc(M4);
  M4->storageType = NM_SPARSE;

  NM_gemv(alpha, M4, x2, beta, y2);

  if(NV_equal(y2, yref2, n, tol))
    printf("Step 3 ( y = alpha*A*x + beta*y, csc storage) ok ...\n");
  else
  {
    printf("Step 3 ( y = alpha*A*x + beta*y, csc storage) failed ...\n");
    info=1;
    goto free_and_return;
  }
  /* end of sparse storage test for M4 */


free_and_return:

  free(x);
  free(x2);
  free(y);
  free(y2);
  free(yref);
  free(yref2);

  printf("== End of test NM_gemv(result = %d\n", info);

  return info;
}



static int NM_gemm_test_all2(void)
{

  printf("========= Starts Numerics tests for NumericsMatrix ========= \n");

  int i, nmm = 4 ;
  NumericsMatrix ** NMM = (NumericsMatrix **)malloc(nmm * sizeof(NumericsMatrix *)) ;

  int info = test_build_first_4_NM(NMM);
  if(info != 0)
  {
    printf("Construction failed ...\n");
    return info;
  }
  printf("Construction ok ...\n");
  info = NM_gemv_test(NMM);
  printf("End of NM_gemv_test ...\n");
  if(info != 0) return info;
  /* free memory */

  for(i = 0 ; i < nmm; i++)
  {
    NM_clear(NMM[i]);
    free(NMM[i]);
  }

  free(NMM);



  printf("========= End Numerics tests for NumericsMatrix ========= \n");
  return info;
}



static int test_NM_row_prod(NumericsMatrix* M1, NumericsMatrix* M2)
{

  printf("== Numerics tests: NM_row_prod(NumericsMatrix,vector) == \n");
  int i, n = M1->size1;
  double * x = (double *)malloc(n * sizeof(double));

  for(i = 0; i < n; i++)
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
  for(i = 0; i < sizeY; i++)
    yref[i] = cblas_ddot(n, &(M1->matrix0[min + i]), incx, x, incy);

  NM_row_prod(n, sizeY, min, M1, x, y, 1);
  double tol = 1e-12;
  int info = 0;
  for(i = 0; i < sizeY; i++)
  {
    if(fabs(y[i] - yref[i]) > tol) info = 1;
    /*        printf("%lf\n", fabs(y[i]-yref[i]));  */
    /*           printf("%lf\n", y[i]); */
    /*           printf("%lf\n", yref[i]); */
  }
  if(info == 0)
    printf("Step 0 ( y = subA*x, double* storage) ok ...\n");
  else
    printf("Step 0 ( y = subA*x, double* storage) failed ...\n");

  /* += */
  NM_row_prod(n, sizeY, min, M1, x, y, 0);
  for(i = 0; i < sizeY; i++)
  {
    if(fabs(y[i] - 2 * yref[i]) > tol) info = 1;
    /*        printf("%lf\n", fabs(y[i]-2*yref[i]));  */
    /*           printf("%lf\n", y[i]); */
    /*           printf("%lf\n", 2*yref[i]); */
  }
  if(info == 0)
    printf("Step 1 ( y += subA*x, double* storage) ok ...\n");
  else
    printf("Step 1 ( y += subA*x, double* storage) failed ...\n");





  free(y);
  sizeY = 2;
  int pos = 1; // pos of the required row of blocks
  y = (double *)malloc(sizeY * sizeof(double));

  for(i = 0; i < sizeY; i++)
  {
    y[i] = 0.0;
    yref[i] = cblas_ddot(n, &(M1->matrix0[4 + i]), incx, x, incy);
  }
  /* Sparse ... */
  NM_row_prod(n, sizeY, pos, M2, x, y, 1);
  for(i = 0; i < sizeY; i++)
  {
    if(fabs(y[i] - yref[i]) > tol) info = 1;
    //  printf("%lf\n", fabs(y[i]-yref[i]));
  }
  for(i = 0; i < sizeY; i++)
    yref[i] = cblas_ddot(n, &(M1->matrix0[6 + i]), incx, x, incy);
  NM_row_prod(n, sizeY, pos + 1, M2, x, y, 1);
  for(i = 0; i < sizeY; i++)
  {
    if(fabs(y[i] - yref[i]) > tol) info = 1;
    //  printf("%lf\n", fabs(y[i]-yref[i]));
  }


  if(info == 0)
    printf("Step 2 ( y = subA*x, sparse storage) ok ...\n");
  else
    printf("Step 2 ( y = subA*x,  sparse storage) failed ...\n");

  /* Sparse, += ... */
  NM_row_prod(n, sizeY, pos + 1, M2, x, y, 0);
  for(i = 0; i < sizeY; i++)
  {
    if(fabs(y[i] - 2 * yref[i]) > tol) info = 1;
    /*       printf("%lf\n", fabs(y[i]-yref[i])); */
  }
  if(info == 0)
    printf("Step 3 ( y += subA*x, sparse storage) ok ...\n");
  else
    printf("Step 3 ( y += subA*x,  sparse storage) failed ...\n");


  free(x);
  free(y);
  printf("== End of test NM_row_prod(NumericsMatrix,vector), result = %d\n", info);

  return info;
}


static int NM_row_prod_test(void)
{

  printf("========= Starts Numerics tests for NumericsMatrix ========= \n");

  int i, nmm = 4 ;
  NumericsMatrix ** NMM = (NumericsMatrix **)malloc(nmm * sizeof(NumericsMatrix *)) ;

  int info = test_build_first_4_NM(NMM);
  if(info != 0)
  {
    printf("Construction failed ...\n");
    return info;
  }
  printf("Construction ok ...\n");
  info = test_NM_row_prod(NMM[0], NMM[1]);
  printf("End of Sub-Prod ...\n");
  if(info != 0) return info;

  /* free memory */

  for(i = 0 ; i < nmm; i++)
  {
    NM_clear(NMM[i]);
    /*    if (NMM[i]->matrix0) */
    /*        free(NMM[i]->matrix0); */
    /*    if (NMM[i]->matrix1) */
    /*        SBM_free(NMM[i]->matrix1); */
    free(NMM[i]);
  }

  free(NMM);



  printf("========= End Numerics tests for NumericsMatrix ========= \n");
  return info;
}



static int test_NM_row_prod_no_diag(NumericsMatrix* M1, NumericsMatrix* M2)
{

  printf("== Numerics tests: NM_row_prod_no_diag(NumericsMatrix,vector) == \n");
  int i, n = M1->size1;
  double * x = (double *)malloc(n * sizeof(double));

  for(i = 0; i < n; i++)
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
  for(i = 0; i < sizeY; i++)
  {
    if(fabs(y[i] - yref[i]) > tol) info = 1;
    //  printf("%lf\n", fabs(y[i]-yref[i]));
  }
  NM_row_prod_no_diag(n, sizeY, pos + 1, SIZE_MAX, M2, x, y, NULL, 1);
  yref[0] = 10;
  yref[1] = 14;
  for(i = 0; i < sizeY; i++)
  {
    if(fabs(y[i] - yref[i]) > tol) info = 1;
  }

  if(info == 0)
    printf("Step 2 ( y = subA*x, sparse storage) ok ...\n");
  else
    printf("Step 2 ( y = subA*x,  sparse storage) failed ...\n");

  /* Sparse, += ... */
  NM_row_prod_no_diag(n, sizeY, pos + 1, SIZE_MAX, M2, x, y, NULL, 0);
  for(i = 0; i < sizeY; i++)
  {
    if(fabs(y[i] - 2 * yref[i]) > tol) info = 1;
    /*       printf("%lf\n", fabs(y[i]-yref[i])); */
  }
  if(info == 0)
    printf("Step 3 ( y += subA*x, sparse storage) ok ...\n");
  else
    printf("Step 3 ( y += subA*x,  sparse storage) failed ...\n");


  free(x);
  free(y);
  printf("== End of test NM_row_prod_no_diag(NumericsMatrix,vector), result = %d\n", info);

  return info;
}
static int NM_row_prod_no_diag_test_all(void)
{

  printf("========= Starts Numerics tests for NumericsMatrix ========= \n");

  int i, nmm = 4 ;
  NumericsMatrix ** NMM = (NumericsMatrix **)malloc(nmm * sizeof(NumericsMatrix *)) ;

  int info = test_build_first_4_NM(NMM);
  if(info != 0)
  {
    printf("Construction failed ...\n");
    return info;
  }
  printf("Construction ok ...\n");
  info = test_NM_row_prod_no_diag(NMM[0], NMM[1]);
  printf("End of Sub-Prod no diag ...\n");
  if(info != 0) return info;

  /* free memory */

  for(i = 0 ; i < nmm; i++)
  {
    NM_clear(NMM[i]);
    free(NMM[i]);
  }

  free(NMM);



  printf("========= End Numerics tests for NumericsMatrix ========= \n");
  return info;
}


static int test_NM_row_prod_no_diag_non_square(NumericsMatrix* M3, NumericsMatrix* M4)
{

  printf("== Numerics tests: NM_row_prod_no_diag_non_square(NumericsMatrix,vector) == \n");
  int i,  m = M3->size1;
  double * x = (double *)malloc(m * sizeof(double));

  for(i = 0; i < m; i++)
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
  for(i = 0; i < sizeY; i++)
  {
    if(fabs(y[i] - yref[i]) > tol) info = 1;
    /*       printf("%lf\n", fabs(y[i]-yref[i])); */
    /*       printf("%lf\n", y[i]); */
    /*       printf("%lf\n", yref[i]); */
  }

  printf("\n");
  yref[0] = 10;
  yref[1] = 14;

  NM_row_prod_no_diag(sizeX, sizeY, pos + 1, SIZE_MAX, M4, x, y, NULL, 1);

  for(i = 0; i < sizeY; i++)
  {
    if(fabs(y[i] - yref[i]) > tol) info = 1;
    /*      printf("%lf\n", fabs(y[i]-yref[i])); */
    /*             printf("%lf\n", y[i]); */
    /*             printf("%lf\n", yref[i]); */
  }

  if(info == 0)
    printf("Step 2 ( y = subA*x, sparse storage) ok ...\n");
  else
    printf("Step 2 ( y = subA*x,  sparse storage) failed ...\n");

  /* Sparse, += ... */
  NM_row_prod_no_diag(sizeX, sizeY, pos + 1, SIZE_MAX, M4, x, y, NULL, 0);
  for(i = 0; i < sizeY; i++)
  {
    if(fabs(y[i] - 2 * yref[i]) > tol) info = 1;
    /*           printf("%lf\n", fabs(y[i]-yref[i])); */
    /*             printf("%lf\n", y[i]); */
    /*             printf("%lf\n", yref[i]); */
  }
  if(info == 0)
    printf("Step 3 ( y += subA*x, sparse storage) ok ...\n");
  else
    printf("Step 3 ( y += subA*x,  sparse storage) failed ...\n");


  free(x);
  free(y);
  printf("== End of test NM_row_prod_no_diag_non_square(NumericsMatrix,vector), result = %d\n", info);

  return info;
}
static int NM_row_prod_no_diag_non_square_test(void)
{

  printf("========= Starts Numerics tests for NumericsMatrix ========= \n");

  int i, nmm = 4 ;
  NumericsMatrix ** NMM = (NumericsMatrix **)malloc(nmm * sizeof(NumericsMatrix *)) ;
  int info = test_build_first_4_NM(NMM);
  if(info != 0)
  {
    printf("Construction failed ...\n");
    return info;
  }
  printf("Construction ok ...\n");

  info = test_NM_row_prod_no_diag_non_square(NMM[2], NMM[3]);
  printf("End of Sub-Prod no diag Non Square...\n");
  if(info != 0) return info;

  /* free memory */

  for(i = 0 ; i < nmm; i++)
  {
    NM_clear(NMM[i]);
    free(NMM[i]);
  }

  free(NMM);



  printf("========= End Numerics tests for NumericsMatrix ========= \n");
  return info;
}

static int test_NM_row_prod_non_square(NumericsMatrix* M3, NumericsMatrix* M4)
{

  printf("== Numerics tests: subRowProd_non_square(NumericsMatrix,vector) == \n");
  int i, n = M3->size0, m = M3->size1;
  double * x = (double *)malloc(m * sizeof(double));

  for(i = 0; i < m; i++)
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
  for(i = 0; i < sizeY; i++)
    yref[i] = cblas_ddot(m, &(M3->matrix0[min + i]), incx, x, incy);

  NM_row_prod(sizeX, sizeY, min, M3, x, y, 1);
  double tol = 1e-12;
  int info = 0;
  for(i = 0; i < sizeY; i++)
  {
    if(fabs(y[i] - yref[i]) > tol) info = 1;
    /*       printf("%lf\n", fabs(y[i]-yref[i]));  */
    /*           printf("%lf\n", y[i]); */
    /*           printf("%lf\n", yref[i]); */
  }
  if(info == 0)
    printf("Step 0 ( y = subA*x, double* storage _non_square) ok ...\n");
  else
    printf("Step 0 ( y = subA*x, double* storage _non_square) failed ...\n");

  /* += */
  NM_row_prod(sizeX, sizeY, min, M3, x, y, 0);
  for(i = 0; i < sizeY; i++)
  {
    if(fabs(y[i] - 2 * yref[i]) > tol) info = 1;
    /*         printf("%lf\n", fabs(y[i]-2*yref[i]));  */
    /*           printf("%lf\n", y[i]); */
    /*           printf("%lf\n", 2*yref[i]); */
  }
  if(info == 0)
    printf("Step 1 ( y += subA*x, double* storage _non_square) ok ...\n");
  else
    printf("Step 1 ( y += subA*x, double* storage _non_square) failed ...\n");


  free(y);


  sizeY = 2;
  int pos = 1; // pos of the required row of blocks
  y = (double *)malloc(sizeY * sizeof(double));

  for(i = 0; i < sizeY; i++)
  {
    y[i] = 0.0;
    yref[i] = cblas_ddot(m, &(M3->matrix0[4 + i]), incx, x, incy);
  }
  /* Sparse ... */
  NM_row_prod(sizeX, sizeY, pos, M4, x, y, 1);
  for(i = 0; i < sizeY; i++)
  {
    if(fabs(y[i] - yref[i]) > tol) info = 1;
    //  printf("%lf\n", fabs(y[i]-yref[i]));
  }
  for(i = 0; i < sizeY; i++)
    yref[i] = cblas_ddot(m, &(M3->matrix0[6 + i]), incx, x, incy);
  NM_row_prod(sizeX, sizeY, pos + 1, M4, x, y, 1);
  for(i = 0; i < sizeY; i++)
  {
    if(fabs(y[i] - yref[i]) > tol) info = 1;
    //  printf("%lf\n", fabs(y[i]-yref[i]));
  }


  if(info == 0)
    printf("Step 2 ( y = subA*x, sparse storage _non_square) ok ...\n");
  else
    printf("Step 2 ( y = subA*x,  sparse storage _non_square) failed ...\n");

  /* Sparse, += ... */
  NM_row_prod(sizeX, sizeY, pos + 1, M4, x, y, 0);
  for(i = 0; i < sizeY; i++)
  {
    if(fabs(y[i] - 2 * yref[i]) > tol) info = 1;
    /*       printf("%lf\n", fabs(y[i]-yref[i])); */
  }
  if(info == 0)
    printf("Step 3 ( y += subA*x, sparse storage) ok ...\n");
  else
    printf("Step 3 ( y += subA*x,  sparse storage) failed ...\n");


  free(x);
  free(y);
  printf("== End of test NM_row_prod(NumericsMatrix,vector), result = %d\n", info);

  return info;
}

static int test_NM_row_prod_non_square_test(void)
{

  printf("========= Starts Numerics tests for NumericsMatrix ========= \n");

  int i, nmm = 4 ;
  NumericsMatrix ** NMM = (NumericsMatrix **)malloc(nmm * sizeof(NumericsMatrix *)) ;
  int info = test_build_first_4_NM(NMM);

  if(info != 0)
  {
    printf("Construction failed ...\n");
    return info;
  }
  printf("Construction ok ...\n");
  info = test_NM_row_prod_non_square(NMM[2], NMM[3]);
  printf("End of Sub-Prod Non Square...\n");
  if(info != 0) return info;

  /* free memory */

  for(i = 0 ; i < nmm; i++)
  {
    NM_clear(NMM[i]);
    free(NMM[i]);
  }
  free(NMM);



  printf("========= End Numerics tests for NumericsMatrix ========= \n");
  return info;
}

static int test_NM_iterated_power_method(void)
{

  printf("========= Starts Numerics tests for NumericsMatrix ========= \n");

  int i, nmm = 4 ;
  NumericsMatrix ** NMM = (NumericsMatrix **)malloc(nmm * sizeof(NumericsMatrix *)) ;
  int info = test_build_first_4_NM(NMM);

  if(info != 0)
  {
    printf("Construction failed ...\n");
    return info;
  }
  printf("Construction ok ...\n");

  NumericsMatrix * Id = NM_eye(50);
  double eig = NM_iterated_power_method(Id, 1e-14, 100);
  printf("eigenvalue = %e\n", eig);
  printf("End of iterated power method...\n");

  if(fabs(eig - 1.0) > 1e-10)
    info =1;
  if(info != 0) return info;

  NumericsMatrix * A = NMM[0];
  NumericsMatrix * Atrans =  NM_transpose(A);
  NumericsMatrix * AAT = NM_add(1/2., A, 1/2., Atrans);
  eig = NM_iterated_power_method(AAT, 1e-14, 100);
  printf("largest eigenvalue = %e\n", eig);
  printf("End of iterated power method...\n");



  if(fabs(eig - 9.983560005532535086558710) > 1e-10)
    info =1;
  if(info != 0) return info;


  NumericsMatrix * B = NMM[1];
  NumericsMatrix * Btrans =  NM_transpose(A);
  NumericsMatrix * BBT = NM_add(1/2., B, 1/2., Btrans);
  eig = NM_iterated_power_method(BBT, 1e-14, 100);
  printf("largest eigenvalue = %e\n", eig);
  printf("End of iterated power method...\n");
  if(fabs(eig - 9.983560005532535086558710) > 1e-10)
    info =1;
  if(info != 0) return info;


  /* free memory */

  for(i = 0 ; i < nmm; i++)
  {
    NM_clear(NMM[i]);
    free(NMM[i]);
  }
  free(NMM);
  NM_clear(Id);
  NM_clear(Atrans);
  NM_clear(Btrans);
  NM_clear(AAT);
  NM_clear(BBT);



  printf("========= End Numerics tests for NumericsMatrix ========= \n");
  return info;
}
static int test_NM_scal(void)
{

  printf("========= Starts Numerics tests for NumericsMatrix NM_scal========= \n");

  int i, nmm = 4 ;
  NumericsMatrix ** NMM = (NumericsMatrix **)malloc(nmm * sizeof(NumericsMatrix *)) ;
  int info = test_build_first_4_NM(NMM);

  if(info != 0)
  {
    printf("Construction failed ...\n");
    return info;
  }
  printf("Construction ok ...\n");

  NumericsMatrix * Id = NM_eye(50);
  NM_scal(1e-03,Id);
  printf("NM_get_value(Id,0,0) =%e \n", NM_get_value(Id,0,0));
  printf("End of NM_scal...\n");


  if(fabs(NM_get_value(Id,0,0) - 1e-03) > 1e-10)
    info =1;
  if(info != 0) return info;

  NumericsMatrix * A = NMM[0];
  NM_scal(1e-03,A);
  printf("End of NM_scal...\n");

  if(fabs(NM_get_value(A,3,0) - 5e-03) > 1e-10)
    info =1;
  if(info != 0) return info;

  NumericsMatrix * A_SBM = NMM[1];
  NM_scal(1e-03, A_SBM);
  printf("End of NM_scal...\n");

  if(fabs(NM_get_value(A_SBM,3,0) - 5e-03) > 1e-10)
    info =1;
  if(info != 0) return info;


  NumericsMatrix * B = test_matrix_5();
  /* NM_display(B); */
  NM_scal(1e-03,B);
  /* NM_display(B); */



  /* free memory */

  for(i = 0 ; i < nmm; i++)
  {
    NM_clear(NMM[i]);
    free(NMM[i]);
  }
  free(NMM);
  NM_clear(Id);
  NM_clear(B);


  printf("========= End Numerics tests for NumericsMatrix NM_scal========= \n");
  return info;
}


static int NM_inv_test_sparse(void)
{
  int size0 =10;
  int size1 =10;
  NumericsMatrix * A  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(A,0);
  A->matrix2->origin= NSM_TRIPLET;
  for(int i =0; i < size0; i++)
  {
    for(int j =i; j < size1; j++)
    {
      NM_zentry(A, i, j, i+j+1);
    }
  }

  //NM_zentry(A, size0-1, size0-1, 10);

  NM_display(A);
  FILE * fileout = fopen("dataA.py", "w");
  NM_write_in_file_python(A, fileout);
  fclose(fileout);

  NumericsMatrix * Ainv  = NM_inv(A);


  NumericsMatrix* AAinv = NM_multiply(A,Ainv);
  //NM_display(AAinv);

  NumericsMatrix * Id  = NM_eye(size0);

  //NM_display(Id);


  return  !NM_equal(AAinv, Id);
}



static int test_NM_inv(void)
{

  printf("========= Starts Numerics tests for NumericsMatrix NM_inv ========= \n");

  int i, nmm = 4 ;
  NumericsMatrix ** NMM = (NumericsMatrix **)malloc(nmm * sizeof(NumericsMatrix *)) ;
  int info = test_build_first_4_NM(NMM);

  if(info != 0)
  {
    printf("Construction failed ...\n");
    return info;
  }
  printf("Construction ok ...\n");

  NumericsMatrix * Id = NM_eye(50);
  NM_csc(Id);
  NumericsMatrix * Iinv = NM_inv(Id);
  NM_csc(Iinv);
  NumericsMatrix* IIinv = NM_multiply(Id,Iinv);
  info = !NM_equal(IIinv, Id);
  printf("info : %i\n", info);
  if(info != 0) return info;
  printf("end if test I  ...\n");

  NumericsMatrix * A = NMM[0];
  NumericsMatrix * Ainv = NM_inv(A);
  NumericsMatrix* AAinv = NM_multiply(A,Ainv);
  NumericsMatrix * IA = NM_eye(A->size0);
  info = !NM_compare(AAinv, IA, 1e-14);
  if(info != 0) return info;
  printf("end if test A dense  ...\n");

  NumericsMatrix * B = NMM[1];
  NumericsMatrix * Binv = NM_inv(B);
  NumericsMatrix* BBinv = NM_multiply(B,Binv);
  NumericsMatrix * IB = NM_eye(B->size0);
  info = !NM_compare(BBinv, IB, 1e-14);
  if(info != 0) return info;
  printf("end if test B  SBM ...\n");

  NumericsMatrix * C = test_matrix_5();
  NumericsMatrix * Cinv = NM_inv(C);
  NumericsMatrix* CCinv = NM_multiply(C,Cinv);
  NumericsMatrix * IC = NM_eye(C->size0);
  info = !NM_compare(CCinv, IC, 1e-14);
  if(info != 0) return info;
  printf("end if test C  Sparse ...\n");

  for(i = 0 ; i < nmm; i++)
  {
    NM_clear(NMM[i]);
    free(NMM[i]);
  }
  free(NMM);
  NM_clear(Id);
  NM_clear(Iinv);
  NM_clear(IIinv);

  NM_clear(AAinv);
  NM_clear(Ainv);
  NM_clear(IA);
  NM_clear(BBinv);
  NM_clear(Binv);
  NM_clear(IB);
  NM_clear(CCinv);
  NM_clear(Cinv);
  NM_clear(IC);
  NM_clear(C);

  info = NM_inv_test_sparse();
  if(info ) return info;

  printf("========= End Numerics tests for NumericsMatrix NM_inv ========= \n");
  return info;
}



static int test_NM_gesv_expert_unit(NumericsMatrix * M1, double * b)
{
  int n = M1->size0;
  int info =-1;
  double * y = (double*)malloc(n* sizeof(double));
  for(int j=0; j < n; j++)
    y[j] = b[j];
  NM_gesv_expert(M1, b, NM_PRESERVE);
  NV_display(b,n);
  NM_gemv(-1.0, M1, b, 1.0, y);
  double res = cblas_dnrm2(n,y,1);
  free(y);
  printf("residual = %e\n", res);
  if(fabs(res) >= sqrt(DBL_EPSILON))
    info = 1;
  else
    info=0;
  return info;
}

static int test_NM_gesv_expert(void)
{

  printf("========= Starts Numerics tests for NumericsMatrix ========= \n");

  int i, nmm = 4 ;
  NumericsMatrix ** NMM = (NumericsMatrix **)malloc(nmm * sizeof(NumericsMatrix *)) ;
  int info = test_build_first_4_NM(NMM);

  if(info != 0)
  {
    printf("Construction failed ...\n");
    return info;
  }
  printf("Construction ok ...\n");


  NumericsMatrix * M1 = NULL;
  double * b = NULL;

  int n =0;

  M1 = NMM[0];
  n = M1->size0;
  b = (double*)malloc(n* sizeof(double));
  for(int j=0; j < n; j++)
    b[j] =1.0;
  info = test_NM_gesv_expert_unit(M1, b);
  if(info != 0) return info;

  M1=NMM[1];
  n = M1->size0;
  b = (double*)malloc(n* sizeof(double));
  for(int j=0; j < n; j++)
    b[j] =1.0;
  info = test_NM_gesv_expert_unit(M1, b);
  if(info != 0) return info;

  M1 = test_matrix_5();
  n = M1->size0;
  b = (double*)malloc(n* sizeof(double));
  for(int j=0; j < n; j++)
    b[j] =1.0;
  info = test_NM_gesv_expert_unit(M1, b);
  if(info != 0) return info;

  free(b);

  printf("End of NM_gesv...\n");
  if(info != 0) return info;

  /* free memory */

  for(i = 0 ; i < nmm; i++)
  {
    NM_clear(NMM[i]);
    free(NMM[i]);
  }
  free(NMM);
  printf("========= End Numerics tests for NumericsMatrix ========= \n");
  return info;
}


static int test_NM_posv_expert_unit(NumericsMatrix * M, double * b)
{
  int n = M->size0;
  int info =-1;
  double * y_save = (double*)malloc(n* sizeof(double));
  for(int j=0; j < n; j++)
    y_save[j] = b[j];


  printf("Cholesky solve preserving matrix\n");
  double * y = (double*)malloc(n* sizeof(double));
  for(int j=0; j < n; j++)
    y[j] = b[j];
  NSM_linear_solver_params* p = NSM_linearSolverParams(M);
#ifdef WITH_MUMPS
  p->solver = NSM_MUMPS;
  NM_MUMPS_set_verbosity(M, 1);
#else
  p->solver = NSM_CSPARSE;
#endif

  double res;
#ifndef WITH_MUMPS
  NM_posv_expert(M, b, NM_PRESERVE);
  NV_display(b,n);
  NM_gemv(-1.0, M, b, 1.0, y);
  res = cblas_dnrm2(n,y,1);

  printf("residual = %e\n", res);
  if(fabs(res) >= sqrt(DBL_EPSILON))
  {
    info = 1;
    return info;
  }
  else
    info=0;
#endif

  printf("Cholesky solve keeping factors\n");
  for(int j=0; j < n; j++)
  {
    b[j] = y_save[j];
    y[j]=b[j];
  }
  NumericsMatrix * M_copy = NM_create(NM_SPARSE,M->size0, M->size1);
  NM_copy(M, M_copy);
  NM_posv_expert(M, b, NM_KEEP_FACTORS);
  NV_display(b,n);
  NM_gemv(-1.0, M_copy, b, 1.0, y);
  res = cblas_dnrm2(n,y,1);
  printf("residual = %e\n", res);
  if(fabs(res) >= sqrt(DBL_EPSILON))
  {
    info = 1;
    return info;
  }
  else
    info=0;

  printf("Cholesky solve with given factors\n");

  for(int j=0; j < n; j++)
  {
    y[j]  = 3.0*y_save[j];
    b[j] = y[j];
  }
  NM_posv_expert(M, b, NM_KEEP_FACTORS);
  NV_display(b,n);
  NM_gemv(-1.0, M_copy, b, 1.0, y);
  res = cblas_dnrm2(n,y,1);
  printf("residual = %e\n", res);
  if(fabs(res) >= sqrt(DBL_EPSILON))
  {
    info = 1;
    return info;
  }
  else
    info=0;


  free(y);
  free(y_save);


  return info;
}

static int test_NM_posv_expert(void)
{

  printf("========= Starts Numerics tests for NumericsMatrix ========= \n");

  int i, nmm = 4 ;
  NumericsMatrix ** NMM = (NumericsMatrix **)malloc(nmm * sizeof(NumericsMatrix *)) ;
  int info = test_build_first_4_NM(NMM);

  if(info != 0)
  {
    printf("Construction failed ...\n");
    return info;
  }
  printf("Construction ok ...\n");


  NumericsMatrix * M1 = NULL;
  double * b = NULL;

  int n =0;

  NumericsMatrix *Id = NM_eye(10);
  //NM_scal(Id, 5.0);
  n = Id->size0;
  b = (double*)malloc(n* sizeof(double));
  for(int j=0; j < n; j++)
  {
    b[j] =2.0*j;
    //NM_set_value(Id, j,j, 2.0*j);
  }
  info = test_NM_posv_expert_unit(Id, b);
  if(info != 0) return info;
  NM_clear(Id);
  free(Id);

  NumericsMatrix * Z = NM_create(NM_SPARSE,2,2);
  NM_triplet_alloc(Z,0);
  Z->matrix2->origin= NSM_TRIPLET;
  NM_zentry(Z,0,0,2.0);
  NM_zentry(Z,1,1,2.0);
  NM_zentry(Z,0,1,1.0);
  NM_zentry(Z,1,0,1.0);
  info = test_NM_posv_expert_unit(Z, b);
  if(info != 0) return info;
  NM_clear(Z);
  free(Z);

  M1 = NMM[0];
  NumericsMatrix * M1T = NM_transpose(M1);
  NumericsMatrix * C = NM_create(NM_DENSE,M1->size0,M1->size1);
  NM_gemm(1.0, M1, M1T, 0.0, C);
  n = M1->size0;
  b = (double*)malloc(n* sizeof(double));
  for(int j=0; j < n; j++)
    b[j] =1.0;
  info = test_NM_posv_expert_unit(C, b);
  if(info != 0) return info;
  NM_clear(M1T);
  NM_clear(C);

  /* M1=NMM[1]; */
  /* M1T = NM_transpose(M1); */
  /* C = NM_create(NM_SPARSE_BLOCK,M1->size0,M1->size1); */
  /* NM_gemm(1.0, M1, M1T, 0.0, C); */
  /* n = M1->size0; */
  /* b = (double*)malloc(n* sizeof(double)); */
  /* for (int j=0; j < n; j++) */
  /*   b[j] =1.0; */
  /* info = test_NM_posv_expert_unit(C, b); */
  /* if (info != 0) return info; */

  M1 = test_matrix_5();
  M1T = NM_transpose(M1);
  C = NM_create(NM_SPARSE,M1->size0,M1->size1);
  NM_triplet_alloc(C,0);
  C->matrix2->origin= NSM_TRIPLET;
  NM_gemm(1.0, M1, M1T, 0.0, C);
  n = M1->size0;
  b = (double*)malloc(n* sizeof(double));
  for(int j=0; j < n; j++)
    b[j] =1.0;
  info = test_NM_posv_expert_unit(C, b);
  if(info != 0) return info;

  free(b);



  printf("End of NM_posv...\n");

  /* free memory */

  for(i = 0 ; i < nmm; i++)
  {
    NM_clear(NMM[i]);
    free(NMM[i]);
  }
  free(NMM);
  printf("========= End Numerics tests for NumericsMatrix ========= \n");
  return info;
}
static int test_NM_LU_solve_unit(NumericsMatrix * M1, double * b)
{
  int n = M1->size0;
  int info =-1;
  double * y = (double*)malloc(n* sizeof(double));
  for(int j=0; j < n; j++)
    y[j] = b[j];
  /* We can also write:
     NM_preserve(M1);
     NM_LU_solve(M1, b, 1); */
  NM_LU_solve(NM_preserve(M1), b, 1);
  NV_display(b,n);
  NM_gemv(-1.0, M1, b, 1.0, y);
  double res = cblas_dnrm2(n,y,1);
  free(y);
  printf("residual = %e\n", res);
  if(fabs(res) >= sqrt(DBL_EPSILON))
    info = 1;
  else
    info=0;
  return info;
}

static int test_NM_LU_solve(void)
{

  printf("========= Starts Numerics tests for NumericsMatrix (test_NM_LU_solve) ========= \n");
  /* numerics_set_verbose(2); */
  int i, nmm = 4 ;
  NumericsMatrix ** NMM = (NumericsMatrix **)malloc(nmm * sizeof(NumericsMatrix *)) ;
  int info = test_build_first_4_NM(NMM);

  if(info != 0)
  {
    printf("Construction failed ...\n");
    return info;
  }
  printf("Construction ok ...\n");


  NumericsMatrix * M1 = NULL;
  double * b = NULL;

  int n =0;
  
  M1 = NMM[0];
  n = M1->size0;
  b = (double*)malloc(n* sizeof(double));
  for(int j=0; j < n; j++)
    b[j] =1.0;
  info = test_NM_LU_solve_unit(M1, b);
  if(info != 0) return info;

  M1=NMM[1];
  n = M1->size0;
  b = (double*)malloc(n* sizeof(double));
  for(int j=0; j < n; j++)
    b[j] =1.0;
  info = test_NM_LU_solve_unit(M1, b);
  if(info != 0) return info;

  M1 = test_matrix_5();
  n = M1->size0;
  b = (double*)malloc(n* sizeof(double));
  for(int j=0; j < n; j++)
    b[j] =1.0;
  info = test_NM_LU_solve_unit(M1, b);
  if(info != 0) return info;

  free(b);

  printf("End of NM_gesv...\n");
  if(info != 0) return info;

  /* free memory */

  for(i = 0 ; i < nmm; i++)
  {
    NM_clear(NMM[i]);
    free(NMM[i]);
  }
  free(NMM);
  printf("========= End Numerics tests for NumericsMatrix  (test_NM_LU_solve) ========= \n");
  return info;
}

static int test_NM_Cholesky_solve_unit(NumericsMatrix * M, double * b)
{
  int n = M->size0;
  int info =-1;
  double * y_save = (double*)malloc(n* sizeof(double));
  for(int j=0; j < n; j++)
    y_save[j] = b[j];


  printf("Cholesky solve preserving matrix\n");
  
  double * y = (double*)malloc(n* sizeof(double));
  for(int j=0; j < n; j++)
    y[j] = b[j];
  NSM_linear_solver_params* p = NSM_linearSolverParams(M);
#ifdef WITH_MUMPS
  p->solver = NSM_MUMPS;
  NM_MUMPS_set_verbosity(M, 1);
#else
  p->solver = NSM_CSPARSE;
#endif

  double res;
#ifndef WITH_MUMPS
  NM_preserve(M);
  NM_Cholesky_solve(M, b, 1);
  NV_display(b,n);
  NM_gemv(-1.0, M, b, 1.0, y);
  res = cblas_dnrm2(n,y,1);

  printf("residual = %e\n", res);
  if(fabs(res) >= sqrt(DBL_EPSILON))
  {
    info = 1;
    return info;
  }
  else
    info=0;
#endif

  printf("Cholesky solve without preserving matrix\n");
  NM_unpreserve(M);
  for(int j=0; j < n; j++)
  {
    b[j] = y_save[j];
    y[j]=b[j];
  }
  NumericsMatrix * M_copy = NM_create(NM_SPARSE,M->size0, M->size1);
  NM_copy(M, M_copy);
  NM_Cholesky_solve(M, b, 1);
  NM_gemv(-1.0, M_copy, b, 1.0, y);
  res = cblas_dnrm2(n,y,1);
  printf("residual = %e\n", res);
  if(fabs(res) >= sqrt(DBL_EPSILON))
  {
    info = 1;
    return info;
  }
  else
    info=0;

  printf("Cholesky solve with given factors\n");

  for(int j=0; j < n; j++)
  {
    y[j]  = 3.0*y_save[j];
    b[j] = y[j];
  }
  NM_Cholesky_solve(M, b, 1);
  NV_display(b,n);
  NM_gemv(-1.0, M_copy, b, 1.0, y);
  res = cblas_dnrm2(n,y,1);
  printf("residual = %e\n", res);
  if(fabs(res) >= sqrt(DBL_EPSILON))
  {
    info = 1;
    return info;
  }
  else
    info=0;
  
  free(y);
  free(y_save);


  return info;
}

static int test_NM_Cholesky_solve(void)
{

  printf("========= Starts Numerics tests for NumericsMatrix (test_NM_Cholesky_solve)  ========= \n");

  int i, nmm = 4 ;
  NumericsMatrix ** NMM = (NumericsMatrix **)malloc(nmm * sizeof(NumericsMatrix *)) ;
  int info = test_build_first_4_NM(NMM);

  if(info != 0)
  {
    printf("Construction failed ...\n");
    return info;
  }
  printf("Construction ok ...\n");


  NumericsMatrix * M1 = NULL;
  double * b = NULL;

  int n =0;
  printf("test 1 ...\n");
  NumericsMatrix *Id = NM_eye(10);
  //NM_scal(Id, 5.0);
  n = Id->size0;
  b = (double*)malloc(n* sizeof(double));
  for(int j=0; j < n; j++)
  {
    b[j] =2.0*j;
    //NM_set_value(Id, j,j, 2.0*j);
  }
  info = test_NM_Cholesky_solve_unit(Id, b);
  if(info != 0) return info;
  NM_clear(Id);
  free(Id);
  printf("test 1 ...ok \n");

  printf("test 2 ...\n");
  NumericsMatrix * Z = NM_create(NM_SPARSE,2,2);
  NM_triplet_alloc(Z,0);
  Z->matrix2->origin= NSM_TRIPLET;
  NM_zentry(Z,0,0,2.0);
  NM_zentry(Z,1,1,2.0);
  NM_zentry(Z,0,1,1.0);
  NM_zentry(Z,1,0,1.0);
  info = test_NM_Cholesky_solve_unit(Z, b);
  if(info != 0) return info;
  NM_clear(Z);
  free(Z);
  printf("test 2 ... ok\n");

  printf("test 3 ...\n");
  M1 = NMM[0];
  NumericsMatrix * M1T = NM_transpose(M1);
  NumericsMatrix * C = NM_create(NM_DENSE,M1->size0,M1->size1);
  NM_gemm(1.0, M1, M1T, 0.0, C);
  //NM_display(C);
  n = M1->size0;
  for(int j=0; j < n; j++)
    b[j] =1.0;
  info = test_NM_Cholesky_solve_unit(C, b);
  if(info != 0) return info;
  NM_clear(M1T);
  NM_clear(C);
  printf("test 3 ...ok \n");


  /* M1=NMM[1]; */
  /* M1T = NM_transpose(M1); */
  /* C = NM_create(NM_SPARSE_BLOCK,M1->size0,M1->size1); */
  /* NM_gemm(1.0, M1, M1T, 0.0, C); */
  /* n = M1->size0; */
  /* b = (double*)malloc(n* sizeof(double)); */
  /* for (int j=0; j < n; j++) */
  /*   b[j] =1.0; */
  /* info = test_Cholesky_solve_unit_unit(C, b); */
  /* if (info != 0) return info; */

  printf("test 5 ...\n");
  M1 = test_matrix_5();
  M1T = NM_transpose(M1);
  C = NM_create(NM_SPARSE,M1->size0,M1->size1);
  NM_triplet_alloc(C,0);
  C->matrix2->origin= NSM_TRIPLET;
  NM_gemm(1.0, M1, M1T, 0.0, C);
  n = M1->size0;
  for(int j=0; j < n; j++)
    b[j] =1.0;
  info = test_NM_Cholesky_solve_unit(C, b);
  if(info != 0) return info;
  printf("test 5 ... ok\n");

  free(b);

  printf("End of NM_Cholesky...\n");

  /* free memory */

  for(i = 0 ; i < nmm; i++)
  {
    NM_clear(NMM[i]);
    free(NMM[i]);
  }
  free(NMM);
  printf("========= End Numerics tests for NumericsMatrix (test_NM_Cholesky_solve) ========= \n");
  return info;
}
#ifdef WITH_OPENSSL
#include <openssl/sha.h>
#include <string.h>
int test_NM_compute_values_sha1()
{
  int info = 0;

  NumericsMatrix* M1;
  M1 = NM_create(NM_DENSE, 2, 2);

  NM_zentry(M1, 0, 0, 2.);
  NM_zentry(M1, 0, 1, -1.);
  NM_zentry(M1, 1, 0, 1.);
  NM_zentry(M1, 1, 1, 1.);

  NM_set_values_sha1(M1);


  char sha1_str[SHA_DIGEST_LENGTH*2];
  for(int i = 0; i < SHA_DIGEST_LENGTH; ++i)
  {
    sprintf(&sha1_str[i*2], "%02x",
            (unsigned char) NM_values_sha1(M1)[i]);
  }

  /* zz.c:
     #include <stdio.h>

     int main()
     {
     FILE *file=fopen("t.bin","wb");
     double buffer[4]= {2,-1,1,1};
     fwrite(buffer,sizeof(double),4,file);
     fclose(file);
     }
  */
  /* cc zz.c -o zz && ./zz && cat t.bin | openssl sha1 */
  char ref_str[SHA_DIGEST_LENGTH*2]
    = "264ac112664a06a83ce53a7b35003c152e4cb8ed";

  info = strncmp(sha1_str, ref_str, SHA_DIGEST_LENGTH*2);

  NM_clear(M1);
  free(M1);
  return info;
}

int test_NM_check_values_sha1()
{
  int info = 0;
  NumericsMatrix* M1;
  M1 = NM_create(NM_DENSE, 2, 2);

  NM_zentry(M1, 0, 0, 2.);
  NM_zentry(M1, 0, 1, -1.);
  NM_zentry(M1, 1, 0, 1.);
  NM_zentry(M1, 1, 1, 1.);

  NM_set_values_sha1(M1);

  if (!NM_check_values_sha1(M1))
    info += 1;

  NM_zentry(M1, 0, 0, 3.);

  NM_clearSparse(M1);

  if (NM_check_values_sha1(M1))
    info += 1;

  return info;
}
#endif


static int test_NM_compute_balancing_matrices(void)
{
  int size0 =50;
  int size1 =50;
  NumericsMatrix * A  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(A,0);
  A->matrix2->origin= NSM_TRIPLET;
  for(int i =0; i < size0; i++)
  {
    for(int j =i; j < size1; j++)
    {
      NM_zentry(A, i, j, i+j+1);
    }
  }

  NM_display(A);

  BalancingMatrices * B = NM_BalancingMatrices_new(A);
  NM_compute_balancing_matrices(A, 1e-03, 10, B);

  printf("D1\n:");
  NM_display(B->D1);
  printf("D2\n:");
  NM_display(B->D2);
  printf("A\n:");
  NM_display(B->A);
  int info=0;
  return  info;
}

static int test_NM_compute_balancing_matrices_sym(void)
{
  int size0 =10;
  int size1 =10;
  NumericsMatrix * A  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(A,0);
  A->matrix2->origin= NSM_TRIPLET;
  for(int i =0; i < size0; i++)
  {
    NM_zentry(A, i, i, i+i+1);
    for(int j =i; j < size1; j++)
    {
      if(i != j)
      {
        NM_zentry(A, i, j, i+j+1);
        NM_zentry(A, j, i, i+j+1);
      }
    }
  }

  NM_display(A);


  BalancingMatrices * B = NM_BalancingMatrices_new(A);
  NM_compute_balancing_matrices(A, 1e-03, 10, B);
  printf("D1\n:");
  NM_display(B->D1);
  printf("D2\n:");
  NM_display(B->D2);
  printf("A\n:");
  NM_display(B->A);
  int info=0;
  return  info;
}

static int test_NM_compute_balancing_matrices_rectangle(void)
{
  int size0 =3;
  int size1 =1;
  NumericsMatrix * A  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(A,0);
  A->matrix2->origin= NSM_TRIPLET;

  for(int i =0; i < size0; i++)
  {
    for(int j =0; j < size1; j++)
    {
      NM_zentry(A, i, j, i+j+1);
    }
  }

  NM_display(A);
  BalancingMatrices * B = NM_BalancingMatrices_new(A);
  NM_compute_balancing_matrices(A, 1e-03, 100, B);
  printf("D1\n:");
  NM_display(B->D1);
  printf("D2\n:");
  NM_display(B->D2);
  printf("A\n:");
  NM_display(B->A);
  int info=0;
  return  info;
}


static int test_NM_max_by_columns_and_rows(void)
{
  int size0 =5;
  int size1 =5;
  NumericsMatrix * A  = NM_create(NM_SPARSE, size0, size1);
  NM_triplet_alloc(A,0);
  A->matrix2->origin= NSM_TRIPLET;
  for(int i =0; i < size0; i++)
  {
    for(int j =i; j < size1; j++)
    {
      NM_zentry(A, i, j, i+j+1);
    }
  }

  NM_display(A);

  double * max = (double *) malloc(size0*sizeof(double));

  int info = NM_max_by_columns(A, max);
  printf("Max by columns:\n");
  NV_display(max,size0);
  printf("Max by rows:\n");
  info += NM_max_by_rows(A, max);
  NV_display(max,size0);

  return  info;
}



int main(int argc, char *argv[])
{

#ifdef SICONOS_HAS_MPI
  MPI_Init(&argc, &argv);
#endif

  int info = NM_read_write_test();

  info += NM_add_to_diag3_test_all();

  info += to_dense_test();

  info += NM_gemm_test_all();

  info += NM_gemm_test_all2();

  info += NM_row_prod_test();

  info += NM_row_prod_no_diag_test_all();

  info += NM_row_prod_no_diag_non_square_test();

  info += test_NM_row_prod_non_square_test();

  info += NM_insert_dense_test();

  info += NM_insert_sparse_test();

  info +=    test_NM_row_prod_non_square_test();


  info +=    test_NM_iterated_power_method();

  info +=    test_NM_scal();

  info +=    test_NM_inv();

  info += test_NM_gesv_expert();

  info += test_NM_posv_expert();

  info += test_NM_LU_solve();

  info += test_NM_Cholesky_solve();

#ifdef WITH_OPENSSL
  info += test_NM_compute_values_sha1();
  info += test_NM_check_values_sha1();
#endif

  info += test_NM_compute_balancing_matrices();
  info += test_NM_compute_balancing_matrices_sym();
  info += test_NM_compute_balancing_matrices_rectangle();
  info += test_NM_max_by_columns_and_rows();

#ifdef SICONOS_HAS_MPI
  MPI_Finalize();
#endif

  return info;

}
