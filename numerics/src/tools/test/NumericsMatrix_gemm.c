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
/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"




/* ============================================================================================================================== */

static int NM_gemm_test(NumericsMatrix** MM, double alpha, double beta)
{


  NumericsMatrix * M1 = MM[0];
  NumericsMatrix * M2 = MM[1];
  NumericsMatrix * M3 = MM[2];
  NumericsMatrix * M4 = MM[3];


  int info = -1;
  printf("== Numerics tests: NM_gemm(NumericsMatrix,NumericsMatrix) == \n");
  int i, j, k;
  
  double tol = 1e-12;

  double * C2ref = 0;

  SparseBlockStructuredMatrix * SBM3 = 0;
  SparseBlockStructuredMatrix * SBM4 = 0;

  NumericsMatrix C;
  NM_null(&C);

  C.storageType = NM_DENSE;
  C.size0 = M1->size0;
  C.size1 = M1->size1;
  C.matrix0 = (double *)calloc(C.size0 * C.size1 , sizeof(double));
  MSAN_INIT_VAR(C.matrix0, C.size0 * C.size1);

  for (i=0; i < 4 ; i++)
  {
    for (j=0; j < 4 ; j++)
      C.matrix0[i + j * C.size0 ] = 1.0;
  }
  for (i=0; i < 4 ; i++)
  {
    for (j=4; j < 6 ; j++)
      C.matrix0[i + j * C.size0 ] = 2.0;
  }
  for (i=4; i < 6 ; i++)
  {
    for (j=4; j < 6 ; j++)
      C.matrix0[i + j * C.size0 ] = 3.0;
  }
  for (i=4; i < 6 ; i++)
  {
    for (j=6; j < 8 ; j++)
      C.matrix0[i + j * C.size0 ] = 4.0;
  }
  for (i=6; i < 8 ; i++)
  {
    for (j=0; j < 4 ; j++)
      C.matrix0[i + j * C.size0 ] = 5.0;
  }
  for (i=6; i < 8 ; i++)
  {
    for (j=6; j < 8 ; j++)
      C.matrix0[i + j * C.size0 ] = 6.0;
  }
  DEBUG_EXPR(NM_display(&C));

  NM_gemm(alpha, M1, M1, beta,  &C);
  
  double * Cref = (double *)calloc(C.size0 * C.size1, sizeof(double));
  for (i=0; i < 4 ; i++)
  {
    for (j=0; j < 4 ; j++)
      Cref[i + j * C.size0 ] = 1.0;
  }
  for (i=0; i < 4 ; i++)
  {
    for (j=4; j < 6 ; j++)
      Cref[i + j * C.size0 ] = 2.0;
  }
  for (i=4; i < 6 ; i++)
  {
    for (j=4; j < 6 ; j++)
      Cref[i + j * C.size0 ] = 3.0;
  }
  for (i=4; i < 6 ; i++)
  {
    for (j=6; j < 8 ; j++)
      Cref[i + j * C.size0 ] = 4.0;
  }
  for (i=6; i < 8 ; i++)
  {
    for (j=0; j < 4 ; j++)
      Cref[i + j * C.size0 ] = 5.0;
  }
  for (i=6; i < 8 ; i++)
  {
    for (j=6; j < 8 ; j++)
      Cref[i + j * C.size0 ] = 6.0;
  }
  double sum;
  for (i = 0; i < C.size0; i++)
  {
    for (j = 0; j < C.size1; j++)
    {
      sum = beta* Cref[i + j * C.size0];
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
  C3.storageType = NM_SPARSE_BLOCK;
  C3.size0 = M2->size0;
  C3.size1 = M2->size1;
  SBM3 = SBM_new();
  C3.matrix1 = SBM3;
  SBM_alloc_for_gemm(M2->matrix1, M2->matrix1, SBM3);

  DEBUG_EXPR(SBM_print(SBM3););
  
  
  NM_gemm(alpha, M2, M2, beta,  &C3);
  DEBUG_EXPR(NM_display(&C3));
  DEBUG_EXPR(NM_dense_display(Cref,M2->size0,M2->size1,M2->size0));
  /*     Check if it is correct */

  /* C3 and CRef must have the same values.*/

  info = SBM_dense_equal(C3.matrix1, Cref, tol);
 
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
  SBM4 = SBM_new();
  C4.matrix1 = SBM4;
  SBM_alloc_for_gemm(M2->matrix1, M4->matrix1, SBM4);
  
  NM_gemm(alpha, M2, M4, beta,  &C4);
  
  DEBUG_EXPR(NM_display(&C4));
  DEBUG_EXPR(NM_dense_display(C2ref,M2->size0,M4->size1,M2->size0));
  /*     Check if it is correct */
  /* C4 and C2Ref must have the same values.*/

  info = SBM_dense_equal(C4.matrix1, C2ref, tol);
  
  if (info == 0)
    printf("Step 3 ( C = alpha*A*B + beta*C, SBM storage) ok ...\n");
  else
  {
    printf("Step 3 ( C = alpha*A*B + beta*C, SBM storage) failed ...\n");
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
  /* if (info != 0) */
  /* { */
  /*   printf("End of ProdNumericsMatrix : unsucessfull\n"); */
  /*   return info; */
  /* } */
  /* info = NM_gemm_test(NMM,1.0,1.0); */
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

