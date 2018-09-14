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
#include "CSparseMatrix_internal.h"
#include "NumericsSparseMatrix.h"
#include "NumericsMatrix.h"
#include "NumericsVector.h"
#include <math.h>
#include "numericsMatrixTestFunction.h"
#include "SiconosLapack.h"


CS_INT cs_print (const cs *A, CS_INT brief);

static int NM_gemv_test(NumericsMatrix** MM)
{
  NumericsMatrix* M1 =  MM[0];
  NumericsMatrix* M2 =  MM[1];
  NumericsMatrix* M3 =  MM[2];
  NumericsMatrix* M4 =  MM[3];

  printf("== Numerics tests: NM_gemv(NumericsMatrix,vector) == \n");
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
  x2[0] = 0.1;
  x2[1] = 0.2;
  x2[2] = 0.3;
  x2[3] = 0.4;
  int incx = 1, incy = 1;
  cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, alpha, M1->matrix0, n, x, incx, beta, yref, incy);

  NM_gemv(alpha, M1, x, beta, y);

  double tol = 1e-12;
  int info = 0;
  
  if (NV_equal(y, yref, n, tol))
    printf("Step 0 ( y = alpha*A*x + beta*y, double* storage) ok ...\n");
  else
  {
    printf("Step 0 ( y = alpha*A*x + beta*y, double* storage) failed ...\n");
    info=1;
    return info;
  }

  /* sparse storage test for M1 */
  for (i=0; i<n; i++) y[i]=0.1*i;
  NM_csc(M1);
  M1->storageType = NM_SPARSE;

  NM_gemv(alpha, M1, x, beta, y);

  if (NV_equal(y, yref, n, tol))
    printf("Step 0 ( y = alpha*A*x + beta*y, csc storage) ok ...\n");
  else
  {
    printf("Step 0 ( y = alpha*A*x + beta*y, csc storage) failed ...\n");
    info=1;
    return info;
  }
  /* end of sparse storage test for M1 */

  cblas_dgemv(CblasColMajor, CblasNoTrans, n, m, alpha, M3->matrix0, n, x2, incx, beta, yref2, incy);

  NM_gemv(alpha, M3, x2, beta, y2);

  if (NV_equal(y2, yref2, n, tol))
    printf("Step 1 ( y = alpha*A*x + beta*y, double* storage, non square) ok ...\n");
  else
  {
    printf("Step 1 ( y = alpha*A*x + beta*y, double* storage, non square) failed ...\n");
    info=1;
    return info;
  }

  /* sparse storage test for M3 */
  for (i=0; i<n; i++) y2[i]=0.1*i;
  NM_csc(M3);
  M3->storageType = NM_SPARSE;

  cs_print(M3->matrix2->csc, 0);

  NM_gemv(alpha, M3, x2, beta, y2);

  if (NV_equal(y2, yref2, n, tol))
    printf("Step 1 ( y = alpha*A*x + beta*y, csc storage, non square) ok ...\n");
  else
  {
    printf("Step 1 ( y = alpha*A*x + beta*y, csc storage, non square) failed ...\n");
    info=1;
    return info;
  }
  /* end of sparse storage test for M3 */



  
  /* Sparse Block... */
  for (i = 0; i < n; i++)
  {
    y[i] = 0.1 * i;
    y2[i] = 0.1 * i;
  }
  NM_gemv(alpha, M2, x, beta, y);
  
  if (NV_equal(y, yref, n, tol))
    printf("Step 2 ( y = alpha*A*x + beta*y, SBM storage) ok ...\n");
  else
  {
    printf("Step 2 ( y = alpha*A*x + beta*y,  SBM  storage) failed ...\n");
    info=1;
    return info;
  }

  /* sparse storage test for M2 */
  for (i=0; i<n; i++) y[i]=0.1*i;
  NM_csc(M2);
  M2->storageType = NM_SPARSE;

  NM_gemv(alpha, M2, x, beta, y);

  if (NV_equal(y, yref, n, tol))
    printf("Step 2 ( y = alpha*A*x + beta*y, csc storage) ok ...\n");
  else
  {
    printf("Step 2 ( y = alpha*A*x + beta*y, csc storage) failed ...\n");
    info=1;
    return info;
  }
  /* end of sparse storage test for M2 */
  
  NM_gemv(alpha, M4, x2, beta, y2);

  if (NV_equal(y2, yref2, n, tol))
    printf("Step 3 ( y = alpha*A*x + beta*y, SBM storage, non square) ok ...\n");
  else
  {
    printf("Step 3 ( y = alpha*A*x + beta*y,  SBM storage, non square) failed ...\n");
    info=1;
    return info;
  }

  /* sparse storage test for M4 */
  for (i=0; i<n; i++) y2[i]=0.1*i;
  NM_csc(M4);
  M4->storageType = NM_SPARSE;

  NM_gemv(alpha, M4, x2, beta, y2);

  if (NV_equal(y2, yref2, n, tol))
    printf("Step 3 ( y = alpha*A*x + beta*y, csc storage) ok ...\n");
  else
  {
    printf("Step 3 ( y = alpha*A*x + beta*y, csc storage) failed ...\n");
    info=1;
    return info;
  }
  /* end of sparse storage test for M4 */

  free(x);
  free(x2);
  free(y);
  free(y2);
  free(yref);
  free(yref2);

  printf("== End of test NM_gemv(result = %d\n", info);

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
  info = NM_gemv_test(NMM);
  printf("End of NM_gemv_test ...\n");
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

