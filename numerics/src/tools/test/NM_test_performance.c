/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
#include "siconos_debug.h"                       // for DEBUG_EXPR, DEBUG_PRINTF
#include "numericsMatrixTestFunction.h"  // for test_build_first_4_NM, NM_de...
#include "numerics_verbose.h"            // for numerics_error
#include "sanitizer.h"                   // for MSAN_INIT_VAR

#include <time.h>


static int test_NM_row_prod_no_diag3(NumericsMatrix * A)
{
  int n = A->size0;
  int nc = n/3;
  
  double * reaction = (double *) malloc (n *sizeof(double));

  for (int k =0 ; k < n; k++)
    {
      reaction[k] = 1.0;
    }
  double * q = (double *) malloc (3 *sizeof(double));
  for (int k =0 ; k < 3; k++)
    {
      q[k] = 1.0;
    }
  
  // Execute a single test
  long clk_tck = CLOCKS_PER_SEC; 
  clock_t t1 = clock();
  for (int i=0 ; i <1000; i++){
    
  for (int contact =0; contact <nc; contact++)
    {
      NM_row_prod_no_diag3(n, contact, 3*contact, A, reaction, q, false);
    }
  }
  clock_t t2 = clock();
  (void)printf("time (s) : %lf \n", (double)(t2-t1)/(double)clk_tck);
}

int main(int argc, char *argv[])
{

  int info =0;

  
  NumericsMatrix * A_sbm = NM_new_from_filename("./data/matrix_sbm.dat");
  test_NM_row_prod_no_diag3(A_sbm);
 
 
  NumericsMatrix * A_sparse = NM_new_from_filename("./data/matrix_sparse.dat");
  //NM_display(A_sparse);
  test_NM_row_prod_no_diag3(A_sparse);


  A_sparse->matrix2->csr = NM_csc_trans(A_sparse);
  A_sparse->matrix2->origin = NSM_CSR;
  test_NM_row_prod_no_diag3(A_sparse);
  
  
  
  return info;
}
