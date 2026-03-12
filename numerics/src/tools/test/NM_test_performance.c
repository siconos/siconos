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
#include "CSparseMatrix.h"
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


static double test_NM_row_prod_no_diag3(NumericsMatrix * A, int number_of_prod)
{
  int n = A->size0;
  int nc = n/3;


  if (A->storageType == NM_SPARSE_BLOCK)
    {
      printf("test_NM_row_prod_no_diag3 for a storage SBM. n =%i, number_of_prod= %i\n", A->size0, number_of_prod);
    }
  else if (A->storageType == NM_SPARSE)
     {
       if (A->matrix2->origin == NSM_CSR)
	 printf("test_NM_row_prod_no_diag3 for a storage CSR. n =%i, number_of_prod= %i\n", A->size0, number_of_prod);
       if (A->matrix2->origin == NSM_CSC)
	 printf("test_NM_row_prod_no_diag3 for a storage CSC. n =%i, number_of_prod= %i\n", A->size0, number_of_prod);
       if (A->matrix2->origin == NSM_TRIPLET)
	 printf("test_NM_row_prod_no_diag3 for a storage COO. n =%i, number_of_prod= %i\n", A->size0, number_of_prod);

     }

  double * reaction = (double *) malloc (n *sizeof(double));

  for (int k =0 ; k < n; k++)
    {
      reaction[k] = 1.0;
    }
  double * q = (double *) malloc (n *sizeof(double));
  for (int k =0 ; k < n; k++)
    {
      q[k] = 1.0;
    }

  // Execute a single test
  long clk_tck = CLOCKS_PER_SEC;
  clock_t t1 = clock();
  for (int i=0 ; i <number_of_prod; i++)

  /* while(1){ */
    {
  for (int contact =0; contact <nc; contact++)
    {
      NM_row_prod_no_diag3(n, contact, 3*contact, A, reaction, &q[contact], false);
    }
  }
  clock_t t2 = clock();
  printf("norm of result = %e\n",cblas_dnrm2(n,q,1) );

  (void)printf("time (s) : %lf \n", (double)(t2-t1)/(double)clk_tck);

  free(reaction);
  free(q);



  return (double)(t2-t1)/(double)clk_tck;
}

static double test_NM_prod_mv_3x3(NumericsMatrix * A, int number_of_prod)
{
  int n = A->size0;
  int nc = n/3;


  if (A->storageType == NM_SPARSE_BLOCK)
    {
      printf("test_NM_prod_mv_3x3 for a storage SBM. n =%i, number_of_prod= %i\n", A->size0, number_of_prod);
    }
  else if (A->storageType == NM_SPARSE)
     {
       if (A->matrix2->origin == NSM_CSR)
	 printf("test_NM_prod_mv_3x3 for a storage CSR. n =%i, number_of_prod= %i\n", A->size0, number_of_prod);
       if (A->matrix2->origin == NSM_CSC)
	 printf("test_NM_prod_mv_3x3 for a storage CSC. n =%i, number_of_prod= %i\n", A->size0, number_of_prod);
       if (A->matrix2->origin == NSM_TRIPLET)
	 printf("test_NM_prod_mv_3x3 for a storage COO. n =%i, number_of_prod= %i\n", A->size0, number_of_prod);

     }

  double * reaction = (double *) malloc (n *sizeof(double));

  for (int k =0 ; k < n; k++)
    {
      reaction[k] = 1.0;
    }
  double * q = (double *) malloc (n *sizeof(double));
  for (int k =0 ; k < n; k++)
    {
      q[k] = 1.0;
    }

  // Execute a single test
  long clk_tck = CLOCKS_PER_SEC;
  clock_t t1 = clock();




  for (int i=0 ; i <number_of_prod; i++)

    NM_prod_mv_3x3(n, n, A, reaction, q);


  clock_t t2 = clock();
  (void)printf("time (s) : %lf \n", (double)(t2-t1)/(double)clk_tck);
  //printf("norm of result = %e\n",cblas_dnrm2(n,q,1) );
  free(reaction);
  free(q);



  return (double)(t2-t1)/(double)clk_tck;
}
int main(int argc, char *argv[])
{

  int info =0;
  int number_of_prod =1000;


  /* NumericsMatrix * A_sbm = NM_new_from_filename("./data/matrix_sbm.dat"); */
  /* test_NM_row_prod_no_diag3(A_sbm,number_of_prod); */

  NumericsMatrix * A_sparse = NM_new_from_filename("./data/matrix_sparse.dat");
  //NM_display(A_sparse);
  test_NM_row_prod_no_diag3(A_sparse,number_of_prod);



  /* NumericsMatrix * A_sparse_csc = NM_new_from_filename("./data/matrix_sparse.dat"); */
  /* A_sparse_csc->matrix2->csc = NM_csc(A_sparse_csc); */
  /* A_sparse_csc->matrix2->origin = NSM_CSC; */
  /* //NM_display(A_sparse_csc); */
  /* test_NM_row_prod_no_diag3(A_sparse_csc,number_of_prod); */


  /* NumericsMatrix * A_sparse_remove = NM_new_from_filename("./data/matrix_sparse.dat"); */

  /* CSparseMatrix * a_remove = NSM_remove_diagonal_blocks(A_sparse_remove, 3); */

  /* NM_clearSparseStorage(A_sparse_remove); */
  /* A_sparse_remove->matrix2->triplet= a_remove; */
  /* A_sparse_remove->matrix2->origin = NSM_TRIPLET; */

  /* //NM_display(A_sparse_remove); */

  /* test_NM_row_prod_no_diag3(A_sparse_remove,number_of_prod); */



  /* printf("\n A sparse remove nnz =%li \n", NM_nnz(A_sparse_remove) ); */
  /* printf("A sparse nnz =%li \n", NM_nnz(A_sparse) ); */
  /* printf("SBM nnz =%li \n", NM_nnz(A_sbm) ); */

  /* printf("Estimation of number of mul operation\n"); */
  /* printf("sparse = %li\n", NM_nnz(A_sparse)*number_of_prod  ); */
  /* printf("sparse remove = %li\n", NM_nnz(A_sparse_remove)*number_of_prod  ); */

  /* printf("sbm = %li\n", (NM_nnz(A_sbm) - (A_sbm->size0/3)*9) *number_of_prod ); */

  /* test_NM_prod_mv_3x3(A_sbm,number_of_prod); */
  /* test_NM_prod_mv_3x3(A_sparse,number_of_prod); */
  /* test_NM_prod_mv_3x3(A_sparse_csc,number_of_prod); */



  return info;
}
