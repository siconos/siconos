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



#ifndef NUMERICSMATRIX_TEST_FUNCTION_H
#define NUMERICSMATRIX_TEST_FUNCTION_H

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif
  NumericsMatrix * test_matrix_1(void);
  NumericsMatrix * test_matrix_2(void);
  NumericsMatrix * test_matrix_3(void);
  NumericsMatrix * test_matrix_4(void);
  NumericsMatrix * test_matrix_5(void);
  NumericsMatrix * test_matrix_6(void);
  
  NumericsMatrix * test_matrix_10(void);
  int SBM_dense_equal(SparseBlockStructuredMatrix * M, double * m, double tol);
  int NM_dense_equal(NumericsMatrix * M, double * m, double tol);
  int test_build_first_4_NM(NumericsMatrix** MM);
  int test_NM_row_prod(NumericsMatrix* M1, NumericsMatrix* M2);
  int test_NM_row_prod_no_diag(NumericsMatrix* M1, NumericsMatrix* M2);
  int test_NM_row_prod_non_square(NumericsMatrix* M3, NumericsMatrix* M4);
  int test_NM_row_prod_no_diag_non_square(NumericsMatrix* M3, NumericsMatrix* M4);
  int test_SBM_row_to_dense(SparseBlockStructuredMatrix *M);
  int test_SBM_row_permutation(SparseBlockStructuredMatrix *M);
  int test_SBM_column_permutation(SparseBlockStructuredMatrix *M);
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif

#endif


