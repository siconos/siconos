/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
  int test_BuildNumericsMatrix(NumericsMatrix** MM);
  int test_prodNumericsMatrix(NumericsMatrix** MM);
  int test_prodNumericsMatrixNumericsMatrix(NumericsMatrix** MM);
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


