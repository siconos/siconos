// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef NUMERICSMATRIX_TEST_FUNCTION_H
#define NUMERICSMATRIX_TEST_FUNCTION_H

#include "NumericsFwd.h"  // for NumericsMatrix, SparseBlockStructuredMatrix

#if defined(__cplusplus)
extern "C" {
#endif
NumericsMatrix* test_matrix_1(void);
NumericsMatrix* test_matrix_2(void);
NumericsMatrix* test_matrix_3(void);
NumericsMatrix* test_matrix_4(void);
NumericsMatrix* test_matrix_5(void);
NumericsMatrix* test_matrix_6(void);

NumericsMatrix* test_matrix_9(void);
NumericsMatrix* test_matrix_10(void);
NumericsMatrix* test_matrix_20(void);

int SBM_dense_equal(SparseBlockStructuredMatrix* M, double* m, double tol);
int NM_dense_equal(NumericsMatrix* M, double* m, double tol);
int test_build_first_4_NM(NumericsMatrix** MM);

int test_SBM_row_to_dense(SparseBlockStructuredMatrix* M);
int test_SBM_row_permutation(SparseBlockStructuredMatrix* M);
int test_SBM_column_permutation(SparseBlockStructuredMatrix* M);
#if defined(__cplusplus)
}
#endif

#endif
