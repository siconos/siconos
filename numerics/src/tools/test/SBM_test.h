// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*
  Tests functions for SBM structure
*/

int SBM_add_test_all(void);

int SBM_entry_all(void);

int SBM_to_dense_all(void);

int SBM_to_sparse_all(void);

int SBM_multiply_test_all(void);

int gemm_square_triplet(void);
int gemm_square_csc(void);
int gemm_square_triplet_into_csc(void);
int gemm_rectangle_triplet(void);

int SBM_gemm_without_allocation_all(void);

int test_SBM_row_to_dense_all(void);

int test_SBM_column_permutation_all(void);

int test_SBM_row_permutation_all(void);

int SBM_extract_component_3x3_all(void);
