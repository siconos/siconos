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
  Tests functions for SBM structure
*/


int SBM_add_test_all(void);

int SBM_zentry_all(void);

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
