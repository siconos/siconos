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

/*
  Tests functions for NumericsMatrix structure

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NumericsMatrix.h"
#include <math.h>
#include "numericsMatrixTestFunction.h"
#include "SparseBlockMatrix.h"

int main(void)
{

  printf("========= Starts SBM tests 7 for SBM ========= \n");
  SparseBlockStructuredMatrix M;
  FILE *file = fopen("data/SBM1.dat", "r");
  SBM_new_from_file(&M, file);
  fclose(file);
  SBM_print(&M);
  
  SparseBlockStructuredMatrix * N = SBM_new();
  unsigned int row_components[1] = {0};
  unsigned int row_components_size =1;
  unsigned int col_components[1] = {0};
  unsigned int col_components_size =1;
  SBM_extract_component_3x3(&M, N, row_components, row_components_size, col_components, col_components_size   );
  SBM_print(N);

  SparseBlockStructuredMatrix * T = SBM_new();
  unsigned int row_components_T[2] = {1,2};
  unsigned int row_components_size_T =2;
  unsigned int col_components_T[2] = {1,2};
  unsigned int col_components_size_T =2;
  SBM_extract_component_3x3(&M, T, row_components_T, row_components_size_T, col_components_T, col_components_size_T   );
  SBM_print(T);

  SparseBlockStructuredMatrix * NT = SBM_new();
  unsigned int row_components_NT[2] = {0};
  unsigned int row_components_size_NT =1;
  
  unsigned int col_components_NT[2] = {1,2};
  unsigned int col_components_size_NT =2;
  SBM_extract_component_3x3(&M, NT, row_components_NT, row_components_size_NT, col_components_NT, col_components_size_NT   );
  SBM_print(NT);
  
  SparseBlockStructuredMatrix * TN = SBM_new();
  unsigned int row_components_TN[2] = {1,2};
  unsigned int row_components_size_TN =2;
  
  unsigned int col_components_TN[2] = {0};
  unsigned int col_components_size_TN =1;
  SBM_extract_component_3x3(&M, TN, row_components_TN, row_components_size_TN, col_components_TN, col_components_size_TN   );
  SBM_print(TN);
  
  
  
  int res = test_SBM_row_to_dense(&M);
  if (res)
  {
    printf("========= Failed SBM tests 7 for SBM  ========= \n");
    return 1;
  }

  SBMfree(&M, NUMERICS_SBM_FREE_BLOCK);
  
  return 0;

}

