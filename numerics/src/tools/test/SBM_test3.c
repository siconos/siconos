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

  printf("========= Starts SBM tests 3 for SBM ========= \n");
  SparseBlockStructuredMatrix M;
  FILE *file = fopen("data/SBM1.dat", "r");
  SBM_new_from_file(&M, file);
  fclose(file);
  /*alloc enough memory */
  int res = test_SBM_column_permutation(&M);
  if (res)
  {
    printf("========= Failed SBM tests 3 for SBM  ========= \n");
    return 1;
  }
  SBMfree(&M, NUMERICS_SBM_FREE_BLOCK);
  file = fopen("data/SBM2.dat", "r");
  SBM_new_from_file(&M, file);
  fclose(file);
  res = test_SBM_column_permutation(&M);
  if (res)
  {
    printf("========= Failed SBM tests 3 for SBM  ========= \n");
    return 1;
  }
  SBMfree(&M, NUMERICS_SBM_FREE_BLOCK);
  printf("\n========= Succed SBM tests 3 for SBM  ========= \n");
  return 0;

}

