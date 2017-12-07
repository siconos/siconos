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
#include "SparseMatrix_internal.h"
#include "NumericsMatrix.h"
#include <math.h>
#include "numericsMatrixTestFunction.h"
#include "SparseMatrix.h"
#include"SparseBlockMatrix.h"

int main(void)
{
  int res;
  printf("========= Starts SBM tests 4 for SBM ========= \n");
  SparseBlockStructuredMatrix M;
  FILE *file = fopen("data/SBM2.dat", "r");
  SBM_new_from_file(&M, file);
  SBM_print(&M);
  fclose(file);
  /*alloc enough memory */
  CSparseMatrix sparseMat;
  res = SBM_to_sparse_init_memory(&M, &sparseMat);
  if (res)
  {
    printf("========= Failed SBM tests 4 for SBM  ========= \n");
    return 1;
  }

  res = SBM_to_sparse(&M, &sparseMat);
  if (res)
  {
    printf("========= Failed SBM tests 4 for SBM  ========= \n");
    return 1;
  }
  cs_print(&sparseMat, 1);
  cs_spfree_on_stack(&sparseMat);

  int n = M.blocksize0[M.blocknumber0 - 1];
  int m = M.blocksize1[M.blocknumber1 - 1];
  double * denseMat = (double *)malloc(n * m * sizeof(double));
  SBM_to_dense(&M, denseMat);
  if (res)
  {
    printf("========= Failed SBM tests 4 for SBM  ========= \n");
    return 1;
  }
  printf("[");
  for (int i = 0; i < n * m; i++)
  {
    printf("%lf ", denseMat[i]);
    if ((i + 1) % m == 0)
      printf("\n");
  }
  printf("]");
  printf("\n (warning: column-major) \n");

  free(denseMat);
  printf("NUMERICS_SBM_FREE_BLOCK value %d", NUMERICS_SBM_FREE_BLOCK);
  SBMfree(&M, NUMERICS_SBM_FREE_BLOCK);
  printf("\n========= Succed SBM tests 4 for SBM  ========= \n");
  return 0;

}

