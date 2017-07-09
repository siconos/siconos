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
#include "SparseMatrix.h"
#include "SparseBlockMatrix.h"

#define nnz 10

int main(void)
{
  unsigned int Ai[nnz] = {2, 1, 0, 3, 4, 5, 6, 7, 9, 8};
  unsigned int Aj[nnz] = {2, 1, 0, 3, 4, 5, 6, 7, 9, 8};
  double Ax[nnz] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

  SparseBlockCoordinateMatrix mc;

  unsigned int blocksize0[nnz] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  unsigned int blocksize1[nnz] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

  double* block[nnz] = {&Ax[0], &Ax[1], &Ax[2], &Ax[3], &Ax[4], &Ax[5], &Ax[6], &Ax[7], &Ax[8], &Ax[9] };

  mc.nbblocks = nnz;
  mc.blocknumber0 = 10;
  mc.blocknumber1 = 10;
  mc.blocksize0 = blocksize0;
  mc.blocksize1 = blocksize1;

  mc.row = Ai;
  mc.column = Aj;

  mc.block = block;

  SparseBlockStructuredMatrix* m = SBCM_to_SBM(&mc);

  SBM_print(m);

  int info1 = SBM_get_value(m, 0, 0) == 2.;
  int info2 = SBM_get_value(m, 8, 8) == 9.;

  SBM_free_from_SBCM(m);

  return 1-info1 + 1-info2;

}

