// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include <assert.h>  // for assert
#include <stdio.h>   // for NULL
#include <stdlib.h>  // for free

#include "FrictionContactProblem.h"  // for FrictionContactProblem, friction...
#include "NumericsFwd.h"             // for SparseBlockStructuredMatrix, Fri...
#include "NumericsMatrix.h"          // for NM_new_SBM, NumericsMatrix, RawN...
#include "SparseBlockMatrix.h"       // for SparseBlockStructuredMatrix, SBC...

int main(void) {
  int info = 0;

  double q[] = {-1, 1, 3, -1, 1, 3, -1, 1, 3};
  double mu[] = {0.1, 0.1, 0.1};

  size_t row[] = {1, 2, 3};
  size_t column[] = {1, 2, 3};
  int m = 3;
  int n = 3;
  double W[] = {1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1,
                0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1};

  SparseBlockCoordinateMatrix* MC = SBCM_new_3x3(m, n, 3, row, column, W);

  SparseBlockStructuredMatrix* M = SBCM_to_SBM(MC);

  NumericsMatrix* NM = NM_new_SBM(m * 3, n * 3, M);

  FrictionContactProblem* FC = frictionContactProblem_new_with_data(3, 3, NM, q, mu);

  //  frictionContactProblem_display(FC);

  assert(FC->M->sparse_block->blocksize0[2] == 9);
  assert(FC->M->sparse_block->blocksize0[1] == 6);
  assert(FC->M->sparse_block->block[0][0] == 1.);

  free(NM);
  NM = NULL;
  free(M->block);
  M->block = NULL;
  free(M->index1_data);
  free(M->index2_data);
  free(M);
  MC = SBCM_free(MC, SBM_FREE_KEEP_BLOCKS);

  free(FC);

  return info;
}
