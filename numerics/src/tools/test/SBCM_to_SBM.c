/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

/*
  Tests functions for NumericsMatrix structure

 */

#include "SparseBlockMatrix.h"  // for SparseBlockCoordinateMatrix, SBM_get_...

#define nnz 10

int main(void) {
  size_t Ai[nnz] = {2, 1, 0, 3, 4, 5, 6, 7, 9, 8};
  size_t Aj[nnz] = {2, 1, 0, 3, 4, 5, 6, 7, 9, 8};
  double Ax[nnz] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

  SparseBlockCoordinateMatrix mc;

  size_t blocksize0[nnz] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  size_t blocksize1[nnz] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

  double* block[nnz] = {&Ax[0], &Ax[1], &Ax[2], &Ax[3], &Ax[4],
                        &Ax[5], &Ax[6], &Ax[7], &Ax[8], &Ax[9]};

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

  m = SBM_free(m, SBM_FREE_KEEP_BLOCKS_AND_SIZES);  // from sbcm -> keep blocks and blocksizes
  return 1 - info1 + 1 - info2;
}
