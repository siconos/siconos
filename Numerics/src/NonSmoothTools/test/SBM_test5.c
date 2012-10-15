/* Siconos-Numerics, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/

/*
  Tests functions for NumericsMatrix structure

 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "NumericsMatrix.h"
#include "LA.h"
#include <math.h>
#include "numericsMatrixTestFunction.h"
#include "SparseMatrix.h"

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

  SparseBlockStructuredMatrix* m = SBCMToSBM(&mc);

  printSBM(m);

  assert(getValueSBM(m, 0, 0) == 2.);
  assert(getValueSBM(m, 8, 8) == 9.);

  freeSBMFromSBCM(m);

}

