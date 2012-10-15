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
#include <stdio.h>
#include <stdlib.h>
#include "NonSmoothDrivers.h"
#include "frictionContact_test_function.h"


int main(void)
{
  int info = 0 ;

  double q[9] = { -1, 1, 3, -1, 1, 3, -1, 1, 3};
  double mu[3] = {0.1, 0.1, 0.1};

  unsigned int row[3] = {1, 2, 3};
  unsigned int column[3] = {1, 2, 3};
  int m = 3;
  int n = 3;
  double W[27] = {1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1};


  SparseBlockCoordinateMatrix* MC = newSparseBlockCoordinateMatrix3x3fortran(m, n, 3, row, column, W);

  SparseBlockStructuredMatrix* M = SBCMToSBM(MC);

  NumericsMatrix* NM = newSparseNumericsMatrix(m * 3, n * 3, M);

  FrictionContactProblem* FC = frictionContactProblem_new(3, 3, NM, q, mu);

  frictionContact_display(FC);

  assert(FC->M->matrix1->blocksize0[2] == 9);
  assert(FC->M->matrix1->blocksize0[1] == 6);
  assert(FC->M->matrix1->block[0][0] == 1.);

  free(NM);
  NM = NULL;
  free(M->block);
  M->block = NULL;
  free(M->index1_data);
  free(M->index2_data);
  free(M);
  freeSparseBlockCoordinateMatrix3x3fortran(MC);

  return info;
}
