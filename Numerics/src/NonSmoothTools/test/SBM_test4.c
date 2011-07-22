/* Siconos-Numerics, Copyright INRIA 2005-2011.
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


int main(void)
{
  int res;
  printf("========= Starts SBM tests 4 for SBM ========= \n");
  SparseBlockStructuredMatrix M;
  FILE *file = fopen("data/SBM2.dat", "r");
  newFromFileSBM(&M, file);
  printSBM(&M);
  fclose(file);
  /*alloc enough memory */
  SparseMatrix sparseMat;
  res = SBMtoSparseInitMemory(&M, &sparseMat);
  if (res)
  {
    printf("========= Failed SBM tests 4 for SBM  ========= \n");
    return 1;
  }

  res = SBMtoSparse(&M, &sparseMat);
  if (res)
  {
    printf("========= Failed SBM tests 4 for SBM  ========= \n");
    return 1;
  }
  printSparse(&sparseMat);

  int n = M.blocksize0[M.blocknumber0 - 1];
  int m = M.blocksize1[M.blocknumber1 - 1];
  double * denseMat = malloc(n * m * sizeof(double));
  SBMtoDense(&M, denseMat);
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
  SBMfree(&M, NUMERICS_SBM_FREE_BLOCK);
  printf("\n========= Succed SBM tests 4 for SBM  ========= \n");
  return 0;

}

