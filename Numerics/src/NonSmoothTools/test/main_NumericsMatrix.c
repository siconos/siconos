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
int main(void)
{

  printf("========= Starts Numerics tests for NumericsMatrix ========= \n");

  int i, nmm = 4 ;
  NumericsMatrix ** NMM = malloc(nmm * sizeof(NumericsMatrix *)) ;


  for (i = 0 ; i < nmm; i++)
  {
    NMM[i] = malloc(sizeof(NumericsMatrix));
  }


  int info = test_BuildNumericsMatrix(NMM);
  if (info != 0)
  {
    printf("Construction failed ...\n");
    return info;
  }
  printf("Construction ok ...\n");
  info = test_prodNumericsMatrix(NMM);
  printf("End of ProdNumericsMatrix ...\n");
  if (info != 0) return info;
  /*   i=1; */
  /*   while (i > 0) */
  /*       { */
  info = test_prodNumericsMatrixNumericsMatrix(NMM);
  printf("End of ProdNumericsMatrixNumericsMatrix ...\n");
  /* i++;} */
  if (info != 0) return info;
  info = test_subRowprod(NMM[0], NMM[1]);
  printf("End of Sub-Prod ...\n");
  if (info != 0) return info;
  info = test_subRowprodNonSquare(NMM[2], NMM[3]);
  printf("End of Sub-Prod Non Square...\n");
  if (info != 0) return info;
  info = test_rowProdNoDiag(NMM[0], NMM[1]);
  printf("End of Sub-Prod no diag ...\n");
  if (info != 0) return info;
  info = test_rowProdNoDiagNonSquare(NMM[2], NMM[3]);
  printf("End of Sub-Prod no diag Non Square...\n");
  if (info != 0) return info;

  /* free memory */

  for (i = 0 ; i < nmm; i++)
  {
    if (NMM[i]->matrix0)
      free(NMM[i]->matrix0);
    if (NMM[i]->matrix1)
      freeSBM(NMM[i]->matrix1);
    free(NMM[i]);
  }

  free(NMM);



  printf("========= End Numerics tests for NumericsMatrix ========= \n");
  return info;
}

