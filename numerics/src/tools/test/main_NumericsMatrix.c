/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
  info = test_NM_row_prod(NMM[0], NMM[1]);
  printf("End of Sub-Prod ...\n");
  if (info != 0) return info;
  info = test_NM_row_prod_non_square(NMM[2], NMM[3]);
  printf("End of Sub-Prod Non Square...\n");
  if (info != 0) return info;
  info = test_NM_row_prod_no_diag(NMM[0], NMM[1]);
  printf("End of Sub-Prod no diag ...\n");
  if (info != 0) return info;
  info = test_NM_row_prod_no_diag_non_square(NMM[2], NMM[3]);
  printf("End of Sub-Prod no diag Non Square...\n");
  if (info != 0) return info;

  /* free memory */

  for (i = 0 ; i < nmm; i++)
  {
    if (NMM[i]->matrix0)
      free(NMM[i]->matrix0);
    if (NMM[i]->matrix1)
      SBM_clear(NMM[i]->matrix1);
    free(NMM[i]);
  }

  free(NMM);



  printf("========= End Numerics tests for NumericsMatrix ========= \n");
  return info;
}

