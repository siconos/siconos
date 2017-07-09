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
#include <stdio.h>
#include <stdlib.h>
#include "NonSmoothDrivers.h"
#include "frictionContact_test_function.h"
#include "SparseBlockMatrix.h"
#include "FrictionContactProblem.h"

int main(void)
{
  int info = 0 ;

  double q[] = { -1, 1, 3, -1, 1, 3, -1, 1, 3};
  double mu[] = {0.1, 0.1, 0.1};

  unsigned int row[] = {1, 2, 3};
  unsigned int column[] = {1, 2, 3};
  int m = 3;
  int n = 3;
  double W[] = {1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1};


  SparseBlockCoordinateMatrix* MC =  SBCM_new_3x3(m, n, 3, row, column, W);

  SparseBlockStructuredMatrix* M = SBCM_to_SBM(MC);

  NumericsMatrix* NM = NM_new_SBM(m * 3, n * 3, M);

  FrictionContactProblem* FC = frictionContactProblem_new(3, 3, NM, q, mu);

//  frictionContact_display(FC);

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
   SBCM_free_3x3(MC);
  free(MC);
  free(FC);

  return info;
}
