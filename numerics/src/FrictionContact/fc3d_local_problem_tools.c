/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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
#include "fc3d_local_problem_tools.h"
#ifndef __cplusplus
#include <stdbool.h>                 // for false
#endif
#include <stdlib.h>                  // for malloc, NULL
#include "FrictionContactProblem.h"  // for FrictionContactProblem, friction...
#include "NumericsMatrix.h"          // for NM_create_from_data, NumericsMatrix


void fc3d_local_problem_compute_q(FrictionContactProblem * problem, FrictionContactProblem * localproblem, double *reaction, int contact)
{

  double *qLocal = localproblem->q;
  int n = 3 * problem->numberOfContacts;


  int in = 3 * contact, it = in + 1, is = it + 1;

  /* qLocal computation*/
  qLocal[0] = problem->q[in];
  qLocal[1] =  problem->q[it];
  qLocal[2] =  problem->q[is];

  NM_row_prod_no_diag3(n, contact, 3*contact, problem->M, reaction, qLocal, false);

}

void fc3d_local_problem_fill_M(FrictionContactProblem * problem, FrictionContactProblem * localproblem, int contact)
{
  NM_extract_diag_block3(problem->M, contact, &localproblem->M->matrix0);
}


FrictionContactProblem* fc3d_local_problem_allocate(FrictionContactProblem* problem)
{
  /* Connect local solver and local problem*/
  FrictionContactProblem* localproblem =
    (FrictionContactProblem*)malloc(sizeof(FrictionContactProblem));
  localproblem->numberOfContacts = 1;
  localproblem->dimension = 3;
  localproblem->q = (double*)malloc(3 * sizeof(double));
  localproblem->mu = (double*)malloc(sizeof(double));

  if(problem->M->storageType != NM_SPARSE_BLOCK)
  {
    localproblem->M = NM_create_from_data(NM_DENSE, 3, 3,
                                          malloc(9 * sizeof(double)));
  }
  else /* NM_SPARSE_BLOCK */
  {
    localproblem->M = NM_create_from_data(NM_DENSE, 3, 3, NULL); /* V.A. 14/11/2016 What is the interest of this line */
  }
  return localproblem;
}

void fc3d_local_problem_free(FrictionContactProblem* localproblem,
                             FrictionContactProblem* problem)
{
  if(problem->M->storageType == NM_SPARSE_BLOCK)
  {
    /* we release the pointer to avoid deallocation of the diagonal blocks of the original matrix of the problem*/
    localproblem->M->matrix0 = NULL;
  }
  frictionContactProblem_free(localproblem);
}
