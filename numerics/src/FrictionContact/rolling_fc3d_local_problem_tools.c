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
#include "RollingFrictionContactProblem.h"
#include "rolling_fc3d_local_problem_tools.h"
#include "NumericsMatrix.h"


void rolling_fc3d_local_problem_compute_q(RollingFrictionContactProblem * problem, RollingFrictionContactProblem * localproblem, double *reaction, int contact)
{

  double *qLocal = localproblem->q;
  int n = 5 * problem->numberOfContacts;


  int in = 5 * contact, it = in + 1, is = it + 1, iv = is + 1, iw = iv + 1;

  /* qLocal computation*/
  qLocal[0] = problem->q[in];
  qLocal[1] =  problem->q[it];
  qLocal[2] =  problem->q[is];
  qLocal[3] =  problem->q[iv];
  qLocal[4] =  problem->q[iw];

  //NM_row_prod_no_diag3(n, contact, 3*contact, problem->M, reaction, qLocal, false);
  NM_row_prod_no_diag(n, 5, contact, 5*contact, problem->M, reaction, qLocal, NULL, false);

}

void rolling_fc3d_local_problem_fill_M(RollingFrictionContactProblem * problem, RollingFrictionContactProblem * localproblem, int contact)
{
  NM_extract_diag_block5(problem->M, contact, &localproblem->M->matrix0);
}


RollingFrictionContactProblem* rolling_fc3d_local_problem_allocate(RollingFrictionContactProblem* problem)
{
  /* Connect local solver and local problem*/
  RollingFrictionContactProblem* localproblem =
    (RollingFrictionContactProblem*)malloc(sizeof(RollingFrictionContactProblem));
  localproblem->numberOfContacts = 1;
  localproblem->dimension = 5;
  localproblem->q = (double*)malloc(5 * sizeof(double));
  localproblem->mu = (double*)malloc(sizeof(double));
  localproblem->mu_r = (double*)malloc(sizeof(double));
  
  if (problem->M->storageType != NM_SPARSE_BLOCK)
  {
    localproblem->M = NM_create_from_data(NM_DENSE, 5, 5,
                                          malloc(25 * sizeof(double)));
  }
  else /* NM_SPARSE_BLOCK */
  {
    localproblem->M = NM_create_from_data(NM_DENSE, 5, 5, NULL); /* V.A. 14/11/2016 What is the interest of this line */
  }
  return localproblem;
}

void rolling_fc3d_local_problem_free(RollingFrictionContactProblem* localproblem,
                      RollingFrictionContactProblem* problem)
{
  if (problem->M->storageType == NM_SPARSE_BLOCK)
  {
    /* we release the pointer to avoid deallocation of the diagonal blocks of the original matrix of the problem*/
    localproblem->M->matrix0 = NULL;
  }
  rollingFrictionContactProblem_free(localproblem);
}
