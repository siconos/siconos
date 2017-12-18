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
#include <string.h>
#include <math.h>
#include <float.h>

#include "Relay_Solvers.h"
#include "LCP_Solvers.h"
#include <assert.h>
#include "NumericsMatrix.h"

void relay_to_lcp(RelayProblem* problem, LinearComplementarityProblem * lcp_problem)
{
  lcp_problem->size = 2 * problem->size ;
  lcp_problem->M = NM_new();
  NM_fill(lcp_problem->M, NM_DENSE, lcp_problem->size, lcp_problem->size, malloc(lcp_problem->size * lcp_problem->size * sizeof(double)));
  lcp_problem->q = (double*)malloc(lcp_problem->size * sizeof(double));

  int i, j;
  for (i = 0; i < problem->size; i++)
  {
    for (j = 0; j < problem->size; j++)
    {
      lcp_problem->M->matrix0[i + j * lcp_problem->size] =  problem->M->matrix0[i + j * problem->size];
    }
  }
  for (i = 0; i < problem->size; i++)
  {
    for (j = problem->size; j < 2 * problem->size; j++)
    {
      lcp_problem->M->matrix0[i + j * lcp_problem->size] =  0.0;
    }
    lcp_problem->M->matrix0[i + (i + problem->size)*lcp_problem->size] =  1.0;
  }
  for (i = problem->size; i < 2 * problem->size; i++)
  {
    for (j = 0; j < 2 * problem->size; j++)
    {
      lcp_problem->M->matrix0[i + j * lcp_problem->size] =  0.0;
    }
    lcp_problem->M->matrix0[i + (i - problem->size)*lcp_problem->size] =  -1.0;
  }

  for (i = 0; i < problem->size; i++)
  {
    lcp_problem->q[i] = problem->q[i];
    lcp_problem->q[i + problem->size] = problem->ub[i] - problem->lb[i];
    for (j = 0; j < problem->size; j++)
    {
      lcp_problem->q[i] += problem->M->matrix0[i + j * (problem->size)] * problem->lb[j];
    }
  }



}
