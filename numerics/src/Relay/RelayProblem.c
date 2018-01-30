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
#ifndef RELAY_PROBLEM_C
#define RELAY_PROBLEM_C


#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "RelayProblem.h"
#include "NumericsMatrix.h"
#include "numerics_verbose.h"

void Relay_display(RelayProblem* p)
{

  assert(p);
  int i, n = p->size;
  printf("RelayProblem Display :\n-------------\n");
  printf("size :%d \n", p->size);
  if (p->M)
  {
    printf("M matrix:\n");
    NM_display(p->M);
  }
  else
    printf("No M matrix:\n");

  if (p->q)
  {
    printf("q vector:\n");
    for (i = 0; i < n; i++) printf("q[ %i ] = %12.8e\n", i, p->q[i]);
  }
  else
    printf("No q vector:\n");

  if (p->lb)
  {
    printf("lb vector:\n");
    for (i = 0; i < n; i++) printf("lb[ %i ] = %12.8e\n", i, p->lb[i]);
  }
  else
    printf("No lb vector:\n");

  if (p->ub)
  {
    printf("ub vector:\n");
    for (i = 0; i < n; i++) printf("ub[ %i ] = %12.8e\n", i, p->ub[i]);
  }
  else
    printf("No ub vector:\n");
}




int relay_printInFile(RelayProblem*  problem, FILE* file)
{
  if (! problem)
  {
    fprintf(stderr, "Numerics, RelayProblem printInFile failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  int i;
  int n = problem->size;
  fprintf(file, "%d\n", n);
  NM_write_in_file(problem->M, file);
  for (i = 0; i < problem->M->size1; i++)
  {
    fprintf(file, "%32.24e ", problem->q[i]);
  }
  fprintf(file, "\n");
  for (i = 0; i < problem->M->size1; i++)
  {
    fprintf(file, "%32.24e ", problem->lb[i]);
  }
  fprintf(file, "\n");
  for (i = 0; i < problem->M->size1; i++)
  {
    fprintf(file, "%32.24e ", problem->ub[i]);
  }
  return 1;
}

int relay_newFromFile(RelayProblem* problem, FILE* file)
{
  int n = 0;
  int i;

  CHECK_IO(fscanf(file, "%d\n", &n));
  problem->size = n;
  problem->M =  NM_new_from_file(file);

  problem->q = (double *) malloc(problem->M->size1 * sizeof(double));
  for (i = 0; i < problem->M->size1; i++)
  {
    CHECK_IO(fscanf(file, "%lf ", &(problem->q[i])));
  }

  problem->lb = (double *) malloc(problem->M->size1 * sizeof(double));
  for (i = 0; i < problem->M->size1; i++)
  {
    CHECK_IO(fscanf(file, "%lf ", &(problem->lb[i])));
  }

  problem->ub = (double *) malloc(problem->M->size1 * sizeof(double));
  for (i = 0; i < problem->M->size1; i++)
  {
    CHECK_IO(fscanf(file, "%lf ", &(problem->ub[i])));
  }
  return 1;
}

void freeRelay_problem(RelayProblem* problem)
{
  assert(problem);
  if (problem->M)
  {
    NM_free(problem->M);
    free(problem->M);
  }
  if (problem->q)  { free(problem->q); }
  if (problem->lb) { free(problem->lb); }
  if (problem->ub) { free(problem->ub); }
  free(problem);
}



#endif

