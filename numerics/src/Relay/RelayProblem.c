/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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


#include "RelayProblem.h"
#include <assert.h>            // for assert
#include <stdio.h>             // for printf, fprintf, fscanf, FILE, stderr
#include <stdlib.h>            // for free, malloc, exit, EXIT_FAILURE
#include "NumericsMatrix.h"    // for NumericsMatrix, NM_display, NM_clear
#include "numerics_verbose.h"  // for CHECK_IO

void Relay_display(RelayProblem* p)
{

  assert(p);
  int i, n = p->size;
  printf("RelayProblem Display :\n-------------\n");
  printf("size :%d \n", p->size);
  if(p->M)
  {
    printf("M matrix:\n");
    NM_display(p->M);
  }
  else
    printf("No M matrix:\n");

  if(p->q)
  {
    printf("q vector:\n");
    for(i = 0; i < n; i++) printf("q[ %i ] = %12.8e\n", i, p->q[i]);
  }
  else
    printf("No q vector:\n");

  if(p->lb)
  {
    printf("lb vector:\n");
    for(i = 0; i < n; i++) printf("lb[ %i ] = %12.8e\n", i, p->lb[i]);
  }
  else
    printf("No lb vector:\n");

  if(p->ub)
  {
    printf("ub vector:\n");
    for(i = 0; i < n; i++) printf("ub[ %i ] = %12.8e\n", i, p->ub[i]);
  }
  else
    printf("No ub vector:\n");
}




int relay_printInFile(RelayProblem*  problem, FILE* file)
{
  if(! problem)
  {
    fprintf(stderr, "Numerics, RelayProblem printInFile failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  int i;
  int n = problem->size;
  fprintf(file, "%d\n", n);
  NM_write_in_file(problem->M, file);
  for(i = 0; i < problem->M->size1; i++)
  {
    fprintf(file, "%32.24e ", problem->q[i]);
  }
  fprintf(file, "\n");
  for(i = 0; i < problem->M->size1; i++)
  {
    fprintf(file, "%32.24e ", problem->lb[i]);
  }
  fprintf(file, "\n");
  for(i = 0; i < problem->M->size1; i++)
  {
    fprintf(file, "%32.24e ", problem->ub[i]);
  }
  return 1;
}

RelayProblem* relayProblem_new(void)
{
  RelayProblem* rp = (RelayProblem*) malloc(sizeof(RelayProblem));
  rp->size = 0;
  rp->M = NULL;
  rp->q = NULL;
  rp->lb = NULL;
  rp->ub = NULL;

  return rp;
}


RelayProblem* relay_newFromFile(FILE* file)
{
  RelayProblem* problem = relayProblem_new();

  int n = 0;
  int i;

  CHECK_IO(fscanf(file, "%d\n", &n));
  problem->size = n;
  problem->M =  NM_new_from_file(file);

  problem->q = (double *) malloc(problem->M->size1 * sizeof(double));
  for(i = 0; i < problem->M->size1; i++)
  {
    CHECK_IO(fscanf(file, "%lf ", &(problem->q[i])));
  }

  problem->lb = (double *) malloc(problem->M->size1 * sizeof(double));
  for(i = 0; i < problem->M->size1; i++)
  {
    CHECK_IO(fscanf(file, "%lf ", &(problem->lb[i])));
  }

  problem->ub = (double *) malloc(problem->M->size1 * sizeof(double));
  for(i = 0; i < problem->M->size1; i++)
  {
    CHECK_IO(fscanf(file, "%lf ", &(problem->ub[i])));
  }
  return problem;
}

RelayProblem * relay_new_from_filename(const char* filename)
{
  RelayProblem* problem = NULL;

  FILE * file = fopen(filename, "r");
  if(file == NULL)
    numerics_error("RelayProblem", "Can not open file ", filename);

  problem = relay_newFromFile(file);

  fclose(file);
  return problem;
}


void freeRelay_problem(RelayProblem* problem)
{
  assert(problem);
  if(problem->M)
  {
    NM_clear(problem->M);
    free(problem->M);
  }
  if(problem->q)
  {
    free(problem->q);
  }
  if(problem->lb)
  {
    free(problem->lb);
  }
  if(problem->ub)
  {
    free(problem->ub);
  }
  free(problem);
}




#endif

