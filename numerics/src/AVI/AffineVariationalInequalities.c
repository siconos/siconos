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
#ifndef AVI_PROBLEM_C
#define AVI_PROBLEM_C


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "NumericsMatrix.h"
#include "SiconosSets.h"
#include "AffineVariationalInequalities.h"
#include "numerics_verbose.h"


void AVI_display(AffineVariationalInequalities* problem)
{

  assert(problem);
  int i, n = problem->size;
  printf("AffineVariationalInequalities Display :\n-------------\n");
  printf("size :%d \n", problem->size);
  if (problem->M)
  {
    printf("M matrix:\n");
    NM_display(problem->M);
  }
  else
    printf("No M matrix.\n");

  if (problem->q)
  {
    printf("q vector:\n");
    for (i = 0; i < n; i++) printf("q[ %i ] = %12.8e\n", i, problem->q[i]);
  }
  else
    printf("No q vector.\n");

  if (problem->d)
  {
    printf("d vector:\n");
    for (i = 0; i < n; i++) printf("d[ %i ] = %12.8e\n", i, problem->d[i]);
  }
  else
    printf("No d vector.\n");

}





int AVI_printInFile(AffineVariationalInequalities*  problem, FILE* file)
{
  if (! problem)
  {
    fprintf(stderr, "Numerics, AffineVariationalInequalities printInFile failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  int i;
  int n = problem->size;
  fprintf(file, "%d\n", n);
  printInFile(problem->M, file);
  for (i = 0; i < problem->M->size1; i++)
  {
    fprintf(file, "%32.24e ", problem->q[i]);
  }
  return 1;
}

int AVI_newFromFile(AffineVariationalInequalities* problem, FILE* file)
{
  int n = 0;
  int i;

  CHECK_IO(fscanf(file, "%d\n", &n));
  problem->size = n;
  problem->M = newNumericsMatrix();

  newFromFile(problem->M, file);

  problem->q = (double *) malloc(problem->M->size1 * sizeof(double));
  for (i = 0; i < problem->M->size1; i++)
  {
    CHECK_IO(fscanf(file, "%lf ", &(problem->q[i])));
  }
  return 1;
}
int AVI_newFromFilename(AffineVariationalInequalities* problem, char* filename)
{
  int info = 0;
  FILE * file = fopen(filename, "r");

  info = AVI_newFromFile(problem, file);

  fclose(file);
  return info;
}

void freeAVI(AffineVariationalInequalities* problem)
{
  freeNumericsMatrix(problem->M);
  free(problem->M);
  free(problem->q);
  if (problem->poly)
  {
    free_polyhedron(problem->poly);
  }
  if (problem->d)
    free(problem->d);
  free(problem);
  problem = NULL;
}



#endif

