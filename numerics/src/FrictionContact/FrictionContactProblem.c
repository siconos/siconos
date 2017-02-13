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
#include <stdlib.h>
#include <assert.h>
#include "FrictionContactProblem.h"
#include "NumericsMatrix.h"
#include <stdio.h>
#include "numerics_verbose.h"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"

void frictionContact_display(FrictionContactProblem* problem)
{

  assert(problem);
  int n = problem->dimension * problem->numberOfContacts;
  printf("FrictionContact Display :\n-------------\n");
  printf("dimension :%d \n", problem->dimension);
  printf("numberOfContacts:%d \n", problem->numberOfContacts);

  if (problem->M)
  {
    printf("M matrix:\n");
    NM_display(problem->M);
  }
  else
    printf("No M matrix:\n");

  if (problem->q)
  {
    printf("q vector:\n");
    NM_vector_display(problem->q,n);
  }
  else
    printf("No q vector:\n");

  if (problem->mu)
  {
    printf("mu vector:\n");
    NM_vector_display(problem->mu,problem->numberOfContacts);
  }
  else
    printf("No mu vector:\n");

}





int frictionContact_printInFile(FrictionContactProblem*  problem, FILE* file)
{
  if (! problem)
  {
    fprintf(stderr, "Numerics, FrictionContactProblem printInFile failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  int i;

  int d  = problem->dimension;
  fprintf(file, "%d\n", d);
  int nc = problem->numberOfContacts;
  fprintf(file, "%d\n", nc);
  printInFile(problem->M, file);
  for (i = 0; i < problem->M->size1; i++)
  {
    fprintf(file, "%32.24e ", problem->q[i]);
  }
  fprintf(file, "\n");
  for (i = 0; i < nc; i++)
  {
    fprintf(file, "%32.24e ", problem->mu[i]);
  }
  fprintf(file, "\n");
  return 0;
}

int frictionContact_printInFilename(FrictionContactProblem* problem, char* filename)
{
  int info = 0;
  FILE * file = fopen(filename, "w");
  
  if (!file)
  {
    return errno;
  }

  info = frictionContact_printInFile(problem, file);

  fclose(file);
  return info;
}

int frictionContact_newFromFile(FrictionContactProblem* problem, FILE* file)
{
  DEBUG_PRINT("Start -- int frictionContact_newFromFile(FrictionContactProblem* problem, FILE* file)\n");
  int nc = 0, d = 0;
  int i;
  CHECK_IO(fscanf(file, "%d\n", &d));
  problem->dimension = d;
  DEBUG_PRINTF("problem->dimension = %i \n",problem->dimension );
  CHECK_IO(fscanf(file, "%d\n", &nc));
  problem->numberOfContacts = nc;
  problem->M = newNumericsMatrix();

  /* fix: problem->M->storageType unitialized ! */

  newFromFile(problem->M, file);

  problem->q = (double *) malloc(problem->M->size1 * sizeof(double));
  for (i = 0; i < problem->M->size1; i++)
  {
    CHECK_IO(fscanf(file, "%lf ", &(problem->q[i])));
  }

  problem->mu = (double *) malloc(nc * sizeof(double));
  for (i = 0; i < nc; i++)
  {
    CHECK_IO(fscanf(file, "%lf ", &(problem->mu[i])));
  }
  DEBUG_PRINT("End --  int frictionContact_newFromFile(FrictionContactProblem* problem, FILE* file)\n");

  return 0;
}

int frictionContact_newFromFilename(FrictionContactProblem* problem, char* filename)
{
  int info = 0;
  FILE * file = fopen(filename, "r");
  
  if (!file)
  {
    return errno;
  }

  info = frictionContact_newFromFile(problem, file);

  fclose(file);
  return info;
}

void freeFrictionContactProblem(FrictionContactProblem* problem)
{
  assert(problem);
  if (problem->M)
  {
    freeNumericsMatrix(problem->M);
    free(problem->M);
    problem->M = NULL;
  }

  if (problem->mu)
  {
  free(problem->mu);
  problem->mu = NULL;
  }

  if (problem->q)
  {
    free(problem->q);
    problem->q = NULL;
  }

  free(problem);

}

FrictionContactProblem* newFCP(void)
{
  FrictionContactProblem* fcp = (FrictionContactProblem*) malloc(sizeof(FrictionContactProblem));
  fcp->dimension = 0;
  fcp->numberOfContacts = 0;
  fcp->M = NULL;
  fcp->q = NULL;
  fcp->mu = NULL;

  return fcp;
}

FrictionContactProblem* frictionContactProblem_new(int dim, int nc, NumericsMatrix* M, double* q, double* mu)
{
  FrictionContactProblem* fcp = (FrictionContactProblem*) malloc(sizeof(FrictionContactProblem));

  fcp->dimension = dim;
  fcp->numberOfContacts = nc;
  fcp->M = M;
  fcp->q = q;
  fcp->mu = mu;

  return fcp;
}
