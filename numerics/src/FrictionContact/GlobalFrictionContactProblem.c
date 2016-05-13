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
#include <errno.h>

#include "GlobalFrictionContactProblem.h"
#include "misc.h"

int globalFrictionContact_printInFile(GlobalFrictionContactProblem*  problem, FILE* file)
{
  if (! problem)
  {
    fprintf(stderr, "Numerics, GlobalFrictionContactProblem printInFile failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  int i;

  int d  = problem->dimension;
  fprintf(file, "%d\n", d);
  int nc = problem->numberOfContacts;
  fprintf(file, "%d\n", nc);
  printInFile(problem->M, file);
  printInFile(problem->H, file);
  for (i = 0; i < problem->M->size1; i++)
  {
    fprintf(file, "%32.24e ", problem->q[i]);
  }
  fprintf(file, "\n");
  for (i = 0; i < problem->H->size1; i++)
  {
    fprintf(file, "%32.24e ", problem->b[i]);
  }
  fprintf(file, "\n");
  for (i = 0; i < nc; i++)
  {
    fprintf(file, "%32.24e ", problem->mu[i]);
  }
  fprintf(file, "\n");
  return 0;
}

int globalFrictionContact_newFromFile(GlobalFrictionContactProblem* problem, FILE* file)
{
  int nc = 0, d = 0;
  int info = 0;
  CHECK_IO(fscanf(file, "%d\n", &d), &info);
  problem->dimension = d;
  CHECK_IO(fscanf(file, "%d\n", &nc), &info);
  problem->numberOfContacts = nc;
  problem->M = newNumericsMatrix();

  info = newFromFile(problem->M, file);
  if (info) goto fail;

  problem->H = newNumericsMatrix();
  info = newFromFile(problem->H, file);
  if (info) goto fail;

  problem->q = (double *) malloc(problem->M->size1 * sizeof(double));
  for (int i = 0; i < problem->M->size1; ++i)
  {
    CHECK_IO(fscanf(file, "%lf ", &(problem->q[i])), &info);
  }
  problem->b = (double *) malloc(problem->H->size1 * sizeof(double));
  for (int i = 0; i < problem->H->size1; ++i)
  {
    CHECK_IO(fscanf(file, "%lf ", &(problem->b[i])), &info);
  }

  problem->mu = (double *) malloc(nc * sizeof(double));
  for (int i = 0; i < nc; ++i)
  {
    CHECK_IO(fscanf(file, "%lf ", &(problem->mu[i])), &info);
  }

fail:
  problem->env = NULL;
  problem->workspace = NULL;
  return info;
}

void freeGlobalFrictionContactProblem(GlobalFrictionContactProblem* problem)
{

  if (problem->M)
  {
    freeNumericsMatrix(problem->M);
    free(problem->M);
    problem->M = NULL;
  }

  if (problem->H)
  {
    freeNumericsMatrix(problem->H);
    free(problem->H);
    problem->H = NULL;
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

  if (problem->b)
  {
    free(problem->b);
    problem->b = NULL;
  }

  if (problem->env) assert(0 && "freeGlobalFrictionContactProblem :: problem->env != NULL, don't know what to do");

  gfc3d_free_workspace(problem);

  free(problem);

}
void globalFrictionContact_display(GlobalFrictionContactProblem* problem)
{

  assert(problem);
  int i, n = problem->dimension * problem->numberOfContacts;
  printf("GlobalFrictionContact Display :\n-------------\n");
  printf("dimension :%d \n", problem->dimension);
  printf("numberOfContacts:%d \n", problem->numberOfContacts);
  int m = problem->M->size0;
  if (problem->M)
  {
    printf("M matrix:\n");
    display(problem->M);
  }
  else
    printf("No M matrix:\n");
  if (problem->H)
  {
    printf("H matrix:\n");
    display(problem->H);
  }
  else
    printf("No H matrix:\n");

  if (problem->q)
  {
    printf("q vector:\n");
    for (i = 0; i < m; i++) printf("q[ %i ] = %12.8e\n", i, problem->q[i]);
  }
  else
    printf("No q vector:\n");

  if (problem->b)
  {
    printf("b vector:\n");
    for (i = 0; i < n; i++) printf("b[ %i ] = %12.8e\n", i, problem->b[i]);
  }
  else
    printf("No q vector:\n");

  if (problem->mu)
  {
    printf("mu vector:\n");
    for (i = 0; i < problem->numberOfContacts; i++) printf("mu[ %i ] = %12.8e\n", i, problem->mu[i]);
  }
  else
    printf("No mu vector:\n");

}

void gfc3d_init_workspace(GlobalFrictionContactProblem* problem)
{
  assert(problem);
  assert(problem->M);

  if (!problem->workspace)
  {
    problem->workspace = (GFC3D_workspace*) malloc(sizeof(GFC3D_workspace));
    problem->workspace->factorized_M = NULL;
    problem->workspace->globalVelocity = NULL;
  }

  if (!problem->workspace->factorized_M)
  {
    problem->workspace->factorized_M = createNumericsMatrix(problem->M->storageType,
                                               problem->M->size0,
                                               problem->M->size1);
    NM_copy(problem->M, problem->workspace->factorized_M);
  }

  if (!problem->workspace->globalVelocity)
  {
    problem->workspace->globalVelocity = (double*)malloc(problem->M->size1 * sizeof(double));
  }
}

void gfc3d_free_workspace(GlobalFrictionContactProblem* problem)
{
  if (problem->workspace)
  {
    if (problem->workspace->factorized_M)
    {
      freeNumericsMatrix(problem->workspace->factorized_M);
      free(problem->workspace->factorized_M);
      problem->workspace->factorized_M = NULL;
    }

    if (problem->workspace->globalVelocity)
    {
      free(problem->workspace->globalVelocity);
      problem->workspace->globalVelocity = NULL;
    }

    free(problem->workspace);
    problem->workspace = NULL;
  }
}
