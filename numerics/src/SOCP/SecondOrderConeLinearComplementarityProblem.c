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
#include "SecondOrderConeLinearComplementarityProblem.h"
#include "numerics_verbose.h"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"

void secondOrderConeLinearComplementarityProblem_display(SecondOrderConeLinearComplementarityProblem* problem)
{

  assert(problem);



  int i, n ;

  n=problem->n;

  printf("SecondOrderConeLinearComplementarityProblem Display :\n-------------\n");
  printf("nc:%d \n", problem->nc);

  if(problem->M)
  {
    printf("M matrix:\n");
    NM_display(problem->M);
  }
  else
    printf("No M matrix:\n");

  if(problem->q)
  {
    printf("q vector:\n");
    for(i = 0; i < n; i++) printf("q[ %i ] = %12.8e\n", i, problem->q[i]);
  }
  else
    printf("No q vector:\n");

  if(problem->coneIndex)
  {
    printf("coneIndex vector:\n");
    for(i = 0; i < problem->nc+1; i++) printf("coneIndex[ %i ] = %i\n", i, problem->coneIndex[i]);
  }
  else
    printf("No mu vector:\n");
  if(problem->tau)
  {
    printf("tau vector:\n");
    for(i = 0; i < problem->nc; i++) printf("mu[ %i ] = %12.8e\n", i, problem->tau[i]);
  }
  else
    printf("No tau vector:\n");

}





int secondOrderConeLinearComplementarityProblem_printInFile(SecondOrderConeLinearComplementarityProblem*  problem, FILE* file)
{
  if(! problem)
  {
    fprintf(stderr, "Numerics, SecondOrderConeLinearComplementarityProblem_printInFile failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  int i;
  int n = problem->n;
  fprintf(file, "%d\n", n);
  int nc = problem->nc;
  fprintf(file, "%d\n", nc);
  NM_write_in_file(problem->M, file);
  for(i = 0; i < problem->M->size1; i++)
  {
    fprintf(file, "%32.24e ", problem->q[i]);
  }
  fprintf(file, "\n");
  for(i = 0; i < nc+1; i++)
  {
    fprintf(file, "%i ", problem->coneIndex[i]);
  }
  fprintf(file, "\n");

  for(i = 0; i < nc; i++)
  {
    fprintf(file, "%32.24e ", problem->tau[i]);
  }
  fprintf(file, "\n");
  return 0;
}

int secondOrderConeLinearComplementarityProblem_printInFilename(SecondOrderConeLinearComplementarityProblem* problem, char* filename)
{
  int info = 0;
  FILE * file = fopen(filename, "w");

  if(!file)
  {
    return errno;
  }

  info = secondOrderConeLinearComplementarityProblem_printInFile(problem, file);

  fclose(file);
  return info;
}

int secondOrderConeLinearComplementarityProblem_newFromFile(SecondOrderConeLinearComplementarityProblem* problem, FILE* file)
{
  DEBUG_PRINT("Start -- int secondOrderConeLinearComplementarityProblem_newFromFile(SecondOrderConeLinearComplementarityProblem* problem, FILE* file)\n");
  int n = 0;
  int nc = 0;
  int i;

  CHECK_IO(fscanf(file, "%d\n", &n));
  problem->n = n;
  CHECK_IO(fscanf(file, "%d\n", &nc));
  problem->nc = nc;
  problem->M = NM_new();

  /* fix: problem->M->storageType unitialized ! */

  NM_new_from_file(problem->M, file);

  problem->q = (double *) malloc(problem->M->size1 * sizeof(double));
  for(i = 0; i < problem->M->size1; i++)
  {
    CHECK_IO(fscanf(file, "%lf ", &(problem->q[i])));
  }
  problem->coneIndex = (unsigned int *) malloc((nc+1) * sizeof(unsigned int));
  for(i = 0; i < nc+1; i++)
  {
    CHECK_IO(fscanf(file, "%d ", &(problem->coneIndex[i])));
  }

  problem->tau = (double *) malloc(nc * sizeof(double));
  for(i = 0; i < nc; i++)
  {
    CHECK_IO(fscanf(file, "%lf ", &(problem->tau[i])));
  }
  DEBUG_PRINT("End --  int secondOrderConeLinearComplementarityProblem_newFromFile(SecondOrderConeLinearComplementarityProblem* problem, FILE* file)\n");

  return 0;
}

int secondOrderConeLinearComplementarityProblem_newFromFilename(SecondOrderConeLinearComplementarityProblem* problem, char* filename)
{
  int info = 0;
  FILE * file = fopen(filename, "r");

  if(!file)
  {
    return errno;
  }

  info = secondOrderConeLinearComplementarityProblem_newFromFile(problem, file);

  fclose(file);
  return info;
}

void freeSecondOrderConeLinearComplementarityProblem(SecondOrderConeLinearComplementarityProblem* problem)
{

  if (problem->M)
  {
    NM_free(problem->M);
    free(problem->M);
    problem->M = NULL;
  }
  if (problem->tau)
  {
    free(problem->tau);
    problem->tau = NULL;
  }
  if (problem->q)
  {
    free(problem->q);
    problem->q = NULL;
  }
  if (problem->coneIndex)
  {
    free(problem->coneIndex);
    problem->coneIndex = NULL;
  }
  free(problem);

}

SecondOrderConeLinearComplementarityProblem* secondOrderConeLinearComplementarityProblem_new(int n,  int nc, NumericsMatrix* M, double* q, unsigned int *coneIndex, double* tau)
{
  SecondOrderConeLinearComplementarityProblem* soclcp = (SecondOrderConeLinearComplementarityProblem*) malloc(sizeof(SecondOrderConeLinearComplementarityProblem));

  soclcp->n = n;
  soclcp->nc = nc;
  soclcp->M = M;
  soclcp->q = q;
  soclcp->coneIndex = coneIndex;
  soclcp->tau = tau;

  return soclcp;
}
