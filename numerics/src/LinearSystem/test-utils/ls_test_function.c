/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.

 * Copyright 2016 INRIA.

 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at

 * http://www.apache.org/licenses/LICENSE-2.0

 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#include <stdio.h>
#include <stdlib.h>
#include "NonSmoothDrivers.h"
#include "ls_test_function.h"


void _LSfillParamWithRespectToSolver(SolverOptions *options, int solverId, LinearSystemProblem* problem)
{
  options->solverId = solverId;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 5;
  options->dSize = 5;
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  /*use dgels ?*/
  options->iparam[4] = 0;
  if (problem)
  {
    options->dWork = (double*) malloc(LinearSystem_getNbDwork(problem, options) * sizeof(double));
    options->iWork = (int*) malloc(LinearSystem_getNbIwork(problem, options) * sizeof(int));
  }
  else
  {
    options->dWork = NULL;
    options->iWork = NULL;
  }
  options->dparam[0] = 1e-12;
}


void _LSfillParamWithRespectToSolver_SBM(SolverOptions *options, int solverId, LinearSystemProblem* problem)
{


}



int ls_test_function(FILE * f, int solverId)
{

  int i, info = 0 ;
  LinearSystemProblem* problem = (LinearSystemProblem *)malloc(sizeof(LinearSystemProblem));

  info = LinearSystem_newFromFile(problem, f);

  SolverOptions * options ;
  options = (SolverOptions *)malloc(sizeof(*options));

  options->solverId = solverId;
  printf("solverName ==> %s\n", idToName(solverId));

  _LSfillParamWithRespectToSolver(options, solverId, problem);

  options->filterOn = 1;
  double * z = (double *)malloc(problem->size * sizeof(double));
  double * w = (double *)malloc(problem->size * sizeof(double));
  for (i = 0; i < problem->size; i++)
  {
    z[i] = 0.0;
    w[i] = 0.0;
  }



  info = LinearSystem_driver(problem, z , w, options);

  for (i = 0 ; i < problem->size ; i++)
  {
    printf("z[%i] = %12.8e\t,w[%i] = %12.8e\n", i, z[i], i, w[i]);
  }

  if (!info)
  {
    printf("test succeeded err = %e \n", options->dparam[1]);
  }
  else
  {
    printf("test unsuccessful err =%e  \n", options->dparam[1]);
  }
  free(z);
  free(w);



  free(options->dWork);
  free(options->iWork);

  free(options->iparam);
  free(options->dparam);


  free(options);

  LinearSystem_freeProblem(problem);

  printf("End of test.\n");

  return info;


}

int ls_test_function_SBM(FILE * f, int solverId)
{

  int info = -1;
  return info;


}


