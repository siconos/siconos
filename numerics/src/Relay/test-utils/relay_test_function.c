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
#include <math.h>
#include "NonSmoothDrivers.h"
#include "relay_test_function.h"
#include "RelayProblem.h"
#include "SolverOptions.h"
#include "Relay_Solvers.h"
#include "SiconosCompat.h"

int relay_test_function(FILE * f, int  solverId)
{

  int i, info = 0 ;
  RelayProblem* problem = (RelayProblem *)malloc(sizeof(RelayProblem));

  info = relay_newFromFile(problem, f);

  FILE * foutput  =  fopen("./relay.verif", "w");
  info = relay_printInFile(problem, foutput);
  SolverOptions * options = (SolverOptions *)malloc(sizeof(SolverOptions));

  relay_setDefaultSolverOptions(problem, options, solverId);

  int maxIter = 50000;
  double tolerance = 1e-8;
  options->iparam[0] = maxIter;
  options->dparam[0] = tolerance;


  double * z = (double *)calloc(problem->size, sizeof(double));
  double * w = (double *)calloc(problem->size, sizeof(double));

  info = relay_driver(problem, z , w, options);

  for (i = 0 ; i < problem->size ; i++)
  {
    printf("z[%i] = %12.8e\t,w[%i] = %12.8e\n", i, z[i], i, w[i]);
    info = info == 0 ? !(isfinite(z[i]) && isfinite(w[i])): info;
  }

  if (!info)
  {
    printf("test succeeded\n");
  }
  else
  {
    printf("test unsuccessful\n");
  }
  free(z);
  free(w);

  solver_options_delete(options);

  free(options);

  freeRelay_problem(problem);

  fclose(foutput);

  return info;


}

