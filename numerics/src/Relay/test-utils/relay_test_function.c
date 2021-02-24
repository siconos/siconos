/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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
#include <math.h>              // for isfinite
#include <stdio.h>             // for printf, fclose, fopen, FILE
#include <stdlib.h>            // for calloc, free, malloc
#include "NonSmoothDrivers.h"  // for relay_driver
#include "NumericsFwd.h"       // for RelayProblem, SolverOptions
#include "RelayProblem.h"      // for RelayProblem, freeRelay_problem, relay...
#include "SolverOptions.h"     // for solver_options_delete, SolverOptions
#include "relay_test_utils.h"  // for relay_test_function
#include "test_utils.h"        // for TestCase

int relay_test_function(TestCase * current)
{
  int i, info = 0 ;
  RelayProblem* problem = relay_new_from_filename(current->filename);

  FILE * foutput  =  fopen("./relay.verif", "w");
  info = relay_printInFile(problem, foutput);
  fclose(foutput);

  int maxIter = 50000;
  double tolerance = 1e-8;
  current->options->iparam[SICONOS_IPARAM_MAX_ITER] = maxIter;
  current->options->dparam[SICONOS_DPARAM_TOL] = tolerance;


  double * z = (double *)calloc(problem->size, sizeof(double));
  double * w = (double *)calloc(problem->size, sizeof(double));

  info = relay_driver(problem, z, w, current->options);

  for(i = 0 ; i < problem->size ; i++)
  {
    printf("z[%i] = %12.8e\t,w[%i] = %12.8e\n", i, z[i], i, w[i]);
    info = info == 0 ? !(isfinite(z[i]) && isfinite(w[i])): info;
  }

  if(!info)
  {
    printf("test succeeded\n");
  }
  else
  {
    printf("test unsuccessful\n");
  }
  free(z);
  free(w);
  freeRelay_problem(problem);

  return info;
}

