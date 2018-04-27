/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#include "NonSmoothDrivers.h"
#include "genericMechanical_test_function.h"
#include "GenericMechanicalProblem.h"
#include "SolverOptions.h"
#include "GenericMechanical_Solvers.h"
int genericMechanical_test_function(FILE * f, SolverOptions * options)
{

  int k, info = -1 ;
  GenericMechanicalProblem* problem = genericMechanical_newFromFile(f);
  double *reaction = (double*)calloc(problem->size, sizeof(double));
  double *velocity = (double*)calloc(problem->size, sizeof(double));
  info = genericMechanical_driver(problem,
                                  reaction , velocity,
                                  options);
  double err = 0;
  GenericMechanical_compute_error(problem, reaction , velocity, options->dparam[0], options, &err);
  printf("\n");
  for (k = 0 ; k < problem->size; k++)
  {
    printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k , reaction[k]);
  }
  printf("\n");

  if (!info)
  {
    if (err > options->dparam[0])
    {
      printf("test failed: err>tol\n");
      printf(" ---> info=%i err=%e and tol=%e\n", info, err, options->dparam[0]);

      return 1;
    }
    else
      printf("test passed: info=%i err=%e and tol=%e\n", info, err, options->dparam[0]);

  }
  else
  {
    printf("test failed.\n");
  }
  free(reaction);
  free(velocity);

  freeGenericMechanicalProblem(problem, NUMERICS_GMP_FREE_MATRIX | NUMERICS_GMP_FREE_GMP);

  return info;

}


