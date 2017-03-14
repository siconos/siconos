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
#include "NonSmoothDrivers.h"
#include "globalFrictionContact_test_function.h"
#include "gfc3d_Solvers.h"
#include "GlobalFrictionContactProblem.h"
#include "NumericsMatrix.h"
#include "numerics_verbose.h"

int globalFrictionContact_test_function(FILE * f, SolverOptions * options)
{

  int k, info = -1 ;
  GlobalFrictionContactProblem* problem = (GlobalFrictionContactProblem *)malloc(sizeof(GlobalFrictionContactProblem));
  numerics_set_verbose(1);

  info = globalFrictionContact_newFromFile(problem, f);
  globalFrictionContact_display(problem);


  FILE * foutput  =  fopen("checkinput.dat", "w");
  info = globalFrictionContact_printInFile(problem, foutput);

  int NC = problem->numberOfContacts;
  int dim = problem->dimension;
  int n = problem->M->size1;


  double *reaction = (double*)malloc(dim * NC * sizeof(double));
  double *velocity = (double*)malloc(dim * NC * sizeof(double));
  double *globalvelocity = (double*)malloc(n * sizeof(double));
  for (k = 0 ; k < dim * NC; k++)
  {
    velocity[k] = 0.0;
    reaction[k] = 0.0;
  }
  for (k = 0 ; k < n; k++)
  {
    globalvelocity[k] = 0.0;
  }

  if (dim == 2)
  {
    info = 1;
  }
  else if (dim == 3)
  {
    info = gfc3d_driver(problem,
			reaction , velocity, globalvelocity,
			options);
  }
  printf("\n");
  for (k = 0 ; k < dim * NC; k++)
  {
    printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k , reaction[k]);
  }
  for (k = 0 ; k < n; k++)
  {
    printf("GlocalVelocity[%i] = %12.8e\n", k, globalvelocity[k]);
  }
  printf("\n");

  for (k = 0; k < dim * NC; ++k)
  {
    info = info == 0 ? !(isfinite(velocity[k]) && isfinite(reaction[k])) : info;
  }

  for (k = 0; k < n; ++k)
  {
    info = info == 0 ? !(isfinite(globalvelocity[k])) : info;
  }

  if (!info)
  {
    printf("test succeeded\n");
  }
  else
  {
    printf("test unsuccessful\n");
  }
  free(reaction);
  free(velocity);
  free(globalvelocity);
  fclose(foutput);

  freeGlobalFrictionContactProblem(problem);


  return info;

}


