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

#define _XOPEN_SOURCE 700
#include <math.h>                           // for isfinite
#include <stdio.h>                          // for printf, fclose, fopen, FILE
#include <stdlib.h>                         // for calloc, free
#include "NonSmoothDrivers.h"               // for rolling_fc3d_driver
#include "NumericsFwd.h"                    // for RollingFrictionContactPro...
#include "NumericsVerbose.h"                // for numerics_set_verbose
#include "RollingFrictionContactProblem.h"  // for rollingFrictionContactPro...
#include "frictionContact_test_utils.h"     // for rollingFrictionContact_te...
#include "test_utils.h"                     // for TestCase

int rollingFrictionContact_test_function(TestCase* current)
{

  int k;
  RollingFrictionContactProblem* problem = rollingFrictionContact_new_from_filename(current->filename);
  numerics_set_verbose(1);

  FILE * foutput  =  fopen("checkinput.dat", "w");
  rollingFrictionContact_printInFile(problem, foutput);

  int NC = problem->numberOfContacts;
  int dim = problem->dimension;

  int info;
  double *reaction = (double*)calloc(dim * NC, sizeof(double));
  double *velocity = (double*)calloc(dim * NC, sizeof(double));
  printf("\n rollingFrictionContact_test_function run\n");
  if(dim == 2)
  {
    info = 1;
  }
  else if(dim == 5)
  {
    info = rolling_fc3d_driver(problem,
                               reaction, velocity,
                               current->options);
  }
  printf("\n");
  for(k = 0 ; k < dim * NC; k++)
  {
    printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k, reaction[k]);
  }
  printf("\n");

  for(k = 0; k < dim * NC; ++k)
  {
    info = info == 0 ? !(isfinite(velocity[k]) && isfinite(reaction[k])) : info;
  }

  if(!info)
  {
    printf("test succeeded\n");
  }
  else
  {
    printf("test unsuccessful\n");
  }
  free(reaction);
  free(velocity);
  fclose(foutput);

  rollingFrictionContactProblem_free(problem);

  return info;
}


