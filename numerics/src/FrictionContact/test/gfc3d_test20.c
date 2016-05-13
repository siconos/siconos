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
#include "globalFrictionContact_test_function.h"
#include "fclib_interface.h"



int main(void)
{

  char filename[51] = "./data/LMGC_GlobalFrictionContactProblem00046.hdf5";

  printf("Test on %s\n", filename);

  SolverOptions * options = (SolverOptions *)malloc(sizeof(SolverOptions));

  gfc3d_setDefaultSolverOptions(options, SICONOS_GLOBAL_FRICTION_3D_NSN_AC);

  int k, info = -1 ;

  GlobalFrictionContactProblem* problem = globalFrictionContact_fclib_read(filename);
  globalFrictionContact_display(problem);

  FILE * foutput  =  fopen("checkinput.dat", "w");
  info = globalFrictionContact_printInFile(problem, foutput);

  NumericsOptions global_options;
  setDefaultNumericsOptions(&global_options);
  global_options.verboseMode = 1; // turn verbose mode to off by default

  int NC = problem->numberOfContacts;
  int dim = problem->dimension;
  int n = problem->M->size1;

  double *reaction = (double*)calloc(dim * NC, sizeof(double));
  double *velocity = (double*)calloc(dim * NC, sizeof(double));
  double *globalvelocity = (double*)calloc(n, sizeof(double));

  if (dim == 2)
  {
    info = 1;
  }
  else if (dim == 3)
  {
    info = gfc3d_driver(problem,
                                          reaction , velocity, globalvelocity,
                                          options, &global_options);
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

  deleteSolverOptions(options);
  free(options);

  printf("End of test on %s\n", filename);


  return info;
}
