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
#include "frictionContact_test_function.h"



int main(void)
{
  int info = 0 ;
  printf("Test on ./data/Example1_Fc3D_SBM.dat\n");

  FILE * finput  =  fopen("./data/Example1_Fc3D_SBM.dat", "r");
  SolverOptions * options = (SolverOptions *) malloc(sizeof(SolverOptions));
  info = fc3d_setDefaultSolverOptions(options, SICONOS_FRICTION_3D_NSGS);
  options->dparam[0] = 1e-16;
  options->internalSolvers->solverId = SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration;
  options->internalSolvers->iparam[0] = 10;
  options->internalSolvers->dparam[0] = 1e-3;
  
  options->iparam[8] = 1;
  options->dparam[8] = 1.70;

  info = frictionContact_test_function(finput, options);

  deleteSolverOptions(options);
  free(options);
  fclose(finput);
  printf("\nEnd of test on ./data/Example1_Fc3D_SBM.dat\n");
  return info;
}
