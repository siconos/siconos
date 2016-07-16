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
#include "frictionContact_test_function.h"


int main(void)
{
  int info = 0 ;
  printf("Test on ./data/BoxesStack1-i100000-32.hdf5.dat\n");

  FILE * finput  =  fopen("./data/BoxesStack1-i100000-32.hdf5.dat", "r");
  SolverOptions * options = (SolverOptions *) malloc(sizeof(SolverOptions));
  info = fc3d_setDefaultSolverOptions(options, SICONOS_FRICTION_3D_NSGS_OPENMP);
  options->dparam[0] = 1e-04;
  options->iparam[0] = 10000; 
  options->iparam[10] = 1; //number of threads 
  options->iparam[11] = 2; // methods u
  options->iparam[12] = 5; // 
  options->iparam[13] = 10; // iteration of interface problem
  options->internalSolvers->solverId = SICONOS_FRICTION_3D_ONECONTACT_NSN_AC;
  options->internalSolvers->iparam[0]=10;
  info = frictionContact_test_function(finput, options);

  deleteSolverOptions(options);
  free(options);
  fclose(finput);
  printf("\nEnd of test on ./data/BoxesStack1-i100000-32.hdf5.dat\n");
  return info;
}
