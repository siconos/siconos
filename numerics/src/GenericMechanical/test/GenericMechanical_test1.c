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


int main(void)
{
  int info = 0 ;

  FILE * finput  =  0;
  SolverOptions * options = (SolverOptions *) malloc(sizeof(SolverOptions));

  genericMechanicalProblem_setDefaultSolverOptions(options, SICONOS_FRICTION_3D_ONECONTACT_QUARTIC);
  options->iparam[0] = 100000000;
  //options->iparam[2]=1;
  printf("Test on ./data/GMP.dat\n");
  finput  =  fopen("./data/GMP.dat", "r");
  info = genericMechanical_test_function(finput, options);
  fclose(finput);



  solver_options_delete(options);
  free(options);
  fclose(finput);
  printf("\nEnd of test on ./data/GMP.dat\n");
  return info;
}
