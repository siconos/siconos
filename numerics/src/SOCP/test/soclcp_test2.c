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
#include <stdio.h>                 // for printf, fclose, fopen, FILE
#include <stdlib.h>                // for free, malloc
#include "NumericsFwd.h"           // for SolverOptions
#include "SOCLCP_Solvers.h"        // for soclcp_setDefaultSolverOptions
#include "SOCLCP_cst.h"            // for SICONOS_SOCLCP_NSGS, SICONOS_SOCLC...
#include "SolverOptions.h"         // for SolverOptions, solver_options_delete
#include "soclcp_test_function.h"  // for soclcp_test_function

int main(void)
{
  int info = 0 ;
  printf("Test on ./data/Capsules-i122-1617.dat \n");

  FILE * finput  =  fopen("./data/Capsules-i122-1617.dat", "r");
  SolverOptions * options = (SolverOptions *) malloc(sizeof(SolverOptions));
  info = soclcp_setDefaultSolverOptions(options, SICONOS_SOCLCP_NSGS);
  options->dparam[0] = 1e-06;
  options->iparam[0] = 2000000;
  //options->iparam[8] = 1;
  //options->dparam[8] = 1.0;
  options->internalSolvers->solverId = SICONOS_SOCLCP_ProjectionOnConeWithLocalIteration;
  options->internalSolvers->dparam[0]=1e-16;
  options->internalSolvers->iparam[0]=100;
  info = soclcp_test_function(finput, options);

  solver_options_delete(options);
  free(options);
  fclose(finput);
  printf("\nEnd of test on ./data/Capsules-i122-1617.dat \n");
  return info;
}
