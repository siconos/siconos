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
#include "soclcp_test_function.h"
#include "SOCLCP_cst.h"
#include "SolverOptions.h"
#include "SOCLCP_Solvers.h"


int main(void)
{
  int info = 0 ;
  printf("Test on ./data/Example1_SOCLCP_SBM.dat\n");

  FILE * finput  =  fopen("./data/Example1_SOCLCP_SBM.dat", "r");
  SolverOptions * options = (SolverOptions *) malloc(sizeof(SolverOptions));
  info = soclcp_setDefaultSolverOptions(options, SICONOS_SOCLCP_NSGS);
  options->dparam[0] = 1e-16;
  options->internalSolvers->solverId = SICONOS_SOCLCP_ProjectionOnCone;

  info = soclcp_test_function(finput, options);

  solver_options_delete(options);
  free(options);
  fclose(finput);
  printf("\nEnd of test on ./data/Example1_SOCLCP_SBM.dat\n");
  return info;
}
