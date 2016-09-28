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
#include "relay_test_function.h"
#include "SolverOptions.h"

int main(void)
{
  int info = 0 ;
  char filename[50] = "./data/step_1x1.dat";

  printf("Test on %s\n", filename);

  FILE * finput  =  fopen(filename, "r");

  char solvername[20] = "RELAY_LEMKE";
  info = relay_test_function(finput, solver_options_name_to_id(solvername));

  fclose(finput);

  printf("End of test on %s\n", filename);

  return info;
}

