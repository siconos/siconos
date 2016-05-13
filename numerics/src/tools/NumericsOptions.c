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
#include "NumericsOptions.h"

/* Default value for verbose mode: turned to off
Warning: global variable
*/
int verbose = 0;
void setNumericsVerbose(int newVerboseMode)
{
  verbose = newVerboseMode;
}

void setNumericsOptions(NumericsOptions* opt)
{
  verbose = opt->verboseMode;
}

void numericsError(char * functionName, char* message)
{
  char output[200] = "Numerics error - ";
  strcat(output, functionName);
  strcat(output, message);
  strcat(output, ".\n");
  fprintf(stderr, "%s", output);
  exit(EXIT_FAILURE);
}

void numericsWarning(char * functionName, char* message)
{
  char output[200] = "Numerics warning - ";
  strcat(output, functionName);
  strcat(output, message);
  strcat(output, ".\n");
  fprintf(stderr, "%s", output);
  exit(EXIT_FAILURE);
}

void setDefaultNumericsOptions(NumericsOptions* opts)
{
  opts->verboseMode = 0;
}
