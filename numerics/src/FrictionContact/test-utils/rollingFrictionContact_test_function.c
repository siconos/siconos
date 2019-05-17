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

#define _XOPEN_SOURCE 700
#include <string.h>

#if (__linux ||  __APPLE__)
#elif _MSC_VER
#define strdup _strdup
#else
static inline char* strdup(char* src)
{
  size_t len = strlen(src) + 1;
  char* dest = (char*)malloc(len * sizeof(char));
  strcpy(dest, src, len);
  return dest;
}
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "CSparseMatrix_internal.h"

// avoid a conflict with old csparse.h in case fclib includes it
#define _CS_H

#include "NonSmoothDrivers.h"
#include "rollingFrictionContact_test_function.h"
#include "gfc3d_Solvers.h"
#include "RollingFrictionContactProblem.h"
#include "NumericsMatrix.h"
#include "numerics_verbose.h"
#include "NumericsVector.h"
#include "SiconosCompat.h"

#include <string.h>
#if defined(WITH_FCLIB)
#include <fclib.h>
#include <fclib_interface.h>
#endif

#ifdef __cplusplus
using namespace std;
#endif

int rollingFrictionContact_test_function(FILE * f, SolverOptions * options)
{

  int k;
  RollingFrictionContactProblem* problem = (RollingFrictionContactProblem *)malloc(sizeof(RollingFrictionContactProblem));
  /* numerics_set_verbose(1); */

  rollingFrictionContact_newFromFile(problem, f);
  rollingFrictionContact_display(problem);


  FILE * foutput  =  fopen("checkinput.dat", "w");
  rollingFrictionContact_printInFile(problem, foutput);

  int NC = problem->numberOfContacts;
  int dim = problem->dimension;

  int info;
  double *reaction = (double*)malloc(dim * NC * sizeof(double));
  double *velocity = (double*)malloc(dim * NC * sizeof(double));
  for (k = 0 ; k < dim * NC; k++)
  {
    velocity[k] = 0.0;
    reaction[k] = 0.0;
  }
  if (dim == 2)
  {
    info = 1;
  }
  else if (dim == 5)
  {
    info = rolling_fc3d_driver(problem,
                               reaction , velocity,
                               options);
  }
  printf("\n");
  for (k = 0 ; k < dim * NC; k++)
  {
    printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k , reaction[k]);
  }
  printf("\n");

  for (k = 0; k < dim * NC; ++k)
  {
    info = info == 0 ? !(isfinite(velocity[k]) && isfinite(reaction[k])) : info;
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
  fclose(foutput);

  rollingFrictionContactProblem_free(problem);


  return info;

}

#if defined(WITH_FCLIB)

int rfc3d_test_function_hdf5(const char* path, SolverOptions* options)
{
  int info =0;

  return info;

}
#endif
