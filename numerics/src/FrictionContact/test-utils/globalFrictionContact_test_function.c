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

#include "CSparseMatrix.h"

// avoid a conflict with old csparse.h in case fclib includes it
#define _CS_H

#include "NonSmoothDrivers.h"
#include "globalFrictionContact_test_function.h"
#include "gfc3d_Solvers.h"
#include "GlobalFrictionContactProblem.h"
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

int globalFrictionContact_test_function(FILE * f, SolverOptions * options)
{

  int k, info = -1 ;
  /* numerics_set_verbose(1); */
  GlobalFrictionContactProblem* problem = globalFrictionContact_newFromFile(f);
  /* globalFrictionContact_display(problem); */


  FILE * foutput  =  fopen("checkinput.dat", "w");
  info = globalFrictionContact_printInFile(problem, foutput);

  int NC = problem->numberOfContacts;
  int dim = problem->dimension;
  int n = problem->M->size1;


  double *reaction = (double*)malloc(dim * NC * sizeof(double));
  double *velocity = (double*)malloc(dim * NC * sizeof(double));
  double *globalvelocity = (double*)malloc(n * sizeof(double));
  for (k = 0 ; k < dim * NC; k++)
  {
    velocity[k] = 0.0;
    reaction[k] = 0.0;
  }
  for (k = 0 ; k < n; k++)
  {
    globalvelocity[k] = 0.0;
  }
  if (dim == 2)
  {
    info = 1;
  }
  else if (dim == 3)
  {
    info = gfc3d_driver(problem,
			reaction , velocity, globalvelocity,
			options);
  }
  int print_size = 10;

  if  (dim * NC >= print_size)
  {
    printf("First values (%i)\n", print_size);
    for (k = 0 ; k < print_size; k++)
    {
      printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k , reaction[k]);
    }
    printf(" ..... \n");
    for (k = 0 ; k < print_size; k++)
  {
    printf("GlocalVelocity[%i] = %12.8e\n", k, globalvelocity[k]);
  }
  }
  else
  {
    for (k = 0 ; k < dim * NC; k++)
    {
      printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k , reaction[k]);
    }
    printf("\n");
    for (k = 0 ; k < dim*NC; k++)
    {
      printf("GlocalVelocity[%i] = %12.8e\n", k, globalvelocity[k]);
    }
  }
  printf("\n");

  for (k = 0; k < dim * NC; ++k)
  {
    info = info == 0 ? !(isfinite(velocity[k]) && isfinite(reaction[k])) : info;
  }

  for (k = 0; k < n; ++k)
  {
    info = info == 0 ? !(isfinite(globalvelocity[k])) : info;
  }

  if (!info)
  {
    printf("test successful with %i iterations and residual = %e\n", options->iparam[SICONOS_IPARAM_ITER_DONE], options->dparam[SICONOS_DPARAM_RESIDU]);
  }
  else
  {
    printf("test unsuccessful with %i iterations and residual = %e\n", options->iparam[SICONOS_IPARAM_ITER_DONE], options->dparam[SICONOS_DPARAM_RESIDU]);
    //getchar();
  }
  free(reaction);
  free(velocity);
  free(globalvelocity);
  fclose(foutput);

  globalFrictionContact_free(problem);


  return info;

}

#if defined(WITH_FCLIB)

int gfc3d_test_function_hdf5(const char* path, SolverOptions* options)
{

  int k, info = -1 ;
  GlobalFrictionContactProblem* problem = globalFrictionContact_fclib_read(path);

  int check_input=1;
  if(check_input)
  {
    int nLen;
    nLen = strlen (path);


    /* remove the extension */
    char * path_copy = strdup(path);
    printf("path_copy = %s \n", path_copy);


    if ((nLen > 0) && (nLen < 400)) {

      while (nLen) {

           // Check for extension character !!!
           if (path_copy [nLen] == '.') {

                path_copy [nLen] = '\0';
                break;
           }

           nLen --;

      }
      printf("path_copy = %s \n", path_copy);
    }
    free(path_copy);

    char * path_out = (char *)calloc((nLen+10), sizeof(char));
    sprintf(path_out, "%s.dat", path); /* finally we keep the extension .hdf5.dat */
    printf("path_out = %s \n", path_out);
    FILE * foutput  =  fopen(path_out, "w");
    info = globalFrictionContact_printInFile(problem, foutput);
    fclose(foutput);
    free(path_out);
  }




  int NC = problem->numberOfContacts;
  int dim = problem->dimension;
  int n = problem->M->size0;

  double *reaction = (double*)calloc(dim * NC, sizeof(double));
  double *velocity = (double*)calloc(dim * NC, sizeof(double));
  double *global_velocity = (double*)calloc(n, sizeof(double));
  /* verbose=1; */
  if (dim == 3)
  {
    info = gfc3d_driver(problem, reaction, velocity, global_velocity, options);
  }
  else
  {
    fprintf(stderr, "gfc3d_test_function_hdf5 :: problem size != 3\n");
    return 1;
  }
  printf("\n");

  int print_size = 10;

  if  (dim * NC >= print_size)
  {
    printf("First values (%i)\n", print_size);
    for (k = 0 ; k < print_size; k++)
    {
      printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k , reaction[k]);
    }
    printf(" ..... \n");
    for (k = 0 ; k < print_size; k++)
  {
    printf("GlocalVelocity[%i] = %12.8e\n", k, global_velocity[k]);
  }
  }
  else
  {
    for (k = 0 ; k < dim * NC; k++)
    {
      printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k , reaction[k]);
    }
    printf("\n");
  }

  /* for (k = 0 ; k < dim * NC; k++) */
  /* { */
  /*   printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k , reaction[k]); */
  /* } */
  /* printf("\n"); */

  if (!info)
  {
    printf("test successful with %i iterations and residual = %e\n", options->iparam[SICONOS_IPARAM_ITER_DONE], options->dparam[SICONOS_DPARAM_RESIDU]);
  }
  else
  {
    printf("test unsuccessful with %i iterations and residual = %e\n", options->iparam[SICONOS_IPARAM_ITER_DONE], options->dparam[SICONOS_DPARAM_RESIDU]);
    //getchar();
  }
  free(reaction);
  free(velocity);
  free(global_velocity);

  globalFrictionContact_free(problem);


  return info;

}
#endif
