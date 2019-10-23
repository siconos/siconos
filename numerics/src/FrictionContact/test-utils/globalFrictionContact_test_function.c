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

#include "globalFrictionContact_test_function.h"
#include <math.h>                          // for isfinite
#include <stdio.h>                         // for printf, fclose, fopen, FILE
#include <stdlib.h>                        // for malloc, free
#include "GlobalFrictionContactProblem.h"  // for GlobalFrictionContactProblem
#include "gfc3d_Solvers.h"  // for GlobalFrictionContactProblem
#include "NonSmoothDrivers.h"              // for gfc3d_driver
#include "NumericsFwd.h"                   // for GlobalFrictionContactProblem
#include "NumericsMatrix.h"                // for NumericsMatrix
#include "SolverOptions.h"                 // for solver_options_delete, sol...
#include "test_utils.h"                    // for TestCase
#include "frictionContact_test_utils.h"
#if defined(WITH_FCLIB)
#include <fclib_interface.h>
#include <string.h>
#endif

int globalFrictionContact_test_function(TestCase* current)
{

  int k, info = -1 ;
  GlobalFrictionContactProblem* problem = globalFrictionContact_new_from_filename(current->filename);
  globalFrictionContact_display(problem);


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
			current->options);
  }
  /* printf("\n"); */
  /* for (k = 0 ; k < dim * NC; k++) */
  /* { */
  /*   printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k , reaction[k]); */
  /* } */
  /* for (k = 0 ; k < n; k++) */
  /* { */
  /*   printf("GlocalVelocity[%i] = %12.8e\n", k, globalvelocity[k]); */
  /* } */
  /* printf("\n"); */

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
  solver_options_delete(current->options);
  solver_options_nullify(current->options);
  freeGlobalFrictionContactProblem(problem);


  return info;

}

#if defined(WITH_FCLIB)

int gfc3d_test_function_hdf5(const char* path, SolverOptions* options)
{

  int k, info = -1 ;
  // Read gfc problem from an hdf5 file, using fclib interface.
  GlobalFrictionContactProblem* problem = globalFrictionContact_fclib_read(path);

  int check_input=1;
  if(check_input)
  {
    int nLen;
    nLen = strlen (path);

    /* remove the extension */
    char * path_copy = (char *) malloc(400*sizeof(char));;
    strcpy(path_copy, path);
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
      free(path_copy);
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
    printf("test successful, residual = %g\n", options->dparam[1]);
  }
  else
  {
    printf("test unsuccessful, residual = %g\n", options->dparam[1]);
  }
  free(reaction);
  free(velocity);
  free(global_velocity);

  freeGlobalFrictionContactProblem(problem);


  return info;

}
#endif

void build_gfc3d_test(const char * filename,
                int solver_id, int* d_ind, double* dparam, int * i_ind, int* iparam,
                int internal_solver_id, int * i_d_ind, double * internal_dparam, int * i_i_ind, int * internal_iparam,
                TestCase* testname)
{
  // reference file name
  testname->filename = filename;
  // By default, test is expected to succeed.
  testname->will_fail = 0;

  // Set solver options to default.
  testname->options = (SolverOptions *)malloc(sizeof(SolverOptions));
  gfc3d_setDefaultSolverOptions(testname->options, solver_id);
  // Fill iparam and dparam in.
  if(iparam)
    for(int i=0; i<i_ind[0]; ++i)
      testname->options->iparam[i_ind[i+1]] = iparam[i];

  // dparam
  if(dparam)
    for(int i=0; i<d_ind[0]; ++i)
      testname->options->dparam[d_ind[i+1]] = dparam[i];

  // Internal solver setup
  if(internal_solver_id>0)
    {
      testname->options->internalSolvers[0].solverId=internal_solver_id;
      // internal iparam
      if(internal_iparam)
        for(int i=0; i<i_i_ind[0]; ++i)
          testname->options->internalSolvers[0].iparam[i_i_ind[i+1]] = internal_iparam[i];
      // internal dparam
      if(internal_dparam)
        for(int i=0; i<i_d_ind[0]; ++i)
          testname->options->internalSolvers[0].dparam[i_d_ind[i+1]] = internal_dparam[i];
    }
}



