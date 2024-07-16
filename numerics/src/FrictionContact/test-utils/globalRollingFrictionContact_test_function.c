/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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
#include <math.h>                                 // for isfinite
#include <stdio.h>                                // for printf, fclose, fopen, FILE
#include <stdlib.h>                               // for calloc, free
#include "GlobalRollingFrictionContactProblem.h"  // for globalRollingFrictionContactPro...
#include "NonSmoothDrivers.h"                     // for g_rolling_fc3d_driver
#include "NumericsFwd.h"                          // for GlobalRollingFrictionContactPro...
#include "NumericsMatrix.h"                       // for NumericsMatrix
#include "SolverOptions.h"                        // for SolverOptions, SICONOS_DPA...
#include "NumericsVerbose.h"                      // for numerics_set_verbose
#include "frictionContact_test_utils.h"           // for globalRollingFrictionContact_te...
#include "test_utils.h"                           // for TestCase
#include <time.h>                                 // for clock
#include "JordanAlgebra.h"                        // for dnrm2l


int globalRollingFrictionContact_test_function(TestCase* current)
{
  int k;
  GlobalRollingFrictionContactProblem* problem = globalRollingFrictionContact_new_from_filename(current->filename);

  problem->name  = malloc(1 + strlen(current->filename));
  strcpy(problem->name, current->filename);

  char *problem_name = NULL;
  char *str = (char *) malloc(200);
  strcpy( str, problem->name );
  const char * separators = "/";
  problem_name = strtok( str, separators );
  for(int i=0; i<5; i++)
  {
    if(problem_name != NULL) problem_name = strtok ( NULL, separators );
  }

  problem_name = strtok ( problem_name, "." );



  numerics_set_verbose(1);

  FILE * foutput  =  fopen("checkinput.dat", "w");
  globalRollingFrictionContact_printInFile(problem, foutput);

  int NC = problem->numberOfContacts;
  int dim = problem->dimension;
  int n = problem->M->size1;

  int info;
  double *reaction = (double*)calloc(dim * NC, sizeof(double));
  double *velocity = (double*)calloc(dim * NC, sizeof(double));
  double *globalvelocity = calloc(n, sizeof(double));

  long clk_tck = CLOCKS_PER_SEC;



  clock_t t1 = clock();

  if(dim == 2)
  {
    info = 1;
  }
  else if(dim == 5)
  {
    info = g_rolling_fc3d_driver(problem,
                               reaction, velocity, globalvelocity,
                               current->options);
  }


  clock_t t2 = clock();

  int print_size = 10;

  printf("Norm velocity:  %Le\n", dnrm2l(NC*dim, velocity));
  printf("Norm reaction:  %Le\n", dnrm2l(NC*dim, reaction));
  printf("Norm GlobalVe:  %Le\n", dnrm2l(n, globalvelocity));

  if(dim * NC >= print_size)
  {
    printf("First values (%i)\n", print_size);
    for(k = 0 ; k < print_size; k++)
    {
      printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k, reaction[k]);
    }
    printf(" ..... \n");
    for(k = 0 ; k < print_size; k++)
    {
      printf("GlocalVelocity[%i] = %12.8e\n", k, globalvelocity[k]);
    }
  }
  else
  {
    for(k = 0 ; k < dim * NC; k++)
    {
      printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k, reaction[k]);
    }
    printf("\n");
    for(k = 0 ; k < dim*NC; k++)
    {
      printf("GlocalVelocity[%i] = %12.8e\n", k, globalvelocity[k]);
    }
  }
  printf("\n");

  for(k = 0; k < dim * NC; ++k)
  {
    info = info == 0 ? !(isfinite(velocity[k]) && isfinite(reaction[k])) : info;
  }

  for(k = 0; k < n; ++k)
  {
    info = info == 0 ? !(isfinite(globalvelocity[k])) : info;
  }

  if(!info)
    printf("test: success\n");
  else
    printf("test: failure\n");


  if (current->options->solverId == SICONOS_GLOBAL_ROLLING_FRICTION_3D_IPM)
  {
    // Take projerr value from test
    double *projerr_ptr = current->options->solverData;
    printf("\nsumry: %d  %.2e  %.2e   %5i %5i    %.6f   %s\n",
          info, current->options->dparam[SICONOS_DPARAM_RESIDU], *projerr_ptr,
          current->options->iparam[SICONOS_IPARAM_ITER_DONE], NC,
          (double)(t2-t1)/(double)clk_tck, problem_name);
  }
  else
  {
    // printf("\nsumry: %d %03.2e %04i %05.4f ", info, current->options->dparam[SICONOS_DPARAM_RESIDU], current->options->iparam[SICONOS_IPARAM_ITER_DONE], (double)(t2-t1)/(double)clk_tck);
    // printf("%02i %05i %05i %s\n\n", dim, NC, n, current->filename);
    printf("\nsumry: %d  %9.2e  %5i  %10.4f", info, current->options->dparam[SICONOS_DPARAM_RESIDU], current->options->iparam[SICONOS_IPARAM_ITER_DONE], (double)(t2-t1)/(double)clk_tck);
    printf("%3i %5i %5i     %s\n\n", dim, NC, n, current->filename);
  }

  free(reaction);
  free(velocity);
  free(globalvelocity);
  fclose(foutput);

  globalRollingFrictionContactProblem_free(problem);

  return info;
}


