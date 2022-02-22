/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#include <math.h>                          // for isfinite
#include <stdio.h>                         // for printf, fclose, fopen, FILE
#include <stdlib.h>                        // for calloc, free
#include "GlobalFrictionContactProblem.h"  // for GlobalFrictionContactProblem
#include "NonSmoothDrivers.h"              // for gfc3d_driver
#include "NumericsFwd.h"                   // for GlobalFrictionContactProblem
#include "NumericsMatrix.h"                // for NumericsMatrix
#include "SolverOptions.h"                 // for SolverOptions, SICONOS_DPA...
#include "frictionContact_test_utils.h"    // for globalFrictionContact_test...
#include "test_utils.h"                    // for TestCase

#include "SiconosConfig.h"                 // for HAVE_GAMS_C_API // IWYU pragma: keep

int globalFrictionContact_test_function(TestCase* current)
{

  int k, info = -1 ;
  GlobalFrictionContactProblem* problem = globalFrictionContact_new_from_filename(current->filename);
  /* globalFrictionContact_display(problem); */


  FILE * foutput  =  fopen("checkinput.dat", "w");
  info = globalFrictionContact_printInFile(problem, foutput);

  int NC = problem->numberOfContacts;
  int dim = problem->dimension;
  int n = problem->M->size1;


  double *reaction = calloc(dim * NC, sizeof(double));
  double *velocity = calloc(dim * NC, sizeof(double));
  double *globalvelocity = calloc(n, sizeof(double));

  // --- Extra setup for options when the solver belongs to GAMS family ---
#ifdef HAVE_GAMS_C_API
  // Do we really need this?
  frictionContact_test_gams_opts(current->options);
#endif


  if(dim == 2)
  {

    info = gfc2d_driver(problem,
                        reaction, velocity, globalvelocity,
                        current->options);
  }
  else if(dim == 3)
  {
    info = gfc3d_driver(problem,
                        reaction, velocity, globalvelocity,
                        current->options);
  }

  int print_size = 10;

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
    printf("test successful, residual = %e\t, number of iterations = %i \n", current->options->dparam[SICONOS_DPARAM_RESIDU], current->options->iparam[SICONOS_IPARAM_ITER_DONE]);
  else
    printf("test unsuccessful, residual = %e, info = %d, nb iter = %d\n", current->options->dparam[SICONOS_DPARAM_RESIDU], info, current->options->iparam[SICONOS_IPARAM_ITER_DONE]);

  free(reaction);
  free(velocity);
  free(globalvelocity);
  fclose(foutput);
  globalFrictionContact_free(problem);
  return info;

}
