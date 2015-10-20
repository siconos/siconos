/* Siconos-Numerics, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#include <stdio.h>
#include <stdlib.h>
#include "NonSmoothDrivers.h"
#include "frictionContact_test_function.h"

#if defined(WITH_FCLIB)
#include <fclib.h>
#include <fclib_interface.h>
#endif


int frictionContact_test_function(FILE * f, SolverOptions * options)
{

  int k, info = -1 ;
  FrictionContactProblem* problem = (FrictionContactProblem *)malloc(sizeof(FrictionContactProblem));

  info = frictionContact_newFromFile(problem, f);

  FILE * foutput  =  fopen("checkinput.dat", "w");
  info = frictionContact_printInFile(problem, foutput);

  NumericsOptions global_options;
  setDefaultNumericsOptions(&global_options);
  global_options.verboseMode = 1; // turn verbose mode to off by default

  int NC = problem->numberOfContacts;
  int dim = problem->dimension;
  //int dim = problem->numberOfContacts;
  
  double *reaction = (double*)calloc(dim * NC, sizeof(double));
  double *velocity = (double*)calloc(dim * NC, sizeof(double));

  if (dim == 2)
  {
    info = fc2d_driver(problem,
                                    reaction , velocity,
                                    options, &global_options);
  }
  else if (dim == 3)
  {
    info = fc3d_driver(problem,
                                    reaction , velocity,
                                    options, &global_options);
  }
  else
  {
    info = 1;
  }
  printf("\n");

  int print_size =10;

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

  freeFrictionContactProblem(problem);
  fclose(foutput);

  return info;

}


#if defined(WITH_FCLIB)
int frictionContact_test_function_hdf5(const char * path, SolverOptions * options)
{

  int k, info = -1 ;
  /* FrictionContactProblem* problem = (FrictionContactProblem *)malloc(sizeof(FrictionContactProblem)); */
  /* info = frictionContact_newFromFile(problem, f); */
  
  FrictionContactProblem* problem = frictionContact_fclib_read(path);
  FILE * foutput  =  fopen("checkinput.dat", "w");
  info = frictionContact_printInFile(problem, foutput);

  
  NumericsOptions global_options;
  setDefaultNumericsOptions(&global_options);
  global_options.verboseMode = 1; // turn verbose mode to off by default

  int NC = problem->numberOfContacts;
  int dim = problem->dimension;
  //int dim = problem->numberOfContacts;
  
  double *reaction = (double*)calloc(dim * NC, sizeof(double));
  double *velocity = (double*)calloc(dim * NC, sizeof(double));

  if (dim == 2)
  {
    info = fc2d_driver(problem,
                                    reaction , velocity,
                                    options, &global_options);
  }
  else if (dim == 3)
  {
    info = fc3d_driver(problem,
                                    reaction , velocity,
                                    options, &global_options);
  }
  else
  {
    info = 1;
  }
  printf("\n");

  int print_size =10;

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

  freeFrictionContactProblem(problem);
  fclose(foutput);

  return info;

}


#endif
