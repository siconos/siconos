/* Siconos-Numerics, Copyright INRIA 2005-2010.
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
#include "frictionContact3D_test_function.h"


void frictionContact3D_fillParamWithRespectToSolver(Solver_Options *options, char * solvername, FrictionContact_Problem* problem)
{


  int maxIter = 1001;
  double tolerance = 1e-8;
  double localtolerance = 1e-8;
  int localsolver = 1;

  if (strcmp(solvername , "NSGS") == 0)
  {
    options->iSize = 5;
    options->dSize = 5;
    options->iparam[0] = maxIter;
    options->iparam[4] = localsolver;
    options->dparam[0] = tolerance;
    options->dparam[2] = localtolerance;
  }

}

int frictionContact3D_test_function(FILE * f, char * solvername, int * iparam, double * dparam)
{

  int i, k, info = -1 ;
  FrictionContact_Problem* problem = (FrictionContact_Problem *)malloc(sizeof(FrictionContact_Problem));

  info = frictionContact_newFromFile(problem, f);

  FILE * foutput  =  fopen("checkinput.dat", "w");
  info = frictionContact_printInFile(problem, foutput);


  Numerics_Options global_options;
  global_options.verboseMode = 1; // turn verbose mode to off by default

  Solver_Options * options ;
  options = malloc(sizeof(*options));
  options->isSet = 1;
  options->filterOn = 1;

  strcpy(options->solverName, solvername);
  printf("solverName ==> %s\n", options->solverName);
  options->iSize = 5;
  options->dSize = 5;
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  options->iWork = NULL;

  if (strcmp(solvername , "NSGS") == 0)
  {
    options->iSize = 5;
    options->dSize = 5;
    for (i = 0; i < 5; i++)
    {
      options->iparam[i] = iparam[i];
      options->dparam[i] = dparam[i];
    }
  }

  int NC = problem->numberOfContacts;

  double *reaction = (double*)malloc(3 * NC * sizeof(double));
  double *velocity = (double*)malloc(3 * NC * sizeof(double));

  info = frictionContact3D_driver(problem,
                                  reaction , velocity,
                                  options, &global_options);
  printf("\n");
  for (k = 0 ; k < 3 * NC; k++)
  {
    printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k , reaction[k]);
  }
  printf("\n");

  if (!info)
  {
    printf("test succeeded\n");
  }
  else
  {
    printf("test unsucceeded\n");
  }
  free(reaction);
  free(velocity);

  free(options->iparam);
  free(options->dparam);


  if (!options->dWork) free(options->dWork);
  if (!options->iWork) free(options->iWork);

  free(options);

  freeFrictionContact_problem(problem);


  return info;


}


