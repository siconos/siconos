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
#include "primalFrictionContact_test_function.h"


int primalFrictionContact_test_function(FILE * f, SolverOptions * options)
{

  int k, info = -1 ;
  PrimalFrictionContactProblem* problem = (PrimalFrictionContactProblem *)malloc(sizeof(PrimalFrictionContactProblem));

  info = primalFrictionContact_newFromFile(problem, f);

  FILE * foutput  =  fopen("checkinput.dat", "w");
  info = primalFrictionContact_printInFile(problem, foutput);


  Numerics_Options global_options;
  global_options.verboseMode = 1; // turn verbose mode to off by default

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
    info = primalFrictionContact3D_driver(problem,
                                          reaction , velocity, globalvelocity,
                                          options, &global_options);
  }
  printf("\n");
  for (k = 0 ; k < dim * NC; k++)
  {
    printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k , reaction[k]);
  }
  for (k = 0 ; k < n; k++)
  {
    printf("GlocalVelocity[%i] = %12.8e\n", k, globalvelocity[k]);
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
  free(globalvelocity);

  freePrimalFrictionContact_problem(problem);


  return info;

}


