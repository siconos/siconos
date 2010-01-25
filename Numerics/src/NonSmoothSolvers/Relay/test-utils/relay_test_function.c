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
#include "relay_test_function.h"


int relay_test_function(FILE * f, char * solvername)
{

  int i, info = 0 ;
  Relay_Problem* problem = (Relay_Problem *)malloc(sizeof(Relay_Problem));

  info = relay_newFromFile(problem, f);

  FILE * foutput  =  fopen("./relay.verif", "w");
  info = relay_printInFile(problem, foutput);


  Numerics_Options global_options;
  global_options.verboseMode = 1;



  Solver_Options * options = malloc(sizeof(Solver_Options));

  relay_setDefaultSolverOptions(problem, options, solvername);

  int maxIter = 50000;
  double tolerance = 1e-8;
  options->iparam[0] = maxIter;
  options->dparam[0] = tolerance;


  double * z = malloc(problem->size * sizeof(double));
  double * w = malloc(problem->size * sizeof(double));

  for (i = 0; i <  problem->size ; i++)
  {
    z[i] = 0.0;
    w[i] = 0.0;
  }


  info = relay_driver(problem, z , w, options, &global_options);

  for (i = 0 ; i < problem->size ; i++)
  {
    printf("z[%i] = %12.8e\t,w[%i] = %12.8e\n", i, z[i], i, w[i]);
  }

  if (!info)
  {
    printf("test succeeded\n");
  }
  else
  {
    printf("test unsucceeded\n");
  }
  free(z);
  free(w);

  deleteSolverOptions(options);

  free(options);

  freeRelay_problem(problem);


  return info;


}

