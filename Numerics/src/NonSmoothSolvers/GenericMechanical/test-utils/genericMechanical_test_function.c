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
#include "genericMechanical_test_function.h"

int genericMechanical_test_function(FILE * f, SolverOptions * options)
{

  int k, info = -1 ;
  GenericMechanicalProblem* problem = genericMechnical_newFromFile(f);

  //NumericsOptions global_options;
  //global_options.verboseMode = 2; // turn verbose mode to off by default



  double *reaction = (double*)calloc(problem->size, sizeof(double));
  double *velocity = (double*)calloc(problem->size, sizeof(double));
  //setNumericsOptions(&global_options);
  info = genericMechanical_driver(problem,
                                  reaction , velocity,
                                  options);
  double err = 0;
  GenericMechanical_compute_error(problem, reaction , velocity, options->dparam[0], options, &err);
  printf("\n");
  for (k = 0 ; k < problem->size; k++)
  {
    printf("Velocity[%i] = %12.8e \t \t Reaction[%i] = %12.8e\n", k, velocity[k], k , reaction[k]);
  }
  printf("\n");

  if (!info)
  {
    printf("test succeeded info=%i err=%e and tol=%e\n", info, err, options->dparam[0]);
    if (err > options->dparam[0])
    {
      printf("but unsucceeded because err>tol\n");
      return 1;
    }
  }
  else
  {
    printf("test unsucceeded\n");
  }
  free(reaction);
  free(velocity);

  freeGenericMechanicalProblem(problem, NUMERICS_GMP_FREE_MATRIX | NUMERICS_GMP_FREE_GMP);

  return info;

}


