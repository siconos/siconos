/* Siconos-Numerics, Copyright INRIA 2005-2011.
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
#include "ls_test_function.h"


void _LSfillParamWithRespectToSolver(SolverOptions *options, int solverId, LinearSystemProblem* problem)
{
  options->solverId = solverId;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 5;
  options->dSize = 5;
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  /*use dgels ?*/
  options->iparam[4] = 0;
  if (problem)
  {
    options->dWork = (double*) malloc(LinearSystem_getNbDwork(problem, options) * sizeof(double));
    options->iWork = (int*) malloc(LinearSystem_getNbIwork(problem, options) * sizeof(int));
  }
  else
  {
    options->dWork = NULL;
    options->iWork = NULL;
  }
  options->dparam[0] = 1e-12;
}


void _LSfillParamWithRespectToSolver_SBM(SolverOptions *options, int solverId, LinearSystemProblem* problem)
{


}



int ls_test_function(FILE * f, int solverId)
{

  int i, info = 0 ;
  LinearSystemProblem* problem = (LinearSystemProblem *)malloc(sizeof(LinearSystemProblem));

  info = LinearSystem_newFromFile(problem, f);


  NumericsOptions global_options;
  global_options.verboseMode = 1;
  SolverOptions * options ;
  options = malloc(sizeof(*options));

  options->solverId = solverId;
  printf("solverName ==> %s\n", idToName(solverId));

  _LSfillParamWithRespectToSolver(options, solverId, problem);

  options->filterOn = 1;
  double * z = malloc(problem->size * sizeof(double));
  double * w = malloc(problem->size * sizeof(double));
  for (i = 0; i < problem->size; i++)
  {
    z[i] = 0.0;
    w[i] = 0.0;
  }



  info = LinearSystem_driver(problem, z , w, options);

  for (i = 0 ; i < problem->size ; i++)
  {
    printf("z[%i] = %12.8e\t,w[%i] = %12.8e\n", i, z[i], i, w[i]);
  }

  if (!info)
  {
    printf("test succeeded err = %e \n", options->dparam[1]);
  }
  else
  {
    printf("test unsucceeded err =%e  \n", options->dparam[1]);
  }
  free(z);
  free(w);



  free(options->dWork);
  free(options->iWork);

  free(options->iparam);
  free(options->dparam);


  free(options);

  LinearSystem_freeProblem(problem);

  printf("End of test.\n");

  return info;


}

int ls_test_function_SBM(FILE * f, int solverId)
{

  int info = -1;
  return info;


}


