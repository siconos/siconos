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


void relay_fillParamWithRespectToSolver(Solver_Options *options, char * solvername, Relay_Problem* problem)
{
  int maxIter = 50000;
  double tolerance = 1e-8;
  if (strcmp(solvername , "PGS") == 0 || strcmp(solvername , "CPG") == 0 || strcmp(solvername , "Lemke") == 0 || strcmp(solvername , "NewtonMin") == 0)
  {
    options->iSize = 2;
    options->dSize = 2;
    options->iparam[0] = maxIter;
    options->dparam[0] = tolerance;
  }

  else if (strcmp(solvername , "RPGS") == 0)
  {
    options->iSize = 2;
    options->dSize = 3;
    options->iparam[0] = maxIter;
    options->dparam[0] = tolerance;
    options->dparam[2] = 1.0;
  }
  else if (strcmp(solvername , "Latin") == 0)
  {
    options->iSize = 2;
    options->dSize = 3;
    options->iparam[0] = maxIter;
    options->dparam[0] = tolerance;
    options->dparam[2] = 0.3;
  }
  else if (strcmp(solvername , "Latin_w") == 0)
  {
    options->iSize = 2;
    options->dSize = 4;
    options->iparam[0] = maxIter;
    options->dparam[0] = tolerance;
    options->dparam[2] = 0.3;
    options->dparam[3] = 1.0;
  }
  else if (strcmp(solvername , "PATH") == 0 || strcmp(solvername , "QP") == 0 || strcmp(solvername , "NSQP") == 0)
  {
    options->iSize = 0;
    options->dSize = 2;
    options->dparam[0] = tolerance;
  }
  else if (strcmp(solvername , "ENUM") == 0)
  {
    options->iSize = 2;
    options->dSize = 2;
    options->dparam[0] = tolerance;
    options->dWork = (double*) malloc((3 * problem->size + problem->size * problem->size) * sizeof(double));
    options->iWork = (int*) malloc(2 * problem->size * sizeof(int));
  }
  else if (strcmp(solvername , "NewtonFB") == 0)
  {
    options->iSize = 2;
    options->dSize = 2;
    options->iparam[0] = maxIter;
    options->dparam[0] = tolerance;

  }


}

int relay_test_function(FILE * f, char * solvername)
{

  int i, info = 0 ;
  Relay_Problem* problem = (Relay_Problem *)malloc(sizeof(Relay_Problem));

  info = relay_newFromFile(problem, f);

  FILE * foutput  =  fopen("./lcp_mmc.verif", "w");
  info = relay_printInFile(problem, foutput);


  Numerics_Options global_options;
  global_options.verboseMode = 1;
  Solver_Options * options ;
  options = malloc(sizeof(*options));

  strcpy(options->solverName, solvername);
  printf("solvername ==> %s\n", options->solverName);
  options->dWork = NULL;
  options->iWork = NULL;

  options->iSize = 10;
  options->dSize = 10;
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  for (i = 0; i < 10; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  relay_fillParamWithRespectToSolver(options, solvername, problem);

  options->isSet = 1;
  options->filterOn = 1;
  double * z = malloc(problem->size * sizeof(double));
  double * w = malloc(problem->size * sizeof(double));


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

  free(options->iparam);
  free(options->dparam);

  if (!options->dWork) free(options->dWork);
  if (!options->iWork) free(options->iWork);

  free(options);

  freeRelay_problem(problem);


  return info;


}

