/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2007.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */
#include "Solver_Options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void readSolverOptions(int driverType, Solver_Options* options)
{
  /* To each problem, corresponds a XXX_parameters.opt file where default parameters can be read, XXX being the problem name (LCP, FrictionContact3D ...) */

  if (verbose > 0)
    printf("\n ========== Numerics Non Smooth Solver - Read default parameters for the solver.\n ==========");

  // Checks if NUMERICSSPATH is set.
  if (getenv("SICONOSPATH") == NULL)
  {
    fprintf(stderr, "Numerics, readSolverOptions error, SICONOSPATH environment variable not set. Can not find default solver options file.\n");
    exit(EXIT_FAILURE);
  }

  FILE * ficin;
  /* Name of the default parameters file */
  char name[64];

  strcpy(name, getenv("SICONOSPATH"));
  strcat(name, "/include/Siconos/Numerics/");

  char buffer[64];
  /* Return value for reading */
  int nval;

  // set default size to 4 ...
  if (options->iparam == NULL)
    options->iparam = (int*)malloc(4 * sizeof(*options->iparam));
  if (options->dparam == NULL)
    options->dparam = (double*)malloc(4 * sizeof(*options->dparam));

  switch (driverType)
  {

  case 0:
    strcat(name, "LCP_parameters.opt");
  case 1:
    strcat(name, "dfc2D_parameters.opt");
  case 2:
    strcat(name, "pfc2D_parameters.opt");
  case 3:
    strcat(name, "FrictionContact3D_parameters.opt");
    ficin = fopen(name, "r");
    if (verbose > 0)
      printf("The default-parameters file is: %s\n", name);
    if (!ficin)
    {
      printf("Numerics, readSolverOptions error. Can not open file %60s", name);
      exit(-1);
    }
    //nval = fscanf(ficin, "%c", &(options->solverName));
    fgets(buffer, 64, ficin);
    fgets(buffer, 64, ficin);
    fgets(buffer, 64, ficin);
    /* Solver name */
    fgets(options->solverName, 64, ficin);
    fgets(buffer, 64, ficin);
    /* iparam */
    nval = fscanf(ficin, "%d%d", &(options->iparam[0]), &(options->iparam[1]));
    if (nval != 4)
    {
      fprintf(stderr, "Numerics, readSolverOptions error, wrong number of parameters for iparam.\n");
      exit(EXIT_FAILURE);

    }
    /* dparam */
    nval = fscanf(ficin, "%lf%lf", &(options->dparam[0]), &(options->dparam[1]));
    if (nval != 2)
    {
      fprintf(stderr, "Numerics, readSolverOptions error, wrong number of parameters for dparam.\n");
      exit(EXIT_FAILURE);

    }
    fclose(ficin);
    break;
  default:
    fprintf(stderr, "Numerics, readSolverOptions error, unknown problem type.\n");
    exit(EXIT_FAILURE);
  }
}

void printSolverOptions(Solver_Options* options)
{
  printf("\n ========== Numerics Non Smooth Solver parameters: \n");
  if (options->isSet == 0)
    printf("The solver parameters have not been set.\n");
  else
  {
    printf("The solver is named %s \n", options->solverName);
    if (options->iparam != NULL)
    {
      printf(" - int parameters (see %s documentation to know what is each parameter): ", options->solverName);
      for (int i = 0; i < options->iSize; ++i)
        printf("%d\t", options->iparam[i]);
      printf("\n");
    }
    if (options->dparam != NULL)
    {
      printf(" - double parameters (see %s documentation to know what is each parameter): ", options->solverName);
      for (int i = 0; i < options->dSize; ++i)
        printf("%.10lf\t", options->dparam[i]);
      printf("\n");
    }
  }
  printf("\n ================================================== \n");
}
