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
#include "Solver_Options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
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
    strcat(name, "FrictionContact2D_parameters.opt");
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
    nval = fscanf(ficin, "%lf%lf%lf", &(options->dparam[0]), &(options->dparam[1]), &(options->dparam[2]));
    if (nval != 3)
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
    printf("The solver parameters have not been set. \t options->isSet = %i \n", options->isSet);
  else
  {
    printf("The solver parameters below have  been set \t options->isSet = %i\n", options->isSet);
    printf("Name of the solver\t\t\t\t options->solverName = %s \n", options->solverName);
    printf("number of internal (or local) solvers \t\t options->numberOfInternalSolvers = %i\n", options->numberOfInternalSolvers);
    if (options->numberOfInternalSolvers > 0)
    {
      for (int i = 1; i < options->numberOfInternalSolvers + 1; i++)
      {
        assert(&options[i]);
        printf("Name internal or local solver\t\t\t options[%i].solverName = %s \t \n", i, options[i].solverName);
      }

    }
    if (options->iparam != NULL)
    {
      printf("int parameters \t\t\t\t\t options->iparam\n");
      printf("size of the int parameters\t\t\t options->iSize = %i\n", options->iSize);
      for (int i = 0; i < options->iSize; ++i)
        printf("\t\t\t\t\t\t options->iparam[%i] = %d\n", i, options->iparam[i]);
    }
    if (options->dparam != NULL)
    {
      printf("double parameters \t\t\t\t options->dparam\n");
      printf("size of the double parameters\t\t\t options->dSize = %i\n", options->dSize);
      for (int i = 0; i < options->iSize; ++i)
        printf("\t\t\t\t\t\t options->dparam[%i] = %.6le\n", i, options->dparam[i]);
    }
  }
  if (options->iWork == NULL)
  {
    printf("integer work array have not been allocated. \t options->iWork = NULL \n");
  }
  else
  {
    printf("integer work array have been allocated. \t options->iWork = %p \n", options->iWork);
    printf("integer work array size \t\t\t options->iSize = %i \n", options->iSize);
  }
  if (options->dWork == NULL)
  {
    printf("double work array have not been allocated. \t options->dWork = NULL \n");
  }
  else
  {
    printf("double work array have been allocated. \t options->dWork = %p \n", options->dWork);
    printf("double work array size \t\t\t options->dSize = %i \n", options->dSize);
  }




  printf("See %s documentation for parameters definition)\n", options->solverName);

  printf("\n");

}

void deleteSolverOptions(Solver_Options* op)
{
  if (op->iparam)
    free(op->iparam);
  if (op->dparam)
    free(op->dparam);
  op->iparam = 0;
  op->dparam = 0;
}
