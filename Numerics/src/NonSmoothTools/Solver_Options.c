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

void recursive_printSolverOptions(Solver_Options* options, int level)
{
  char* marge;
  int i;
  marge = (char*) malloc((level + 1) * sizeof(char));
  for (i = 0; i < level; i++)
    marge[i] = ' ';
  marge[level] = '\0';

  printf("%s\n ========== Numerics Non Smooth Solver parameters: \n", marge);
  if (options->isSet == 0)
    printf("%sThe solver parameters have not been set. \t options->isSet = %i \n", marge, options->isSet);
  else
  {
    printf("%sThe solver parameters below have  been set \t options->isSet = %i\n", marge, options->isSet);
    printf("%sName of the solver\t\t\t\t options->solverName = %s \n", marge, options->solverName);
    if (options->iparam != NULL)
    {
      printf("%sint parameters \t\t\t\t\t options->iparam\n", marge);
      printf("%ssize of the int parameters\t\t\t options->iSize = %i\n", marge, options->iSize);
      for (int i = 0; i < options->iSize; ++i)
        printf("%s\t\t\t\t\t\t options->iparam[%i] = %d\n", marge, i, options->iparam[i]);
    }
    if (options->dparam != NULL)
    {
      printf("%sdouble parameters \t\t\t\t options->dparam\n", marge);
      printf("%ssize of the double parameters\t\t\t options->dSize = %i\n", marge, options->dSize);
      for (int i = 0; i < options->iSize; ++i)
        printf("%s\t\t\t\t\t\t options->dparam[%i] = %.6le\n", marge, i, options->dparam[i]);
    }
  }
  if (options->iWork == NULL)
  {
    printf("%sinteger work array have not been allocated. \t options->iWork = NULL \n", marge);
  }
  else
  {
    printf("%sinteger work array have been allocated. \t options->iWork = %p \n", marge, options->iWork);
    printf("%sinteger work array size \t\t\t options->iSize = %i \n", marge, options->iSize);
  }
  if (options->dWork == NULL)
  {
    printf("%sdouble work array have not been allocated. \t options->dWork = NULL \n", marge);
  }
  else
  {
    printf("%sdouble work array have been allocated. \t options->dWork = %p \n", marge, options->dWork);
    printf("%sdouble work array size \t\t\t options->dSize = %i \n", marge, options->dSize);
  }




  printf("%sSee %s documentation for parameters definition)\n", marge, options->solverName);

  printf("\n");

  printf("%snumber of internal (or local) solvers \t\t options->numberOfInternalSolvers = %i\n", marge, options->numberOfInternalSolvers);
  for (i = 0; i < options->numberOfInternalSolvers; i++)
  {
    recursive_printSolverOptions(options->internalSolvers + i, level + 1);
  }
  free(marge);

}
void printSolverOptions(Solver_Options* options)
{
  recursive_printSolverOptions(options, 0);
}
void recursive_deleteSolverOptions(Solver_Options* op)
{

  for (int i = 0; i < op->numberOfInternalSolvers; i++)
    recursive_deleteSolverOptions(&(op->internalSolvers[i]));

  if (op->numberOfInternalSolvers && op->internalSolvers)
    free(op->internalSolvers);
  op->internalSolvers = 0;
  if (op->iparam)
    free(op->iparam);
  op->iparam = 0;
  if (op->dparam)
    free(op->dparam);
  op->dparam = 0;
  if (op->iWork)
    free(op->iWork);
  op->iWork = 0;
  if (op->dWork)
    free(op->dWork);
  op->dWork = 0;
}


void deleteSolverOptions(Solver_Options* op)
{
  for (int i = 0; i < op->numberOfInternalSolvers; i++)
    recursive_deleteSolverOptions(&(op->internalSolvers[i]));
  if (op->numberOfInternalSolvers && op->internalSolvers)
    free(op->internalSolvers);
  op->internalSolvers = 0;
  if (op->iparam)
    free(op->iparam);
  op->iparam = 0;
  if (op->dparam)
    free(op->dparam);
  op->dparam = 0;
  if (op->iWork)
    free(op->iWork);
  op->iWork = 0;
  if (op->dWork)
    free(op->dWork);
  op->dWork = 0;



}

