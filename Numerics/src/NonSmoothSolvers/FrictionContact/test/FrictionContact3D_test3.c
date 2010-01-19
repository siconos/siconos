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
#include "frictionContact_test_function.h"

int setSolver_Options(Solver_Options * options)
{

  int i;

  strcpy(options->solverName, "NSGS");
  printf("solverName ==> %s\n", options->solverName);
  options->numberOfInternalSolvers = 1;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 5;
  options->dSize = 5;
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  options->iWork = NULL;
  for (i = 0; i < 5; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->iparam[0] = 1001;
  options->dparam[0] = 1e-16;


  strcpy(options[1].solverName, "ProjectionOnConeWithLocalIteration");
  options[1].numberOfInternalSolvers = 0;
  options[1].isSet = 1;
  options[1].filterOn = 1;
  options[1].iSize = 5;
  options[1].dSize = 5;
  options[1].iparam = (int *)malloc(options->iSize * sizeof(int));
  options[1].dparam = (double *)malloc(options->dSize * sizeof(double));
  options[1].dWork = NULL;
  options[1].iWork = NULL;

  for (i = 0; i < 5; i++)
  {
    options[1].iparam[i] = 0;
    options[1].dparam[i] = 0.0;
  }

  options[1].iparam[0] = 10;
  options[1].dparam[0] = 1e-3;


  return 0;
}

void freeSolver_Options(Solver_Options * options)
{
  for (int i = 0; i < options->numberOfInternalSolvers + 1; i++)
  {
    free(options[i].iparam);
    free(options[i].dparam);
    if (!options[i].dWork) free(options[i].dWork);
    if (!options[i].iWork) free(options[i].iWork);
  }
  free(options);
}


int main(void)
{
  int info = 0 ;
  printf("Test on ./data/Example1_Fc3D_SBM.dat\n");

  FILE * finput  =  fopen("./data/Example1_Fc3D_SBM.dat", "r");
  int nbsolvers = 2;
  Solver_Options * options = (Solver_Options *)malloc(nbsolvers * sizeof(*options));  ;
  info = setSolver_Options(options);
  printf("solverName ==> %s\n", options->solverName);
  info = frictionContact_test_function(finput, options);
  freeSolver_Options(options);
  fclose(finput);
  printf("End of test on ./data/Example1_Fc3D_SBM.dat\n");
  return info;
}
