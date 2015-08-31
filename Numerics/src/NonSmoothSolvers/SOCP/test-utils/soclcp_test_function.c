/* Siconos-Numerics, Copyright INRIA 2005-2012.
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

#include "soclcp_test_function.h"



int soclcp_test_function(FILE * f, SolverOptions * options)
{

  int k, info = -1 ;
  SecondOrderConeLinearComplementarityProblem* problem = (SecondOrderConeLinearComplementarityProblem *)malloc(sizeof(SecondOrderConeLinearComplementarityProblem));

  info = secondOrderConeLinearComplementarityProblem_newFromFile(problem, f);

  FILE * foutput  =  fopen("checkinput.dat", "w");

  info = secondOrderConeLinearComplementarityProblem_printInFile(problem, foutput);

  /* secondOrderConeLinearComplementarityProblem_display(problem); */

  NumericsOptions global_options;
  setDefaultNumericsOptions(&global_options);
  global_options.verboseMode = 1; // turn verbose mode to off by default


  int n = problem->n;

  double *r = (double*)malloc(n * sizeof(double));
  double *v = (double*)malloc(n * sizeof(double));
  for(k = 0 ; k <n; k++)
  {
    v[k] = 0.0;
    r[k] = 0.0;
  }
  info = soclcp_driver(problem,
                       r , v,
                       options, &global_options);

  printf("\n");
  for(k = 0 ; k < n; k++)
  {
    printf("v[%i] = %12.8e \t \t r[%i] = %12.8e\n", k, v[k], k , r[k]);
  }
  printf("\n");

  if(!info)
  {
    printf("test succeeded\n");
  }
  else
  {
    printf("test unsuccessful\n");
  }
  free(r);
  free(v);

  freeSecondOrderConeLinearComplementarityProblem(problem);
  fclose(foutput);

  return info;

}
