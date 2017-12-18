/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "NonSmoothDrivers.h"
#include "SecondOrderConeLinearComplementarityProblem.h"
#include "SOCLCP_Solvers.h"
#include "soclcp_test_function.h"
#include "SiconosCompat.h"

#ifdef __cplusplus
using namespace std;
#endif

int soclcp_test_function(FILE * f, SolverOptions * options)
{

  int k, info = -1 ;
  SecondOrderConeLinearComplementarityProblem* problem = (SecondOrderConeLinearComplementarityProblem *)malloc(sizeof(SecondOrderConeLinearComplementarityProblem));

  assert(f);
  assert(problem);
  
  info = secondOrderConeLinearComplementarityProblem_newFromFile(problem, f);

  FILE * foutput  =  fopen("checkinput.dat", "w");

  info = secondOrderConeLinearComplementarityProblem_printInFile(problem, foutput);

  /* secondOrderConeLinearComplementarityProblem_display(problem); */
  int n = problem->n;

  double *r = (double*)calloc(n, sizeof(double));
  double *v = (double*)calloc(n, sizeof(double));

  info = soclcp_driver(problem, r , v, options);

  printf("\n");
  for(k = 0 ; k < n; k++)
  {
    printf("v[%i] = %12.8e \t \t r[%i] = %12.8e\n", k, v[k], k , r[k]);
    info = info == 0 ? !(isfinite(v[k]) && isfinite(r[k])) : info;
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
