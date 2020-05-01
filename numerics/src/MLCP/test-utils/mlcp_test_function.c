/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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
#define _XOPEN_SOURCE 700
#include <math.h>                          // for isfinite
#include <stdio.h>                         // for printf
#include <stdlib.h>                        // for calloc, free, malloc
#include "MixedLinearComplementarityProblem.h"  // for MixedLinearComplementarityProblem
#include "NonSmoothDrivers.h"              // for mixedlinearComplementarity_driver
#include "MLCP_Solvers.h"
//#include "NumericsFwd.h"                   // for LinearComplementarityProblem
#include "SolverOptions.h"                 // for SICONOS_DPARAM_RESIDU, Sol...
#include "mlcp_test_utils.h"                // for lcp_test_function
#include "test_utils.h"                    // for TestCase



int mlcp_test_function(TestCase * current)
{
  //numerics_set_verbose(2);
  int i, info = 0 ;
  MixedLinearComplementarityProblem* problem = (MixedLinearComplementarityProblem *)malloc(sizeof(MixedLinearComplementarityProblem));
  info = mixedLinearComplementarity_newFromFilename(problem, current->filename);


  /* mixedLinearComplementarity_display(problem); */

  double * z = (double *)calloc(problem->n+problem->m, sizeof(double));
  double * w = (double *)calloc(problem->n+problem->m, sizeof(double));

  mlcp_driver_init(problem, current->options);
  info = mlcp_driver(problem, z, w, current->options);
  mlcp_driver_reset(problem, current->options);

  int size_max = problem->n+problem->m;
  bool brief= false;
  if (size_max > 10)
  {
    size_max=10;
    brief= true;
  }


  for(i = 0 ; i < size_max; i++)
  {
    printf("z[%i] = %12.8e\t,w[%i] = %12.8e\n", i, z[i], i, w[i]);
  }
  if (brief)
  {
    printf(".....\n");
    printf("z[%i] = %12.8e\t,w[%i] = %12.8e\n", problem->n+problem->m-1, z[problem->n+problem->m-1], problem->n+problem->m-1, w[problem->n+problem->m-1]);
  }
  for(i = 0 ; i < problem->n+problem->m ; i++)
  {
    info = info == 0 ? !(isfinite(z[i]) && isfinite(w[i])): info;
  }

  if(!info)
    printf("test succeeded err = %e \n", current->options->dparam[SICONOS_DPARAM_RESIDU]);
  else
    printf("test unsuccessful err = %e  \n", current->options->dparam[SICONOS_DPARAM_RESIDU]);
  free(z);
  free(w);
  mixedLinearComplementarity_free(problem);
  printf("End of test.\n");
  return info;
}
