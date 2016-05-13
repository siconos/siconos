/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.

 * Copyright 2016 INRIA.

 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at

 * http://www.apache.org/licenses/LICENSE-2.0

 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "Relay_Solvers.h"
#include "LCP_Solvers.h"
#include <assert.h>
#include "relay_cst.h"

void relay_enum(RelayProblem* problem, double *z, double *w, int *info, SolverOptions* options, NumericsOptions* global_options)
{
  int i;
  // conversion into LCP
  LinearComplementarityProblem* lcp_problem = (LinearComplementarityProblem*)malloc(sizeof(LinearComplementarityProblem));


  /* Relay_display(problem); */

  relay_to_lcp(problem, lcp_problem);

  /* linearComplementarity_display(lcp_problem);  */

  double *zlcp = (double*)malloc(lcp_problem->size * sizeof(double));
  double *wlcp = (double*)malloc(lcp_problem->size * sizeof(double));

  /*  FILE * fcheck = fopen("lcp_relay.dat","w"); */
  /*  info = linearComplementarity_printInFile(lcp_problem,fcheck); */


  // Call the lcp_solver
  options->solverId = SICONOS_LCP_ENUM;
  lcp_enum_init(lcp_problem, options, 1);

  * info = linearComplementarity_driver(lcp_problem, zlcp , wlcp, options, global_options);
  if (options->filterOn > 0)
    lcp_compute_error(lcp_problem, zlcp, wlcp, options->dparam[0], &(options->dparam[1]));

  lcp_enum_reset(lcp_problem, options, 1);
  options->solverId = SICONOS_RELAY_ENUM;

  // Conversion of result
  for (i = 0; i < problem->size; i++)
  {
    /* z[i] = 1.0/2.0*(zlcp[i]-wlcp[i+problem->size]); works only for ub=1 and lb=-1 */
    z[i] = zlcp[i] +  problem->lb[i];

    w[i] = wlcp[i] - zlcp[i + problem->size];
    //printf("w[ %i]=%12.10e\n", i, w[i]);
  }

  /* for (i=0; i< lcp_problem->size; i++){ */
  /*   printf("zlcp[ %i]=%12.10e,\t wlcp[ %i]=%12.10e \n", i, zlcp[i],i, wlcp[i]); */
  /* } */
  /* printf("\n"); */

  /* for (i=0; i< problem->size; i++){ */
  /*   printf("z[ %i]=%12.10e,\t w[ %i]=%12.10e\n", i, z[i],i, w[i]); */
  /* } */


  /*        printf("\n"); */
  free(zlcp);
  free(wlcp);
  freeLinearComplementarityProblem(lcp_problem);

}

int relay_enum_setDefaultSolverOptions(RelayProblem* problem, SolverOptions* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the ENUM Solver\n");
  }

  /*  strcpy(options->solverName,"ENUM");*/
  options->solverId = SICONOS_RELAY_ENUM;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 5;
  options->dSize = 5;
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  options->dWork = NULL ;/* (double*) malloc((3*problem->size +problem->size*problem->size)*sizeof(double)); */
  options->iWork = NULL ; /* (int*) malloc(2*problem->size*sizeof(int)); */
  for (i = 0; i < 5; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->dparam[0] = 1e-12;




  return 0;
}

