
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

#include <math.h>                                         // for fabs
#include <stdio.h>                                        // for printf
#include <stdlib.h>                                       // for free, malloc
#include <string.h>                                       // for NULL, memcpy
#include "FrictionContactProblem.h"                       // for FrictionCon...
#include "NumericsFwd.h"                                  // for SecondOrder...
#include "SOCLCP_Solvers.h"                               // for soclcp_nsgs
#include "SOCLCP_cst.h"                                   // for SICONOS_SOC...
#include "SecondOrderConeLinearComplementarityProblem.h"  // for SecondOrder...
#include "SiconosBlas.h"                                  // for cblas_dnrm2
#include "SolverOptions.h"                                // for SolverOptions
#include "fc3d_Solvers.h"                                 // for fc3d_SOCLCP
#include "fc3d_compute_error.h"                           // for fc3d_comput...
#include "numerics_verbose.h"                             // for verbose


/** pointer to function used to call internal solver for proximal point solver */
typedef void (*soclcp_InternalSolverPtr)(SecondOrderConeLinearComplementarityProblem*, double*, double*, int *, SolverOptions *);

void fc3d_SOCLCP(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options)
{
  /* int and double parameters */
  double* dparam = options->dparam;

  /* Number of contacts */
  int nc = problem->numberOfContacts;

  /* Tolerance */
  double tolerance = dparam[SICONOS_DPARAM_TOL];
  double norm_q = cblas_dnrm2(nc*3, problem->q, 1);
  /*****  Fixed Point Iterations *****/
  double error = 1.; /* Current error */

  soclcp_InternalSolverPtr internalsolver;

  SecondOrderConeLinearComplementarityProblem* soclcp = (SecondOrderConeLinearComplementarityProblem *)malloc(sizeof(SecondOrderConeLinearComplementarityProblem));
  soclcp->n = problem->numberOfContacts * problem->dimension;
  soclcp->nc= problem->numberOfContacts;
  soclcp->M = problem->M;
  soclcp->q = (double *) malloc(soclcp->n * sizeof(double));
  soclcp->tau = problem->mu;
  soclcp->coneIndex = (unsigned int *) malloc((soclcp->nc+1) * sizeof(unsigned int));

  memcpy(soclcp->q, problem->q, (soclcp->n) * sizeof(double));

  for(int i=0; i <= soclcp->nc; ++i)
  {
    soclcp->coneIndex[i] = 3*i;
  }

  options->solverId = SICONOS_SOCLCP_NSGS;
  // if (options->solverId == SICONOS_SOCLCP_NSGS)
  // This is the only allowed option, at the time
  {
    if(verbose == 1)
      printf(" ========================== Call NSGS solver SOCLCP problem ==========================\n");
    internalsolver = &soclcp_nsgs;
    //internalsolver_options->internalSolvers->dWork = options->dWork;
  }
  /* else */
  /* { */
  /*   fprintf(stderr, "Numerics, fc3d_SOCLCP failed. Unknown internal solver.\n"); */
  /*   exit(EXIT_FAILURE); */
  /* } */
  (*internalsolver)(soclcp, reaction, velocity, info, options);

  error = options->dparam[SICONOS_DPARAM_RESIDU];
  double real_error=0.0;

  fc3d_compute_error(problem, reaction, velocity, tolerance, options, norm_q, &real_error);

  if(options->callback)
  {
    options->callback->collectStatsIteration(options->callback->env, nc * 3,
        reaction, velocity, error, NULL);
  }

  if(verbose > 0)
  {
    printf("--------------- FC3D - SOCLCP - # Iteration %i Final Residual = %14.7e\n", options->iparam[SICONOS_IPARAM_ITER_DONE], error);
    printf("--------------- FC3D - SOCLCP - #              error of the real problem = %14.7e\n", real_error);
    printf("--------------- FC3D - SOCLCP - #              gap with the real problem = %14.7e\n", fabs(real_error-error));
  }

  free(soclcp->q);
  free(soclcp->coneIndex);
  free(soclcp);
  dparam[SICONOS_DPARAM_RESIDU] = error;
  dparam[2] = fabs(real_error-error);

}
