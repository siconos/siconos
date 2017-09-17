
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

#include "fc3d_Solvers.h"
#include "fc3d_compute_error.h"
#include "SOCLCP_Solvers.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

//#define VERBOSE_DEBUG
#include "Friction_cst.h"
#include "SiconosBlas.h"

#define DEBUG_MESSAGES
#define DEBUG_STDOUT
#include "debug.h"
#include "numerics_verbose.h"


/** pointer to function used to call internal solver for proximal point solver */
typedef void (*soclcp_InternalSolverPtr)(SecondOrderConeLinearComplementarityProblem*, double*, double*, int *, SolverOptions *);

void fc3d_SOCLCP(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options)
{
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;

  /* Number of contacts */
  int nc = problem->numberOfContacts;

  /* Tolerance */
  double tolerance = dparam[0];
  double norm_q = cblas_dnrm2(nc*3 , problem->q , 1);
 


  if (options->numberOfInternalSolvers < 1)
  {
    numerics_error("fc3d_SOCLCP", "The SOCLCP for FrictionContactProblem needs options for the internal solvers, options[0].numberOfInternalSolvers should be >1");
  }

  SolverOptions * internalsolver_options = options->internalSolvers;

  if (verbose > 0)
  {
    printf("Local solver data :");
    solver_options_print(internalsolver_options);
  }


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

  for (int i=0; i <= soclcp->nc; ++i)
  {
    soclcp->coneIndex[i] = 3*i;
  }

  if (internalsolver_options->solverId == SICONOS_SOCLCP_NSGS)
  {
    if (verbose == 1)
      printf(" ========================== Call NSGS solver SOCLCP problem ==========================\n");
    internalsolver = &soclcp_nsgs;
    //internalsolver_options->internalSolvers->dWork = options->dWork;
  }
  else
  {
    fprintf(stderr, "Numerics, fc3d_SOCLCP failed. Unknown internal solver.\n");
    exit(EXIT_FAILURE);
  }
  internalsolver_options->dparam[0]=options->dparam[0];
  internalsolver_options->iparam[0]=options->iparam[0];

  (*internalsolver)(soclcp, reaction , velocity , info , internalsolver_options);

  error = internalsolver_options->dparam[1];


  double real_error=0.0;

  fc3d_compute_error(problem, reaction , velocity, tolerance, options, norm_q, &real_error);

  if (options->callback)
  {
    options->callback->collectStatsIteration(options->callback->env, nc * 3,
                                             reaction, velocity, error, NULL);
  }

  if (verbose > 0)
  {
    printf("----------------------------------- FC3D - SOCLCP - # Iteration %i Final Residual = %14.7e\n", internalsolver_options->iparam[7], error);
    printf("----------------------------------- FC3D - SOCLCP - #              error of the real problem = %14.7e\n", real_error );
    printf("----------------------------------- FC3D - SOCLCP - #              gap with the real problem = %14.7e\n", fabs(real_error-error) );
  }

  free(soclcp->q);
  free(soclcp->coneIndex);
  free(soclcp);


  if (internalsolver_options->internalSolvers != NULL)
    internalsolver_options->internalSolvers->dWork = NULL;

  dparam[0] = tolerance;
  dparam[1] = error;
  dparam[2] = fabs(real_error-error);
  iparam[7] = internalsolver_options->iparam[7];

}



int fc3d_SOCLCP_setDefaultSolverOptions(SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the SOCLCP Solver\n");
  }

  options->solverId = SICONOS_FRICTION_3D_SOCLCP;
  options->numberOfInternalSolvers = 1;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 8;
  options->dSize = 8;
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  options->dWork = NULL;
  solver_options_nullify(options);
  options->iparam[0] = 1000;
  options->dparam[0] = 1e-4;
  options->internalSolvers = (SolverOptions *)malloc(sizeof(SolverOptions));

  soclcp_nsgs_setDefaultSolverOptions(options->internalSolvers);
  options->internalSolvers->iparam[0] =10000;
  return 0;
}
