
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
#include "VI_cst.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
//#define VERBOSE_DEBUG
#include "Friction_cst.h"
#include "sanitizer.h"

#define DEBUG_MESSAGES
#define DEBUG_STDOUT
#include "debug.h"
#include "numerics_verbose.h"



/** pointer to function used to call internal solver for proximal point solver */
typedef void (*soclcp_InternalSolverPtr)(SecondOrderConeLinearComplementarityProblem*, double*, double*, int *, SolverOptions *);

void fc3d_ACLMFixedPoint(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options)
{
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;

  /* Number of contacts */
  int nc = problem->numberOfContacts;


  /* Maximum number of iterations */
  int itermax = iparam[0];
  /* Tolerance */
  double tolerance = dparam[0];
  double norm_q = cblas_dnrm2(nc*3 , problem->q , 1);



  if (options->numberOfInternalSolvers < 1)
  {
    numerics_error("fc3d_ACLMFixedpoint", "The ACLM Fixed Point method needs options for the internal solvers, options[0].numberOfInternalSolvers should be >1");
  }

  SolverOptions * internalsolver_options = options->internalSolvers;

  if (verbose > 0)
  {
    solver_options_print(options);
  }


  /*****  Fixed Point Iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;
  soclcp_InternalSolverPtr internalsolver;

  SecondOrderConeLinearComplementarityProblem* soclcp = (SecondOrderConeLinearComplementarityProblem *)malloc(sizeof(SecondOrderConeLinearComplementarityProblem));
  soclcp->n = problem->numberOfContacts * problem->dimension;
  soclcp->nc= problem->numberOfContacts;
  soclcp->M = problem-> M;
  soclcp->q = (double *) malloc(soclcp->n * sizeof(double));
  soclcp->tau = problem->mu;
  soclcp->coneIndex = (unsigned int *) malloc((soclcp->nc+1) * sizeof(unsigned int));

  for (int i=0; i < soclcp->n; i++)
  {
    soclcp->q[i]=problem->q[i];
  }
  for (int i=0; i < soclcp->nc+1; i++)
  {
    soclcp->coneIndex[i] = 3*i;
  }

  if (internalsolver_options->solverId == SICONOS_SOCLCP_NSGS)
  {
    if (verbose == 1)
      printf(" ========================== Call NSGS solver SOCLCP problem ==========================\n");
    internalsolver = &soclcp_nsgs;
  }
  else if (internalsolver_options->solverId == SICONOS_SOCLCP_VI_FPP)
  {
    if (verbose == 1)
      printf(" ========================== Call VI_FPP solver SOCLCP problem ==========================\n");
    internalsolver = &soclcp_VI_FixedPointProjection;
  }
    else if (internalsolver_options->solverId == SICONOS_SOCLCP_VI_EG)
  {
    if (verbose == 1)
      printf(" ========================== Call VI_EG solver SOCLCP problem ==========================\n");
    internalsolver = &soclcp_VI_ExtraGradient;
  }
  else
  {
    fprintf(stderr, "Numerics, fc3d_ACLMFixedPoint failed. Unknown internal solver.\n");
    exit(EXIT_FAILURE);
  }

  double normUT;
  int cumul_iter =0;
  while ((iter < itermax) && (hasNotConverged > 0))
  {
    ++iter;
    // internal solver for the regularized problem

    /* Compute the value of the initial value of q */
    for (int ic = 0 ; ic < nc ; ic++)
    {
      normUT = sqrt(velocity[ic*3+1] * velocity[ic*3+1] + velocity[ic*3+2] * velocity[ic*3+2]);
      soclcp->q[3*ic] = problem->q[3*ic] + problem->mu[ic]*normUT;
    }

    fc3d_FixedPoint_set_internalsolver_tolerance(problem,options,internalsolver_options, error);


    (*internalsolver)(soclcp, reaction , velocity , info , internalsolver_options);
    cumul_iter +=  internalsolver_options->iparam[7];
    /* **** Criterium convergence **** */

    fc3d_compute_error(problem, reaction , velocity, tolerance, options, norm_q, &error);

    if (options->callback)
    {
      options->callback->collectStatsIteration(options->callback->env, nc * 3,
                                      reaction, velocity, error, NULL);
    }

    if (verbose > 0)
      printf("------------------------ FC3D - ACLMFP - Iteration %i Residual = %14.7e\n", iter, error);

    if (error < tolerance) hasNotConverged = 0;
    *info = hasNotConverged;
  }
  if (verbose > 0)
  {
    printf("----------------------------------- FC3D - ACLMFP - # Iteration %i Final Residual = %14.7e\n", iter, error);
    printf("----------------------------------- FC3D - ACLMFP - #              internal iteration = %i\n", cumul_iter);
  }
  free(soclcp->q);
  free(soclcp->coneIndex);
  free(soclcp);


  if (internalsolver_options->internalSolvers != NULL)
    internalsolver_options->internalSolvers->dWork = NULL;
  dparam[SICONOS_VI_EG_DPARAM_RHO] = internalsolver_options->dparam[SICONOS_VI_EG_DPARAM_RHO];
  dparam[SICONOS_DPARAM_RESIDU] = error;
  iparam[SICONOS_IPARAM_ITER_DONE] = iter;

}



int fc3d_ACLMFixedPoint_setDefaultSolverOptions(SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the ACLMFP Solver\n");
  }

  options->solverId = SICONOS_FRICTION_3D_ACLMFP;
  options->numberOfInternalSolvers = 1;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 8;
  options->dSize = 8;
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  solver_options_nullify(options);

  options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
  options->iparam[SICONOS_FRICTION_3D_FP_ERROR_STRATEGY] =  SICONOS_FRICTION_3D_FP_ERROR_STRATEGY_ADAPTIVE;
  options->dparam[SICONOS_DPARAM_TOL] = 1e-4;
  options->internalSolvers = (SolverOptions *)malloc(sizeof(SolverOptions));

  soclcp_nsgs_setDefaultSolverOptions(options->internalSolvers);
  options->internalSolvers->iparam[SICONOS_IPARAM_MAX_ITER] =10000;
  return 0;
}
