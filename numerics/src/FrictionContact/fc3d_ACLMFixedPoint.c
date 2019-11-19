
/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#include <assert.h>                                       // for assert
#include <math.h>                                         // for sqrt
#include <stdio.h>                                        // for printf, fpr...
#include <stdlib.h>                                       // for free, malloc
#include "FrictionContactProblem.h"                       // for FrictionCon...
#include "Friction_cst.h"                                 // for SICONOS_FRI...
#include "NumericsFwd.h"                                  // for SolverOptions
#include "SOCLCP_Solvers.h"                               // for soclcp_VI_E...
#include "SOCLCP_cst.h"                                   // for SICONOS_SOC...
#include "SecondOrderConeLinearComplementarityProblem.h"  // for SecondOrder...
#include "SiconosBlas.h"                                  // for cblas_dnrm2
#include "SolverOptions.h"                                // for SolverOptions
#include "VI_cst.h"                                       // for SICONOS_VI_...
#define DEBUG_MESSAGES
#define DEBUG_STDOUT
#include "debug.h"  // lines 32-32
#include "fc3d_Solvers.h"                                 // for fc3d_set_in...
#include "fc3d_compute_error.h"                           // for fc3d_comput...
#include "numerics_verbose.h"                             // for verbose


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
  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];
  /* Tolerance */
  double tolerance = dparam[SICONOS_DPARAM_TOL];
  double norm_q = cblas_dnrm2(nc*3 , problem->q , 1);



  if (options->numberOfInternalSolvers < 1)
  {
    numerics_error("fc3d_ACLMFixedpoint", "The ACLM Fixed Point method needs options for the internal solvers, please check your options.");
  }

  SolverOptions * internalsolver_options = options->internalSolvers[0];

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

    fc3d_set_internalsolver_tolerance(problem,options,internalsolver_options, error);


    (*internalsolver)(soclcp, reaction , velocity , info , internalsolver_options);
    cumul_iter +=  internalsolver_options->iparam[SICONOS_IPARAM_ITER_DONE];
    /* **** Criterium convergence **** */

    fc3d_compute_error(problem, reaction , velocity, tolerance, options, norm_q, &error);

    if (options->callback)
    {
      options->callback->collectStatsIteration(options->callback->env, nc * 3,
                                      reaction, velocity, error, NULL);
    }

    if (verbose > 0)
      printf("---- FC3D - ACLMFP - Iteration %i Residual = %14.7e\n", iter, error);

    if (error < tolerance) hasNotConverged = 0;
    *info = hasNotConverged;
  }
  if (verbose > 0)
  {
    printf("--------------- FC3D - ACLMFP - # Iteration %i Final Residual = %14.7e\n", iter, error);
    printf("--------------- FC3D - ACLMFP - #              internal iteration = %i\n", cumul_iter);
  }
  free(soclcp->q);
  free(soclcp->coneIndex);
  free(soclcp);

  dparam[SICONOS_VI_DPARAM_RHO] = internalsolver_options->dparam[SICONOS_VI_DPARAM_RHO];
  dparam[SICONOS_DPARAM_RESIDU] = error;
  iparam[SICONOS_IPARAM_ITER_DONE] = iter;

}



void fc3d_aclmfp_set_options(SolverOptions* options)
{
  options->iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] =  SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE;
  options->dparam[SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO] =10.0;
  
  assert(options->numberOfInternalSolvers == 1);
  options->internalSolvers[0] = solver_options_create(SICONOS_SOCLCP_NSGS);
  options->internalSolvers[0]->iparam[SICONOS_IPARAM_MAX_ITER] =10000;
}
