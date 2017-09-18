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
#include "SiconosBlas.h"
#include "Friction_cst.h"
#include "numerics_verbose.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>



void fc3d_FixedPoint_set_internalsolver_tolerance(FrictionContactProblem* problem,
                                                  SolverOptions* options,
                                                  SolverOptions* internalsolver_options,
                                                  double error)
{
  int* iparam = options->iparam;
  if (iparam[SICONOS_FRICTION_3D_FP_ERROR_STRATEGY] == SICONOS_FRICTION_3D_FP_ERROR_STRATEGY_ADAPTIVE )
  {
    internalsolver_options->dparam[0] = fmax(error/10.0, options->dparam[0]/problem->numberOfContacts);
    numerics_printf("fc3d_FixedPoint_set_internalsolver_tolerance - Internal solver tolerance is set to %e\n", internalsolver_options->dparam[0] );
  }
  else if (iparam[SICONOS_FRICTION_3D_FP_ERROR_STRATEGY] == SICONOS_FRICTION_3D_FP_ERROR_STRATEGY_FRACTION)
  {
    internalsolver_options->dparam[0] = options->dparam[0]/2.0;
    numerics_printf("fc3d_FixedPoint_set_internalsolver_tolerance - Internal solver tolerance is set to %e\n", internalsolver_options->dparam[0] );
  }
  else if (iparam[SICONOS_FRICTION_3D_FP_ERROR_STRATEGY] == SICONOS_FRICTION_3D_FP_ERROR_STRATEGY_GIVEN_VALUE)
  {
    // We use the user value for the error of the local solver
    numerics_printf("fc3d_FixedPoint_set_internalsolver_tolerance - Internal solver tolerance is set to %e\n", internalsolver_options->dparam[0] );
  }
  else
  {
    numerics_error("fc3d_FixedPoint_set_internalsolver_tolerance","Unknown strategy for driving the tolerance");
  }

}

void fc3d_TrescaFixedPoint(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options)
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
    numerics_error("fc3d_TrescaFixedpoint", "The Tresca Fixed Point method needs options for the internal solvers, options[0].numberOfInternalSolvers should be >1");
  }

  SolverOptions * internalsolver_options = options->internalSolvers;

  if (verbose) solver_options_print(options);

  /*****  Fixed Point Iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;

  internalSolverPtr internalsolver;
  options->dWork = (double *) malloc(nc * sizeof(double));
  options->dWorkSize = nc;
  double * mu = options->dWork;
  internalsolver_options->dWork = options->dWork;


  if (internalsolver_options->solverId == SICONOS_FRICTION_3D_NSGS)
  {
    if (verbose > 0)
      printf(" ========================== Call NSGS solver for Friction-Contact 3D problem ==========================\n");
    internalsolver = &fc3d_nsgs;
  }
  else if (internalsolver_options->solverId == SICONOS_FRICTION_3D_ConvexQP_PG_Cylinder)
  {
    if (verbose > 0)
      printf(" ========================== Call ConvexQP PG solver for Friction-Contact 3D problem ==========================\n");
    internalsolver = &fc3d_ConvexQP_ProjectedGradient_Cylinder;
  }
 else if (internalsolver_options->solverId == SICONOS_FRICTION_3D_VI_FPP_Cylinder)
  {
    if (verbose > 0)
      printf(" ========================== Call VI FPP solver for Friction-Contact 3D problem ==========================\n");
    internalsolver = &fc3d_VI_FixedPointProjection_Cylinder;
  }
  else
  {
    numerics_error("fc3d_TrescaFixedpoint", "Unknown internal solver.");
  }

  int cumul_internal=0;

  while ((iter < itermax) && (hasNotConverged > 0))
  {
    ++iter;

    /* Compute the value of the initial value friction threshold*/
    for (int ic = 0 ; ic < nc ; ic++) mu[ic] = fmax(0.0, problem->mu[ic] *  reaction [ic * 3]);

    if (verbose>0)
      printf("norm of mu = %10.5e \n", cblas_dnrm2(nc , mu , 1));

    fc3d_FixedPoint_set_internalsolver_tolerance(problem,options,internalsolver_options, error);

    (*internalsolver)(problem, reaction , velocity , info , internalsolver_options);

    cumul_internal += internalsolver_options->iparam[SICONOS_IPARAM_ITER_DONE];

    /* **** Criterium convergence **** */

    fc3d_compute_error(problem, reaction , velocity, tolerance, options, norm_q,  &error);

    if (options->callback)
    {
      options->callback->collectStatsIteration(options->callback->env, nc * 3,
                                      reaction, velocity, error, NULL);
    }

    if (error < tolerance) hasNotConverged = 0;
    *info = hasNotConverged;

    if (verbose > 0)
    {
      if (hasNotConverged)
      {
        printf("----------------------------------- FC3D - TFP - Iteration %i error = %14.7e > %10.5e\n", iter, error, tolerance);
      }
      else
      {
        printf("----------------------------------- FC3D - TFP - Iteration %i error = %14.7e < %10.5e\n", iter, error, tolerance);
        printf("----------------------------------- FC3D - TFP - #              Internal iteration = %i\n", cumul_internal);
      }
    }
  }

  free(options->dWork);
  options->dWork = NULL;
  internalsolver_options->dWork = NULL;

  if (internalsolver_options->internalSolvers != NULL)
    internalsolver_options->internalSolvers->dWork = NULL;

  dparam[SICONOS_DPARAM_RESIDU] = error;
  iparam[SICONOS_IPARAM_ITER_DONE] = iter;

}



int fc3d_TrescaFixedPoint_setDefaultSolverOptions(SolverOptions* options)
{

  numerics_printf("fc3d_TrescaFixedPoint_setDefaultSolverOptions", "set default options"); 

  options->solverId = SICONOS_FRICTION_3D_TFP;
  options->numberOfInternalSolvers = 1;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 8;
  options->dSize = 8;
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  options->dWork = NULL;
  solver_options_nullify(options);

  options->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
  options->iparam[SICONOS_FRICTION_3D_FP_ERROR_STRATEGY] =  SICONOS_FRICTION_3D_FP_ERROR_STRATEGY_ADAPTIVE;
  options->dparam[SICONOS_DPARAM_TOL] = 1e-4;

  options->internalSolvers = (SolverOptions *)malloc(sizeof(SolverOptions));

  fc3d_nsgs_setDefaultSolverOptions(options->internalSolvers);
  options->internalSolvers->iparam[SICONOS_IPARAM_MAX_ITER]=1000;

  SolverOptions * subsubsolver = options->internalSolvers->internalSolvers;

  subsubsolver->iparam[SICONOS_IPARAM_MAX_ITER] = 50;
  subsubsolver->dparam[SICONOS_DPARAM_TOL] = 1e-14;

  subsubsolver->solverId = SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinderWithLocalIteration;

  return 0;
}
