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
#include <assert.h>                  // for assert
#include <math.h>                    // for fmax
#include <stdio.h>                   // for printf, NULL
#include <stdlib.h>                  // for malloc
#include "FrictionContactProblem.h"  // for FrictionContactProblem
#include "Friction_cst.h"            // for SICONOS_FRICTION_3D_NSGS, SICONO...
#include "NumericsFwd.h"             // for SolverOptions, FrictionContactPr...
#include "SiconosBlas.h"             // for cblas_dnrm2
#include "SolverOptions.h"           // for SolverOptions, solver_options_cr...
#include "fc3d_Solvers.h"            // for fc3d_set_internalsolver_tolerance
#include "fc3d_compute_error.h"      // for fc3d_compute_error
#include "numerics_verbose.h"        // for numerics_error, verbose

void fc3d_TrescaFixedPoint(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options)
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
  double norm_q = cblas_dnrm2(nc*3, problem->q, 1);


  if(options->numberOfInternalSolvers < 1)
  {
    numerics_error("fc3d_TrescaFixedpoint", "The Tresca Fixed Point method needs options for the internal solvers; please check your options.");
  }

  SolverOptions * internalsolver_options = options->internalSolvers[0];

  if(verbose) solver_options_print(options);

  /*****  Fixed Point Iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;

  internalSolverPtr internalsolver;
  // dwork will be used to save friction threshold
  options->dWork = (double *) malloc(nc * sizeof(double));
  options->dWorkSize = nc;
  double * mu = options->dWork;
  // Warning : same dwork for current and internal solver !!
  internalsolver_options->dWork = mu;


  if(internalsolver_options->solverId == SICONOS_FRICTION_3D_NSGS)
  {
    if(verbose > 0)
      printf(" ========================== Call NSGS solver for Friction-Contact 3D problem ==========================\n");
    internalsolver = &fc3d_nsgs;
  }
  else if(internalsolver_options->solverId == SICONOS_FRICTION_3D_CONVEXQP_PG_CYLINDER)
  {
    if(verbose > 0)
      printf(" ========================== Call ConvexQP PG solver for Friction-Contact 3D problem ==========================\n");
    internalsolver = &fc3d_ConvexQP_ProjectedGradient_Cylinder;
  }
  else if(internalsolver_options->solverId == SICONOS_FRICTION_3D_VI_FPP_Cylinder)
  {
    if(verbose > 0)
      printf(" ========================== Call VI FPP solver for Friction-Contact 3D problem ==========================\n");
    internalsolver = &fc3d_VI_FixedPointProjection_Cylinder;
  }
  else
  {
    numerics_error("fc3d_TrescaFixedpoint", "Unknown internal solver.");
  }

  int cumul_internal=0;

  while((iter < itermax) && (hasNotConverged > 0))
  {
    ++iter;

    /* Compute the value of the initial value friction threshold*/
    for(int ic = 0 ; ic < nc ; ic++) mu[ic] = fmax(0.0, problem->mu[ic] *  reaction [ic * 3]);

    if(verbose>0)
      printf("norm of mu = %10.5e \n", cblas_dnrm2(nc, mu, 1));

    fc3d_set_internalsolver_tolerance(problem,options,internalsolver_options, error);

    (*internalsolver)(problem, reaction, velocity, info, internalsolver_options);

    cumul_internal += internalsolver_options->iparam[SICONOS_IPARAM_ITER_DONE];

    /* **** Criterium convergence **** */

    fc3d_compute_error(problem, reaction, velocity, tolerance, options, norm_q,  &error);

    if(options->callback)
    {
      options->callback->collectStatsIteration(options->callback->env, nc * 3,
          reaction, velocity, error, NULL);
    }

    if(error < tolerance) hasNotConverged = 0;
    *info = hasNotConverged;

    if(verbose > 0)
    {
      if(hasNotConverged)
      {
        printf("--------------- FC3D - TFP - Iteration %i error = %14.7e > %10.5e\n", iter, error, tolerance);
      }
      else
      {
        printf("--------------- FC3D - TFP - Iteration %i error = %14.7e < %10.5e\n", iter, error, tolerance);
        printf("--------------- FC3D - TFP - #              Internal iteration = %i\n", cumul_internal);
      }
    }
  }
  mu = NULL;
  internalsolver_options->dWork = NULL;
  dparam[SICONOS_DPARAM_RESIDU] = error;
  iparam[SICONOS_IPARAM_ITER_DONE] = iter;

}



void fc3d_tfp_set_default(SolverOptions* options)
{

  options->iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] =  SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE;
  options->dparam[SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO] =10.0;

  // internal solver
  assert(options->numberOfInternalSolvers == 1);
  options->internalSolvers[0] = solver_options_create(SICONOS_FRICTION_3D_NSGS);
  options->internalSolvers[0]->iparam[SICONOS_IPARAM_MAX_ITER]=1000;

  // internal solver of the internal solver
  assert(options->internalSolvers[0]->numberOfInternalSolvers == 1);
  options->internalSolvers[0]->internalSolvers[0] = solver_options_create(SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinderWithLocalIteration);
  options->internalSolvers[0]->internalSolvers[0]->iparam[SICONOS_IPARAM_MAX_ITER] = 50;
  options->internalSolvers[0]->internalSolvers[0]->dparam[SICONOS_DPARAM_TOL] = 1e-14;
}
