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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#define VERBOSE_DEBUG
#include "Friction_cst.h"
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
  double normq = cblas_dnrm2(nc*3 , problem->q , 1);
 


  if (options->numberOfInternalSolvers < 1)
  {
    numericsError("fc3d_TrescaFixedpoint", "The Tresca Fixed Point method needs options for the internal solvers, options[0].numberOfInternalSolvers should be >1");
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

  internalSolverPtr internalsolver;
  options->dWork = (double *) malloc(nc * sizeof(double));
  options->dWorkSize = nc;
  double * mu = options->dWork;
  internalsolver_options->dWork = options->dWork;


  if (internalsolver_options->solverId == SICONOS_FRICTION_3D_NSGS)
  {
    if (verbose == 1)
      printf(" ========================== Call NSGS solver for Friction-Contact 3D problem ==========================\n");
    internalsolver = &fc3d_nsgs;
    //internalsolver_options->internalSolvers->dWork = options->dWork;
  }
  else if (internalsolver_options->solverId == SICONOS_FRICTION_3D_PGoC)
  {
    if (verbose == 1)
      printf(" ========================== Call PGoC solver for Friction-Contact 3D problem ==========================\n");
    internalsolver = &fc3d_ProjectedGradientOnCylinder;
  }
  else
  {
    fprintf(stderr, "Numerics, fc3d_TrescaFixedPoint failed. Unknown internal solver.\n");
    exit(EXIT_FAILURE);
  }



  int cumul_internal=0;

  while ((iter < itermax) && (hasNotConverged > 0))
  {
    ++iter;

    // internal solver for the regularized problem

    /* Compute the value of the initial value friction threshold*/
    for (int ic = 0 ; ic < nc ; ic++) mu[ic] = fmax(0.0, problem->mu[ic] *  reaction [ic * 3]);
    if (iparam[1] == 0 )
    {
      internalsolver_options->dparam[0] = fmax(error/10.0, options->dparam[0]/problem->numberOfContacts);
    }
    else if (iparam[1] ==1)
    {
      internalsolver_options->dparam[0] = options->dparam[0]/2.0;
    }
    else
    {
      fprintf(stderr, "Numerics, fc3d_TrescaFixedPoint failed. Unknown startegy for driving tolerence of internal.\n");
    exit(EXIT_FAILURE);
    }
    (*internalsolver)(problem, reaction , velocity , info , internalsolver_options);

    cumul_internal += internalsolver_options->iparam[7];

    /* **** Criterium convergence **** */

    fc3d_compute_error(problem, reaction , velocity, tolerance, options, normq,  &error);

    if (options->callback)
    {
      options->callback->collectStatsIteration(options->callback->env, nc * 3, 
                                      reaction, velocity, error, NULL);
    }

    if (verbose > 0)
      printf("------------------------ FC3D - TFP - Iteration %i Residual = %14.7e\n", iter, error);

    if (error < tolerance) hasNotConverged = 0;
    *info = hasNotConverged;
  }
  if (verbose > 0){
    printf("----------------------------------- FC3D - TFP - # Iteration %i Final Residual = %14.7e\n", iter, error);
    printf("----------------------------------- FC3D - TFP - #              Internal iteration = %i\n", cumul_internal);
  }
  free(options->dWork);
  options->dWork = NULL;
  internalsolver_options->dWork = NULL;

  if (internalsolver_options->internalSolvers != NULL)
    internalsolver_options->internalSolvers->dWork = NULL;
  dparam[0] = tolerance;
  dparam[1] = error;
  iparam[7] = iter;

}



int fc3d_TrescaFixedPoint_setDefaultSolverOptions(SolverOptions* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the TFP Solver\n");
  }


  options->solverId = SICONOS_FRICTION_3D_TFP;
  options->numberOfInternalSolvers = 1;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 8;
  options->dSize = 8;
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  solver_options_nullify(options);
  for (i = 0; i < 8; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->iparam[0] = 1000;
  options->dparam[0] = 1e-4;

  options->internalSolvers = (SolverOptions *)malloc(sizeof(SolverOptions));
  fc3d_nsgs_setDefaultSolverOptions(options->internalSolvers);
  options->internalSolvers->iparam[0]=1000;

  SolverOptions * subsubsolver = options->internalSolvers->internalSolvers;

  subsubsolver->iparam[0] = 50;
  subsubsolver->dparam[0] = 1e-14;

  subsubsolver->solverId = SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinderWithLocalIteration;

  return 0;
}
