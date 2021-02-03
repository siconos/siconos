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

#include <assert.h>                                   // for assert
#include <math.h>                                     // for fabs, pow
#include <stdio.h>                                    // for printf, fprintf
#include <stdlib.h>                                   // for exit, free, malloc
#include "FrictionContactProblem.h"                   // for FrictionContact...
#include "Friction_cst.h"                             // for SICONOS_FRICTIO...
#include "NumericsFwd.h"                              // for SolverOptions
#include "NumericsMatrix.h"                           // for NM_add_to_diag3
#include "SiconosBlas.h"                              // for cblas_daxpy
#include "SolverOptions.h"                            // for SolverOptions
#include "fc3d_Solvers.h"                             // for fc3d_set_intern...
#include "fc3d_compute_error.h"                       // for fc3d_compute_error
#include "fc3d_nonsmooth_Newton_AlartCurnier.h"       // for fc3d_nonsmooth_...
#include "fc3d_nonsmooth_Newton_FischerBurmeister.h"  // for fc3d_nonsmooth_...
#include "numerics_verbose.h"                         // for numerics_error
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "siconos_debug.h"                                    // for DEBUG_PRINTF

void fc3d_proximal(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options)
{
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;

  /* Number of contacts */
  int nc = problem->numberOfContacts;
  NumericsMatrix* M = problem->M;

  /* Dimension of the problem */
  int n = 3 * nc;

  /* Maximum number of iterations */
  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];
  /* Tolerance */
  double tolerance = dparam[SICONOS_DPARAM_TOL];
  double norm_q = cblas_dnrm2(nc*3, problem->q, 1);
  if(verbose > 0)
  {
    solver_options_print(options);
  }

  if(options->numberOfInternalSolvers < 1)
  {
    numerics_error("fc3d_proximal", "The PROX method needs options for the internal solvers, options[0].numberOfInternalSolvers should be >1");
  }
  SolverOptions *internalsolver_options = options->internalSolvers[0];


  /*****  PROXIMAL Iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;

  int isVariable = 1;
  double alpha = dparam[SICONOS_FRICTION_3D_PROXIMAL_DPARAM_ALPHA];
  double sigma= 0.0;
  double nu =0.0;



  if(iparam[SICONOS_FRICTION_3D_PROXIMAL_IPARAM_STRATEGY] ==
      SICONOS_FRICTION_3D_PROXIMAL_REGULARIZATION) /* Only regularization */
  {
    if(fabs(dparam[SICONOS_FRICTION_3D_PROXIMAL_DPARAM_ALPHA]) < 1e-12)
    {
      fprintf(stderr, "Numerics,  fc3d_proximal. Initial alpha parameters equal 0 \n");
      exit(EXIT_FAILURE);
    }
    else if(dparam[SICONOS_FRICTION_3D_PROXIMAL_DPARAM_ALPHA] < -1e-12)
    {
      internalsolver_options->dparam[SICONOS_DPARAM_TOL]=options->dparam[SICONOS_DPARAM_TOL];
      alpha = - dparam[3];
      isVariable = 0;
    }
    else
    {
      internalsolver_options->dparam[SICONOS_DPARAM_TOL]=options->dparam[SICONOS_DPARAM_TOL];
      isVariable = 1;
      alpha = dparam[SICONOS_FRICTION_3D_PROXIMAL_DPARAM_ALPHA];
    }
  }
  else if(iparam[SICONOS_FRICTION_3D_PROXIMAL_IPARAM_STRATEGY] ==
          SICONOS_FRICTION_3D_PROXIMAL_PROX)
  {
    if(fabs(dparam[SICONOS_FRICTION_3D_PROXIMAL_DPARAM_ALPHA]) < 1e-12)
    {
      fprintf(stderr, "Numerics,  fc3d_proximal. Initial alpha parameters equal 0 \n");
      exit(EXIT_FAILURE);
    }
    else if(dparam[SICONOS_FRICTION_3D_PROXIMAL_DPARAM_ALPHA] < -1e-12)
    {
      alpha = - dparam[SICONOS_FRICTION_3D_PROXIMAL_DPARAM_ALPHA];
      isVariable = 0;
    }
    else
    {
      isVariable = 1;

      sigma = options->dparam[SICONOS_FRICTION_3D_PROXIMAL_DPARAM_SIGMA];
      nu = options->dparam[SICONOS_FRICTION_3D_PROXIMAL_DPARAM_NU];

      fc3d_compute_error(problem, reaction, velocity, tolerance, options, norm_q,  &error);

      alpha = sigma*pow(error,nu);
      DEBUG_PRINTF("Initial error = %g\n",error);
      if(error < tolerance) hasNotConverged = 0;
    }
  }
  else
    numerics_error("fc3d_proximal", "Proximal strategy is unknown");


  double * reactionold = (double *)malloc(n * sizeof(double));
  cblas_dcopy(n, reaction, 1, reactionold, 1);

  internalSolverPtr internalsolver =0;
  options->iparam[SICONOS_FRICTION_3D_PROXIMAL_IPARAM_CUMULATIVE_ITER_DONE]= 0;

  if(internalsolver_options)
  {
    if(internalsolver_options->solverId == SICONOS_FRICTION_3D_NSGS)
    {
      internalsolver = &fc3d_nsgs;
    }
    else if(internalsolver_options->solverId == SICONOS_FRICTION_3D_DSFP)
    {
      internalsolver = &fc3d_DeSaxceFixedPoint;
    }
    else if(internalsolver_options->solverId == SICONOS_FRICTION_3D_EG)
    {
      internalsolver = &fc3d_ExtraGradient;
    }
    else if(internalsolver_options->solverId == SICONOS_FRICTION_3D_NSN_AC)
    {
      internalsolver = &fc3d_nonsmooth_Newton_AlartCurnier;
    }
    else if(internalsolver_options->solverId == SICONOS_FRICTION_3D_NSN_FB)
    {
      internalsolver = &fc3d_nonsmooth_Newton_FischerBurmeister;
    }
    else
      numerics_error("fc3d_proximal", "unknown internal solver");
  }
  else
  {
    numerics_error("fc3d_proximal", "The PROX method needs options for the internal solvers, soptions->internalSolvers should be different from NULL");
  }



  DEBUG_PRINTF("isVariable = %i\n",isVariable);
  DEBUG_PRINTF("options->iparam[SICONOS_FRICTION_3D_PROXIMAL_IPARAM_STRATEGY] = %i\n",options->iparam[SICONOS_FRICTION_3D_PROXIMAL_IPARAM_STRATEGY]);

  if(verbose > 0)
  {
    printf("---- FC3D - PROXIMAL - Start with alpha = %12.8e\n", alpha);

    if(isVariable)
      printf("---- FC3D - PROXIMAL - Variable alpha strategy\n");
    else
      printf("---- FC3D - PROXIMAL - Fixed alpha strategy\n");

    if(options->iparam[SICONOS_FRICTION_3D_PROXIMAL_IPARAM_STRATEGY]==SICONOS_FRICTION_3D_PROXIMAL_REGULARIZATION)
      printf("---- FC3D - PROXIMAL - Only regularization \n");
  }

  if(iparam[SICONOS_FRICTION_3D_PROXIMAL_IPARAM_STRATEGY] ==
      SICONOS_FRICTION_3D_PROXIMAL_REGULARIZATION)
  {

    DEBUG_PRINTF("hasNotConverged = %i \n",hasNotConverged);
    while((iter < itermax) && (hasNotConverged > 0))
    {
      ++iter;
      cblas_dcopy(n, reaction, 1, reactionold, 1);
      /* add proximal regularization on q */
      cblas_daxpy(n, -alpha, reactionold, 1, problem->q, 1) ;

      /* add proximal regularization on M */
      NM_add_to_diag3(M, alpha);
      numerics_printf("---- FC3D - PROXIMAL - alpha = %8.4e\n",alpha);

      fc3d_set_internalsolver_tolerance(problem,options,internalsolver_options, error);
      DEBUG_PRINTF("internal solver tolerance = %21.8e \n",internalsolver_options->dparam[SICONOS_DPARAM_TOL]);

      /* call internal solver */
      (*internalsolver)(problem, reaction, velocity, info, internalsolver_options);

      if(verbose >0 && *info)
        printf("---- FC3D - PROXIMAL - internalsolver no convergence\n");

      int iter_internalsolver = internalsolver_options->iparam[SICONOS_IPARAM_ITER_DONE];
      options->iparam[SICONOS_FRICTION_3D_PROXIMAL_IPARAM_CUMULATIVE_ITER_DONE ] +=iter_internalsolver;
      DEBUG_PRINTF("iter_internalsolver = %i\n",iter_internalsolver);
      DEBUG_PRINTF("info = %i\n",*info);
      DEBUG_PRINTF("options->iparam[SICONOS_IPARAM_ITER_DONE] = %i\n",options->iparam[SICONOS_IPARAM_ITER_DONE]);


      /* **** Criterium convergence **** */

      /* substract proximal regularization on q */
      cblas_daxpy(n, alpha, reactionold, 1, problem->q, 1) ;
      /* substract proximal regularization on M */
      NM_add_to_diag3(M, -alpha);

      fc3d_compute_error(problem, reaction, velocity, tolerance, options, norm_q, &error);
      if(verbose > 0)
        printf("---- FC3D - PROXIMAL - Iteration %i Residual = %14.7e with alpha = %12.8e\n\n", iter, error, alpha);

      /* update alpha */
      if(isVariable)
      {
        if(*info)
        {
          alpha = alpha*5;
        }
        else
        {
          alpha = alpha/10.0;
        }
      }
      DEBUG_PRINTF("alpha = %8.4e\n",alpha);

      if(options->callback)
      {
        options->callback->collectStatsIteration(options->callback->env, nc * 3,
            reaction, velocity,
            error, NULL);
      }


      if(error < tolerance) hasNotConverged = 0;
      *info = hasNotConverged;
    }
  }
  /* { */
  /*   double alpha_old = alpha; */
  /*   while ((iter < itermax) && (hasNotConverged > 0)) */
  /*   { */
  /*     ++iter; */
  /*     //Add proximal regularization on M */
  /*     //This code looked weird before. Do we add (alpha - alpha_old) ? */
  /*     double pert = alpha; */
  /*     if (iter > 0) */
  /*     { */
  /*       pert -= alpha_old; */
  /*     } */
  /*     NM_add_to_diag3(M, pert); */

  /*     DEBUG_PRINTF("internal solver tolerance = %21.8e \n",internalsolver_options->dparam[SICONOS_DPARAM_TOL]); */

  /*     (*internalsolver)(problem, reaction , velocity , info , internalsolver_options); */

  /*     /\* **** Criterium convergence **** *\/ */
  /*     fc3d_compute_error(problem, reaction , velocity, tolerance, options, norm_q, &error); */

  /*     int iter_internalsolver = internalsolver_options->iparam[SICONOS_IPARAM_ITER_DONE]; */
  /*     options->iparam[SICONOS_FRICTION_3D_PROXIMAL_IPARAM_CUMULATIVE_ITER_DONE] */
  /*       +=iter_internalsolver; */

  /*     DEBUG_PRINTF("iter_internalsolver = %i\n",iter_internalsolver); */
  /*     DEBUG_PRINTF("info = %i\n",*info); */
  /*     DEBUG_PRINTF("options->iparam[SICONOS_IPARAM_ITER_DONE] = %i\n",options->iparam[SICONOS_IPARAM_ITER_DONE]); */
  /*     DEBUG_PRINTF("alpha = %8.4e\n",alpha); */

  /*     if (options->callback) */
  /*     { */
  /*       options->callback->collectStatsIteration(options->callback->env, nc * 3, */
  /*                                                reaction, velocity, */
  /*                                                error, NULL); */
  /*     } */

  /*     if (verbose > 0) */
  /*       printf("---- FC3D - PROXIMAL - Iteration %i Residual = %14.7e with alpha = %12.8e\n\n", iter, error, alpha); */
  /*     if (isVariable) */
  /*     { */
  /*       alpha_old =alpha; */
  /*       alpha = alpha*10; */
  /*     } */
  /*     if (error < tolerance) hasNotConverged = 0; */
  /*     *info = hasNotConverged; */
  /*   } */
  /* } */
  else if(iparam[SICONOS_FRICTION_3D_PROXIMAL_IPARAM_STRATEGY] ==
          SICONOS_FRICTION_3D_PROXIMAL_PROX) // Real PROX iparam[9] == 0
  {

    DEBUG_PRINTF("hasNotConverged = %i \n",hasNotConverged);
    while((iter < itermax) && (hasNotConverged > 0))
    {
      ++iter;
      cblas_dcopy(n, reaction, 1, reactionold, 1);
      /* add proximal regularization on q */
      cblas_daxpy(n, -alpha, reactionold, 1, problem->q, 1) ;

      /* add proximal regularization on M */
      NM_add_to_diag3(M, alpha);


      fc3d_set_internalsolver_tolerance(problem,options,internalsolver_options, error);
      DEBUG_PRINTF("internal solver tolerance = %21.8e \n",internalsolver_options->dparam[SICONOS_DPARAM_TOL]);

      /* call internal solver */
      (*internalsolver)(problem, reaction, velocity, info, internalsolver_options);

      if(verbose >0 && *info)
        printf("---- FC3D - PROXIMAL - internalsolver no convergence\n");

      int iter_internalsolver = internalsolver_options->iparam[SICONOS_IPARAM_ITER_DONE];
      options->iparam[SICONOS_FRICTION_3D_PROXIMAL_IPARAM_CUMULATIVE_ITER_DONE ] +=iter_internalsolver;
      DEBUG_PRINTF("iter_internalsolver = %i\n",iter_internalsolver);
      DEBUG_PRINTF("info (internal solver)= %i\n",*info);




      /* **** Criterium convergence **** */

      /* substract proximal regularization on q */
      cblas_daxpy(n, alpha, reactionold, 1, problem->q, 1) ;
      /* substract proximal regularization on M */
      NM_add_to_diag3(M, -alpha);

      fc3d_compute_error(problem, reaction, velocity, tolerance, options, norm_q, &error);

      /* update alpha */
      if(isVariable)
      {
        if(*info)
        {
          alpha = alpha*10;
        }
        else
        {
          alpha = sigma*pow(error,nu);
        }
      }
      DEBUG_PRINTF("alpha = %8.4e\n",alpha);

      if(options->callback)
      {
        options->callback->collectStatsIteration(options->callback->env, nc * 3,
            reaction, velocity,
            error, NULL);
      }

      if(verbose > 0)
        printf("---- FC3D - PROXIMAL - Iteration %i Residual = %14.7e with alpha = %12.8e\n\n", iter, error, alpha);

      if(error < tolerance) hasNotConverged = 0;
      *info = hasNotConverged;
    }
  }
  else
    numerics_error("fc3d_proximal", "Proximal strategy is unknown");


  if(verbose > 0)
  {
    printf("---- FC3D - PROXIMAL - # Iteration %i Final Residual = %14.7e  \n", iter, error);
    printf("---- FC3D - PROXIMAL - # Iteration of internal solver %i \n", options->iparam[SICONOS_FRICTION_3D_PROXIMAL_IPARAM_CUMULATIVE_ITER_DONE]);
  }

  iparam[SICONOS_IPARAM_ITER_DONE] = iter;
  dparam[SICONOS_DPARAM_RESIDU] = error;

  free(reactionold);

}


void fc3d_proximal_set_default(SolverOptions* options)
{
  options->iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] =
    SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE ;
  /* no overrelaxation by default */
  //   options->iparam[SICONOS_FRICTION_3D_PROXIMAL_IPARAM_RELAXATION] = 0;

  /* fixed regularization or proximal */
  options->iparam[SICONOS_FRICTION_3D_PROXIMAL_IPARAM_STRATEGY ] =
    SICONOS_FRICTION_3D_PROXIMAL_PROX;
  options->dparam[SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO] =10.0;
  /* default value for proximal parameter alpha; */
  options->dparam[SICONOS_FRICTION_3D_PROXIMAL_DPARAM_ALPHA] = 1.e4;
  /* default value for sigma; */
  options->dparam[SICONOS_FRICTION_3D_PROXIMAL_DPARAM_SIGMA] = 5.0;
  /* default value for nu; */
  options->dparam[SICONOS_FRICTION_3D_PROXIMAL_DPARAM_NU] = 1.0;
  /* default value for relaxation parameter omega */
  // options->dparam[SICONOS_FRICTION_3D_PROXIMAL_DPARAM_RELAXATION] = 1.5;


  assert(options->numberOfInternalSolvers == 1);
  options->internalSolvers[0] = solver_options_create(SICONOS_FRICTION_3D_NSN_AC);
}
