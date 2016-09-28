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
#include "NCP_Solvers.h"
#include "SiconosBlas.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
//#define DEBUG_STDOUT */
//#define DEBUG_MESSAGES */
#include "debug.h"
#include <math.h>
#include "numerics_verbose.h"

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
  int itermax = iparam[0];
  /* Tolerance */
  double tolerance = dparam[0];
  double normq = cblas_dnrm2(nc*3 , problem->q , 1);
  if (verbose > 0){
    solver_options_print(options);
  }

  if (options->numberOfInternalSolvers < 1)
  {
    numerics_error("fc3d_proximal", "The PROX method needs options for the internal solvers, options[0].numberOfInternalSolvers should be >1");
  }
  SolverOptions *internalsolver_options = options->internalSolvers;


  /*****  PROXIMAL Iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;
  
  int isVariable = 1;
  double alpha = dparam[3];
  double sigma= 0.0;
  double nu =0.0;

  

  if (iparam[9]) /* Only regularization */
  {
    if (fabs(dparam[3]) < 1e-12)
    {
      fprintf(stderr, "Numerics,  fc3d_proximal. Initial alpha parameters equal 0 \n");
      exit(EXIT_FAILURE);
    }
    else if (dparam[3] < -1e-12)
    {
      internalsolver_options->dparam[0]=options->dparam[0];
      alpha = - dparam[3];
      isVariable = 0;
    }
    else
    {
      internalsolver_options->dparam[0]=options->dparam[0];
      isVariable = 1;
      alpha = dparam[3];
    }
  }
  else
  {
    if (fabs(dparam[3]) < 1e-12)
    {
      fprintf(stderr, "Numerics,  fc3d_proximal. Initial alpha parameters equal 0 \n");
      exit(EXIT_FAILURE);
    }
    else if (dparam[3] < -1e-12)
    {
      alpha = - dparam[3];
      isVariable = 0;
    }
    else
    {
      isVariable = 1;
      /* parameters to set the value of alpha */
      sigma = options->dparam[4];
      nu = options->dparam[5];
      fc3d_compute_error(problem, reaction , velocity, tolerance, options, normq,  &error);
      internalsolver_options->dparam[0] = options->dparam[0];
      internalsolver_options->dparam[0] = error;
      alpha = sigma*pow(error,nu);
      DEBUG_PRINTF("error = %g\n",error);
    
      DEBUG_PRINTF("tolerance = %g\n",tolerance);
    
      if (error < tolerance) hasNotConverged = 0;
      DEBUG_PRINTF(" hasNotConverged = %i \n",hasNotConverged);
    }
  }
  
  DEBUG_PRINTF("options->iparam[2] = %i\n",options->iparam[2]);
  DEBUG_PRINTF("options->iparam[3] = %i\n",options->iparam[3]);
  DEBUG_PRINTF("options->dparam[4] = %i\n",options->dparam[4]);
  DEBUG_PRINTF("options->dparam[5] = %i\n",options->dparam[5]);

  double * reactionold = (double *)malloc(n * sizeof(double));
  cblas_dcopy(n , reaction , 1 , reactionold , 1);

  internalSolverPtr internalsolver =0;
  int iter_iparam =7;
  options->iparam[6]= 0;

  if (internalsolver_options) /* we assume that the solver option has been coherently set*/
  {
    if (internalsolver_options->solverId == SICONOS_FRICTION_3D_NSGS)
    {
      internalsolver = &fc3d_nsgs;
    }
    else if (internalsolver_options->solverId == SICONOS_FRICTION_3D_DeSaxceFixedPoint)
    {
      internalsolver = &fc3d_DeSaxceFixedPoint;
    }
    else if (internalsolver_options->solverId == SICONOS_FRICTION_3D_EG)
    {
      internalsolver = &fc3d_ExtraGradient;
    }
    else if (internalsolver_options->solverId == SICONOS_FRICTION_3D_NSN_AC)
    {
      internalsolver = &fc3d_nonsmooth_Newton_AlartCurnier;
      iter_iparam =1;
    }
    else if (internalsolver_options->solverId == SICONOS_FRICTION_3D_NSN_FB)
    {

      internalsolver = &fc3d_nonsmooth_Newton_FischerBurmeister;
      iter_iparam =1;
    }
  }
  else
  {
    numerics_error("fc3d_proximal", "The PROX method needs options for the internal solvers, soptions->internalSolvers should be diffrent from NULL");
  }




  
  DEBUG_PRINTF("options->iparam[2] = %i\n",options->iparam[2]);
  DEBUG_PRINTF("options->iparam[3] = %i\n",options->iparam[3]);


  DEBUG_PRINTF("isVariable = %i\n",isVariable);
  DEBUG_PRINTF("(only regularization) options->iparam[9] = %i\n",options->iparam[9]);

  if (verbose > 0){
    printf("------------------------ FC3D - PROXIMAL - Start with alpha = %12.8e\n", alpha);
    if (isVariable)
      printf("------------------------ FC3D - PROXIMAL - Variable alpha strategy\n");
    else
      printf("------------------------ FC3D - PROXIMAL - Fixed alpha strategy\n");
    if (options->iparam[9])
      printf("------------------------ FC3D - PROXIMAL - Only regularization \n");
  }

  if (iparam[9]){
    double alpha_old;
    while ((iter < itermax) && (hasNotConverged > 0))
    {
      ++iter;
      //Add proximal regularization on M
      if (M->storageType == 0)
      {
        for (int i = 0 ; i < n ; i++) M->matrix0[i + i * n] += alpha;
      }
      else if (M->storageType == 1)
      {
        for (int ic = 0 ; ic < nc ; ic++)
        {
          int diagPos = getDiagonalBlockPos(M->matrix1, ic);
          if (iter >0 )
            for (int i = 0 ; i < 3 ; i++) M->matrix1->block[diagPos][i + 3 * i] -= alpha_old ;
          for (int i = 0 ; i < 3 ; i++) M->matrix1->block[diagPos][i + 3 * i] += alpha ;
        }
      }

      DEBUG_PRINTF("internal solver tolerance = %21.8e \n",internalsolver_options->dparam[0]);

      (*internalsolver)(problem, reaction , velocity , info , internalsolver_options);

      /* **** Criterium convergence **** */
      fc3d_compute_error(problem, reaction , velocity, tolerance, options, normq, &error);

      int iter_internalsolver = internalsolver_options->iparam[iter_iparam];
      options->iparam[6] +=iter_internalsolver;

      DEBUG_PRINTF("iter_internalsolver = %i\n",iter_internalsolver);
      DEBUG_PRINTF("info = %i\n",*info);
      DEBUG_PRINTF("options->iparam[1] = %i\n",options->iparam[1]);
      DEBUG_PRINTF("alpha = %8.4e\n",alpha);

      if (options->callback)
      {
        options->callback->collectStatsIteration(options->callback->env, nc * 3,
                                                 reaction, velocity,
                                                 error, NULL);
      }

      if (verbose > 0)
        printf("------------------------ FC3D - PROXIMAL - Iteration %i Residual = %14.7e with alpha = %12.8e\n\n", iter, error, alpha);
      if (isVariable)
      {
        alpha_old =alpha;
        alpha = alpha*10;
      }
      if (error < tolerance) hasNotConverged = 0;
      *info = hasNotConverged;
    }
  }
  else // Real PROX iparam[9] == 0
  {

    DEBUG_PRINTF(" hasNotConverged = %i \n",hasNotConverged);
    while ((iter < itermax) && (hasNotConverged > 0))
    {
      ++iter;
      cblas_dcopy(n , reaction , 1 , reactionold , 1);
      //Add proximal regularization on q
      cblas_daxpy(n, -alpha, reactionold, 1, problem->q , 1) ;

      //Add proximal regularization on M
      if (M->storageType == 0)
      {
        for (int i = 0 ; i < n ; i++) M->matrix0[i + i * n] += alpha;
      }
      else if (M->storageType == 1)
      {
        for (int ic = 0 ; ic < nc ; ic++)
        {
          int diagPos = getDiagonalBlockPos(M->matrix1, ic);
          for (int i = 0 ; i < 3 ; i++) M->matrix1->block[diagPos][i + 3 * i] += alpha ;
        }
      }
      // internal solver for the regularized problem
      /*       fc3d_nsgs(problem, reaction , velocity , info , internalsolver_options); */

      /* internalsolver_options->dparam[0] = max(error/10.0, options->dparam[0]); */
      //internalsolver_options->dparam[0] = options->dparam[0];

      internalsolver_options->dparam[0] = alpha*error;
      DEBUG_PRINTF("internal solver tolerance = %21.8e \n",internalsolver_options->dparam[0]);
      (*internalsolver)(problem, reaction , velocity , info , internalsolver_options);

      /* **** Criterium convergence **** */

      //substract proximal regularization on q

      cblas_daxpy(n, alpha, reactionold, 1, problem->q, 1) ;

      //substract proximal regularization on M
      if (M->storageType == 0)
      {
        for (int i = 0 ; i < n ; i++) M->matrix0[i + i * n] -= alpha;
      }
      else if (M->storageType == 1)
      {
        for (int ic = 0 ; ic < nc ; ic++)
        {
          int diagPos = getDiagonalBlockPos(M->matrix1, ic);
          for (int i = 0 ; i < 3 ; i++) M->matrix1->block[diagPos][i + 3 * i] -= alpha ;
        }
      }

      fc3d_compute_error(problem, reaction , velocity, tolerance, options, normq, &error);
      //update the alpha with respect to the number of internal iterations.

      int iter_internalsolver = internalsolver_options->iparam[iter_iparam];
      options->iparam[6] +=iter_internalsolver;
      DEBUG_PRINTF("iter_internalsolver = %i\n",iter_internalsolver);
      DEBUG_PRINTF("info = %i\n",*info);
      DEBUG_PRINTF("options->iparam[1] = %i\n",options->iparam[1]);
      if (verbose >0 && *info)
        printf("------------------------ FC3D - PROXIMAL - internalsolver not converged\n");
      if (isVariable)
      {
        if (*info){
          alpha = alpha*10;
        }
        else{
          alpha = sigma*pow(error,nu);
        } 
      }

      DEBUG_PRINTF("alpha = %8.4e\n",alpha);

      if (options->callback)
      {
        options->callback->collectStatsIteration(options->callback->env, nc * 3,
                                                 reaction, velocity,
                                                 error, NULL);
      }

      if (verbose > 0)
        printf("------------------------ FC3D - PROXIMAL - Iteration %i Residual = %14.7e with alpha = %12.8e\n\n", iter, error, alpha);

      if (error < tolerance) hasNotConverged = 0;
      *info = hasNotConverged;
    }
  }




  if (verbose > 0)
  {
    printf("------------------------ FC3D - PROXIMAL - # Iteration %i Final Residual = %14.7e  \n", iter, error);
    printf("------------------------ FC3D - PROXIMAL - # Iteration of internal solver %i \n", options->iparam[6]);
  }

  iparam[7] = iter;
  dparam[0] = tolerance;
  dparam[1] = error;

  free(reactionold);

}


int fc3d_proximal_setDefaultSolverOptions(SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the PROX Solver\n");
  }

  /*  strcpy(options->solverName,"PROX");*/
  options->solverId = SICONOS_FRICTION_3D_PROX;
  options->numberOfInternalSolvers = 1;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 10;
  options->dSize = 10;
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  options->dWork = NULL;
  solver_options_nullify(options);
  options->iparam[0] = 1000;

  options->iparam[2] = 5;    /* Default mimimun iteration of the internal
                                solver for decreasing alpha */
  options->iparam[3] = 15;   /* Default maximum iteration of the internal
                                solver for increasing alpha */
  options->iparam[8] = 0;    /* no overrelaxation by default */
  options->iparam[9] = 0;    /* fixed regularization and not proximal */

  options->iparam[3] = 15;   /* Default maximum iteration of the internal
                                solver for increasing alpha */



  options->dparam[0] = 1e-4;
  options->dparam[3] = 1.e4; /* default value for proximal parameter alpha; */
  options->dparam[4] = 5.0;  /* default value for sigma; */
  options->dparam[5] = 1.0;  /* default value for nu; */
  options->dparam[8] = 1.5;  /* default value for relaxation parameter omega */

  options->internalSolvers = (SolverOptions *)malloc(sizeof(SolverOptions));
  fc3d_nonsmooth_Newton_AlartCurnier_setDefaultSolverOptions(options->internalSolvers);

  return 0;
}
