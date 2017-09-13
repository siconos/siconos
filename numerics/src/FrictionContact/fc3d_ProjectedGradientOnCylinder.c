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

#include "projectionOnCylinder.h"
#include "fc3d_Solvers.h"
#include "fc3d_compute_error.h"
#include "SiconosBlas.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Friction_cst.h"
#include "numerics_verbose.h"

//#define VERBOSE_DEBUG

void fc3d_ProjectedGradientOnCylinder(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options)
{
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;
  /* Number of contacts */
  int nc = problem->numberOfContacts;
  double* q = problem->q;
  NumericsMatrix* M = problem->M;
  /* Dimension of the problem */
  int n = 3 * nc;
  /* Maximum number of iterations */
  int itermax = iparam[0];
  /* Tolerance */
  double tolerance = dparam[0];

  /*****  Projected Gradient iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;
  int contact; /* Number of the current row of blocks in M */
  int nLocal = 3;
  
  double rho = 0.0;
  double rhomin = 0.0;
  int isVariable = 0;

  if (dparam[SICONOS_FRICTION_3D_PGOC_RHO] > 0.0)
  {
    rho = dparam[SICONOS_FRICTION_3D_PGOC_RHO];
  }
  else
  {
    /* Variable step in fixed*/
    isVariable = 1;
    if (verbose > 0)
       printf("Variable step (line search) in Projected Gradient iterations\n");
    rho = -dparam[SICONOS_FRICTION_3D_PGOC_RHO];
    rhomin = dparam[SICONOS_FRICTION_3D_PGOC_RHOMIN];
  }

  if (rho == 0.0)
    numerics_error("fc3d_ProjectedGradientOnCylinder", "dparam[SICONOS_FRICTION_3D_PGOC_RHO] must be nonzero");
 

  
  double * reaction_k;
  double * direction;
  double * velocity_k;
  if (isVariable)
  {
    reaction_k = (double *)malloc(n * sizeof(double));
    direction = (double *)malloc(n * sizeof(double));
    velocity_k = (double *)malloc(n * sizeof(double));
  }
  double alpha = 1.0;
  double beta = 1.0;

  double norm_q = cblas_dnrm2(nc*3 , problem->q , 1);

  if (!isVariable)
  {
    double minusrho  = -1.0*rho;

    while ((iter < itermax) && (hasNotConverged > 0))
    {
      ++iter;

      /* q --> velocity */
      cblas_dcopy(n , q , 1 , velocity, 1);

      /* M r + q --> velocity */
      NM_gemv(alpha, M, reaction, beta, velocity);

      /* projection for each contact of reaction - rho * velocity  */
      cblas_daxpy(n, minusrho, velocity, 1, reaction , 1);
      for (contact = 0 ; contact < nc ; ++contact)
        projectionOnCylinder(&reaction[ contact * nLocal],
                             options->dWork[contact]);

      /* **** Criterium convergence **** */
      fc3d_Tresca_compute_error(problem, reaction , velocity, tolerance, options, norm_q, &error);

      if (options->callback)
      {
        options->callback->collectStatsIteration(options->callback->env, nc * 3,
                                        reaction, velocity,
                                        error, NULL);
      }

      if (verbose > 0)
        printf("----------------------------------- FC3D - Projected Gradient On Cylinder (PGoC) - Iteration %i rho = %14.7e \tError = %14.7e\n", iter, rho, error);

      if (error < tolerance) hasNotConverged = 0;
      *info = hasNotConverged;
    }
  }
  else
  {
    
    double tau = dparam[SICONOS_FRICTION_3D_PGOC_LINESEARCH_TAU] ;
    if ((tau <= 0.0 ) || ( tau >= 1.0))
      numerics_error("fc3d_ProjectedGradientOnCylinder", "dparam[SICONOS_FRICTION_3D_PGOC_LINESEARCH_TAU] must in (0,1)");
    

    double mu = dparam[SICONOS_FRICTION_3D_PGOC_LINESEARCH_MU];
    if ((mu <= 0.0 ) || ( mu >= 1.0))
      numerics_error("fc3d_ProjectedGradientOnCylinder", "dparam[SICONOS_FRICTION_3D_PGOC_LINESEARCH_MU] must in (0,1)");

    int ls_iter_max = iparam[SICONOS_FRICTION_3D_PGOC_LINESEARCH_MAXITER];

    double theta = 0.0;
    double theta_old = 0.0;
    double criterion =1.0;

    int ls_iter =0;

    while ((iter < itermax) && (hasNotConverged > 0))
    {
      //rho =  rhoinit;
      ++iter;

      // store the old reaction
      cblas_dcopy(n , reaction , 1 , reaction_k , 1);

      /* q --> velocity_k */
      cblas_dcopy(n , q , 1 , velocity_k, 1);
      /* M r + q --> velocity_k */
      NM_gemv(1.0, M, reaction_k, 1.0, velocity_k);

      /* Armijo line-search */
      /* Compute the new criterion (we use velocity as tmp) */
      /* velocity_k --> velocity */
      cblas_dcopy(n , velocity_k , 1 , velocity, 1);
      /* M r + 2* q --> velocity */
      cblas_daxpy(n, 1.0, q, 1, velocity, 1);
      theta_old = 0.5*cblas_ddot(n, reaction, 1, velocity, 1);

      criterion =1.0;
      ls_iter=0;
      while ((criterion > 0.0) && (ls_iter < ls_iter_max))
      {
        if (verbose > 1 ) printf("ls_iter = %i\t rho = %e ", ls_iter, rho);

        // Compute the new trial reaction
        // r = P_D(r_k - rho v_k)
        cblas_dcopy(n , reaction_k , 1 , reaction , 1);
        cblas_daxpy(n, -1.0*rho, velocity_k, 1, reaction, 1);
        for (contact = 0 ; contact < nc ; ++contact)
          projectionOnCylinder(&reaction[contact * nLocal], options->dWork[contact]);

        //  Compute the new criterion (we use velocity as tmp)
        /* q --> velocity */
        cblas_dcopy(n , q , 1 , velocity, 1);
        /* M r + q --> velocity */
        NM_gemv(1.0, M, reaction, 1.0, velocity);
        /* M r + 2*q --> velocity */
        cblas_daxpy(n, 1.0, q, 1, velocity, 1);
        /* theta = 0.5 R^T M r + r^T q --> velocity */
        theta = 0.5*  cblas_ddot(n, velocity, 1, reaction, 1);

        // Compute direction d _k = r - r_k
        cblas_dcopy(n , reaction , 1 , direction , 1);
        cblas_daxpy(n, -1, reaction_k, 1, direction, 1);

        criterion =    theta - theta_old - mu * cblas_ddot(n, velocity_k, 1, direction, 1);
        if (verbose > 1 ) printf("theta = %e\t criterion = %e\n  ", theta, criterion);
        rho = rho*tau;
        ls_iter++;
        if (rho < rhomin)
          break;
      }
      rho  = rho/tau;
      /* **** Criterium convergence **** */
      fc3d_Tresca_compute_error(problem, reaction , velocity, tolerance, options,  norm_q, &error);

      if (verbose > 0)
        printf("----------------------------------- FC3D - Projected Gradient On Cylinder (PGoC) - Iteration %i rho =%10.5e error = %10.5e < %10.5e\n", iter, rho, error, tolerance);


      /* Try to increase rho is the ls_iter succeeds in one iteration */
      if (ls_iter == 1)
        rho  = rho/tau;

      if (error < tolerance) hasNotConverged = 0;
      *info = hasNotConverged;
    }
  }




  if (verbose > 0)
  printf("----------------------------------- FC3D - Projected Gradient On Cylinder (PGoC)- #Iteration %i Final Residual = %14.7e\n", iter, error);
  dparam[SICONOS_DPARAM_RESIDU] = error;
  iparam[SICONOS_IPARAM_ITER_DONE] = iter;
  if (isVariable)
  {
    free(reaction_k);
    free(direction);
    free(velocity_k);
  }
}


int fc3d_ProjectedGradientOnCylinder_setDefaultSolverOptions(SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the PGoC Solver\n");
  }

  options->solverId = SICONOS_FRICTION_3D_PGoC;

  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 10;
  options->dSize = 10;

  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  options->dWork = NULL;
  solver_options_nullify(options);
  
  options->iparam[SICONOS_IPARAM_MAX_ITER] = 20000;
  
  options->iparam[SICONOS_FRICTION_3D_PGOC_LINESEARCH_MAXITER] =10;
  
  options->dparam[SICONOS_DPARAM_TOL] = 1e-6;
  options->dparam[SICONOS_FRICTION_3D_PGOC_RHO] = -1e-3; /* rho is variable by default */
  options->dparam[SICONOS_FRICTION_3D_PGOC_RHOMIN] = 1e-9; 
  options->dparam[SICONOS_FRICTION_3D_PGOC_LINESEARCH_MU] =0.9;
  options->dparam[SICONOS_FRICTION_3D_PGOC_LINESEARCH_TAU]  = 2.0/3.0;

  options->internalSolvers = NULL;


  return 0;
}
