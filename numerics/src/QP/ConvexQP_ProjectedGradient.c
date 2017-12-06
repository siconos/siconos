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


#include "ConvexQP_Solvers.h"
#include "ConvexQP_computeError.h"
#include "NumericsMatrix.h"
#include "SiconosCompat.h"


#include "SiconosBlas.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ConvexQP_cst.h"
#include "numerics_verbose.h"

//#define VERBOSE_DEBUG
const char* const   SICONOS_CONVEXQP_PG_STR = "CONVEXQP PG";


void convexQP_ProjectedGradient(ConvexQP* problem, double *z, double *w, int* info, SolverOptions* options)
{
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;



  double* q = problem->q;
  NumericsMatrix* M = problem->M;

  NumericsMatrix* A = problem->A;
  if (A)
  {
    numerics_error("ConvexQP_ProjectedGradient", "This solver does not support a specific matrix A different from the identity");
  }

  /* Dimension of the problem */
  int n =  problem->size;

  /* Maximum number of iterations */
  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];
  /* Tolerance */
  double tolerance = dparam[SICONOS_DPARAM_TOL];

  /*****  Projected Gradient iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;

  double rho = 0.0;
  double rhomin = 0.0;
  int isVariable = 0;

  if (dparam[SICONOS_CONVEXQP_PGOC_RHO] > 0.0)
  {
    rho = dparam[SICONOS_CONVEXQP_PGOC_RHO];
  }
  else
  {
    /* Variable step in fixed*/
    isVariable = 1;
    if (verbose > 0)
       printf("Variable step (line search) in Projected Gradient iterations\n");
    rho = -dparam[SICONOS_CONVEXQP_PGOC_RHO];
    rhomin = dparam[SICONOS_CONVEXQP_PGOC_RHOMIN];
  }

  if (rho == 0.0)
    numerics_error("ConvexQP_ProjectedGradient", "dparam[SICONOS_CONVEXQP_PGOC_RHO] must be nonzero");

  double rho_k =rho;
  double * z_tmp= (double *)malloc(n * sizeof(double));
  double * z_k;
  double * direction;
  double * w_k;
  if (isVariable)
  {
    z_k = (double *)malloc(n * sizeof(double));
    direction = (double *)malloc(n * sizeof(double));
    w_k = (double *)malloc(n * sizeof(double));
  }
  double alpha = 1.0;
  double beta = 1.0;

  if (!isVariable)
  {
    double minusrho  = -1.0*rho;

    while ((iter < itermax) && (hasNotConverged > 0))
    {
      ++iter;

      /* q --> w */
      cblas_dcopy(n , q , 1 , w, 1);

      /* M z + q --> w */
      NM_gemv(alpha, M, z, beta, w);

      /* projection for each contact of z - rho * w  */
      cblas_daxpy(n, minusrho, w, 1, z , 1);

      cblas_dcopy(n , z , 1 , z_tmp, 1);
      problem->ProjectionOnC(problem,z_tmp,z);

      /* **** Criterium convergence **** */
      convexQP_computeError(problem, z , w, tolerance, options, &error);

      if (verbose > 0)
        printf("--------------- ConvexQP - Projected Gradient (PG) - Iteration %i rho = %14.7e \tError = %14.7e\n", iter, rho, error);

      if (error < tolerance) hasNotConverged = 0;
      *info = hasNotConverged;
    }
  }
  else
  {

    double tau = dparam[SICONOS_CONVEXQP_PGOC_LINESEARCH_TAU] ;
    if ((tau <= 0.0 ) || ( tau >= 1.0))
      numerics_error("fc3d_ProjectedGradientOnCylinder", "dparam[SICONOS_CONVEXQP_PGOC_LINESEARCH_TAU] must in (0,1)");


    double mu = dparam[SICONOS_CONVEXQP_PGOC_LINESEARCH_MU];
    if ((mu <= 0.0 ) || ( mu >= 1.0))
      numerics_error("fc3d_ProjectedGradientOnCylinder", "dparam[SICONOS_CONVEXQP_PGOC_LINESEARCH_MU] must in (0,1)");

    int ls_iter_max = iparam[SICONOS_CONVEXQP_PGOC_LINESEARCH_MAXITER];

    double theta = 0.0;
    double theta_k = 0.0;
    double correction =0.0;
    double criterion =1.0;

    int ls_iter =0;
    //verbose=2;

    cblas_dcopy(n , z , 1 , z_tmp, 1);
    problem->ProjectionOnC(problem,z_tmp,z);

    /* store the old z */
    cblas_dcopy(n , z , 1 , z_k , 1);

    /* q --> velocity_k */
    cblas_dcopy(n , q , 1 , w_k, 1);
    /* M r + q --> velocity_k */
    NM_gemv(1.0, M, z, 1.0, w_k);

    /* Compute the value fo the cost function (we use w as tmp) */
    /* w_k --> w */
    cblas_dcopy(n , w_k , 1 , w, 1);
    /* M z + 2* q --> w */
    cblas_daxpy(n, 1.0, q, 1, w, 1);
    theta_k = 0.5*cblas_ddot(n, z, 1, w, 1);

    while ((iter < itermax) && (hasNotConverged > 0))
    {
      rho =  rho_k;
      ++iter;

      /* Armijo line-search */
      criterion =1.0;
      ls_iter=0;
      while ((criterion > 0.0) && (ls_iter < ls_iter_max))
      {
        if (verbose > 1 ) printf("ls_iter = %i\t rho = %e\t", ls_iter, rho);

        // Compute the new trial reaction
        // r = P_D(r_k - rho v_k)
        cblas_dcopy(n , z_k , 1 , z , 1);
        cblas_daxpy(n, -1.0*rho, w_k, 1, z, 1);

        cblas_dcopy(n , z , 1 , z_tmp, 1);
        problem->ProjectionOnC(problem,z_tmp,z);

        //  Compute the new value of the cost function (we use velocity as tmp)
        /* q --> w */
        cblas_dcopy(n , q , 1 , w, 1);

        cblas_daxpy(n, 1.0, q, 1, w, 1);
        NM_gemv(1.0, M, z, 1.0, w);

        /* M r + 2*q --> w */
        //NM_gemv(1.0, M, z, 2.0, w);


        /* theta = 0.5 R^T M r + r^T q */
        theta = 0.5*  cblas_ddot(n, w, 1, z, 1);

        // Compute direction d_k = z - z_k
        cblas_dcopy(n , z , 1 , direction , 1);
        cblas_daxpy(n, -1, z_k, 1, direction, 1);

        correction = mu * cblas_ddot(n, w_k, 1, direction, 1);
        criterion =    theta - theta_k - correction;
        if (verbose > 1 ) printf("theta = %e\t theta_k = %e\t  correction = %e\t criterion = %e\n", theta, theta_k, correction, criterion);
        rho = rho*tau;
        ls_iter++;
        if (rho < rhomin)
          break;
      }
      /* store the old z */
      cblas_dcopy(n , z , 1 , z_k , 1);
      /* Compute the value fo the cost function (we use w as tmp) */
      /* q --> velocity_k */
      cblas_dcopy(n , q , 1 , w_k, 1);
      /* M r + q --> velocity_k */
      NM_gemv(1.0, M, z, 1.0, w_k);
      /* w_k --> w */
      cblas_dcopy(n , w_k , 1 , w, 1);
      /* M z + 2* q --> w */
      cblas_daxpy(n, 1.0, q, 1, w, 1);
      theta_k = 0.5*cblas_ddot(n, z, 1, w, 1);
      rho  = rho/tau;
      /* **** Criterium convergence **** */
      convexQP_computeError(problem, z , w, tolerance, options,  &error);

      if (verbose > 0)
        printf("--------------- ConvexQP - Projected Gradient (PG) - Iteration %i rho_k = %10.5e  rho =%10.5e error = %10.5e < %10.5e\n", iter, rho_k, rho, error, tolerance);


      /* Try to increase rho is the ls_iter succeeds in one iteration */
      if (ls_iter == 1){
        if (verbose > 1)
          printf("--------------- ConvexQP - Projected Gradient (PG) -  rho_k is increased from %10.5e to %10.5e \n", rho_k, rho_k/tau);
        rho_k = rho_k/tau;
      }
      if (error < tolerance) hasNotConverged = 0;
      *info = hasNotConverged;
    }
  }


  //verbose=1;

  if (verbose > 0)
  printf("---------------  ConvexQP - Projected Gradient (PG) - #Iteration %i Final Residual = %14.7e\n", iter, error);
  dparam[SICONOS_DPARAM_RESIDU] = error;
  iparam[SICONOS_IPARAM_ITER_DONE] = iter;

  free(z_tmp);
  if (isVariable)
  {
    free(z_k);
    free(direction);
    free(w_k);
  }
}


int convexQP_ProjectedGradient_setDefaultSolverOptions(SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the PGoC Solver\n");
  }

  options->solverId = SICONOS_CONVEXQP_PG;

  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 20;
  options->dSize = 20;

  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  options->dWork = NULL;
  solver_options_nullify(options);

  options->iparam[SICONOS_IPARAM_MAX_ITER] = 20000;

  options->iparam[SICONOS_CONVEXQP_PGOC_LINESEARCH_MAXITER] =20;

  options->dparam[SICONOS_DPARAM_TOL] = 1e-6;
  options->dparam[SICONOS_CONVEXQP_PGOC_RHO] = -1.e-3; /* rho is variable by default */
  options->dparam[SICONOS_CONVEXQP_PGOC_RHOMIN] = 1e-9;
  options->dparam[SICONOS_CONVEXQP_PGOC_LINESEARCH_MU] =0.9;
  options->dparam[SICONOS_CONVEXQP_PGOC_LINESEARCH_TAU]  = 2.0/3.0;

  options->internalSolvers = NULL;


  return 0;
}
