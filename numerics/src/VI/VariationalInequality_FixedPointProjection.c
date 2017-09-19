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
//#include "projectionOnCone.h"
#include "VariationalInequality_Solvers.h"
#include "VariationalInequality_computeError.h"
#include "SiconosBlas.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
/* #define DEBUG_STDOUT   */
/* #define DEBUG_MESSAGES   */
#include "debug.h"
#include "numerics_verbose.h"

/* #define min(a,b) (a<=b?a:b) */


static
int determine_convergence(double error, double *tolerance, int iter,
                          SolverOptions *options,
                          VariationalInequality* problem,
                          double *z , double *w, double rho)

{
  int hasNotConverged = 1;
  if (error < *tolerance)
  {
    if (verbose > 0)
      printf("----------------------------------- VI - Fixed Point Projection (FPP) - Iteration %i "
             "Residual = %14.7e < %7.3e\n", iter, error, *tolerance);
    double absolute_error =0.0;
    variationalInequality_computeError(problem, z , w, *tolerance, options, &absolute_error);
    if (verbose > 0)
      printf("----------------------------------  VI - Fixed Point Projection (FPP)- Full error criterion =  %e\n", absolute_error);


    hasNotConverged = 0;


    if (options->iparam[SICONOS_VI_IPARAM_ERROR_EVALUATION]==SICONOS_VI_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL)
    {
      if (absolute_error > options->dparam[SICONOS_DPARAM_TOL])
      {
        *tolerance = error/absolute_error*options->dparam[SICONOS_DPARAM_TOL];
        if (verbose > 0)
          printf("----------------------------------  VI - Fixed Point Projection (FPP)- We modify the required incremental precision to reach accuracy to %e\n", *tolerance);
        hasNotConverged = 1;
      }
      else
      {
        if (verbose > 0)
          printf("---------------------------------- VI - Fixed Point Projection (FPP) - The incremental precision is sufficient to reach accuracy to %e\n", *tolerance);
      }
    }
  }
  else
  {
    if (verbose > 0)
      printf("----------------------------------- VI - Fixed Point Projection (FPP) - Iteration %i rho = %14.7e error = %14.7e > %10.5e \n", iter, rho, error, *tolerance);
  }
  return hasNotConverged;
}

static
double compute_error(VariationalInequality* problem,
                     double *z , double *w, double norm_z_z_k,
                     double tolerance,
                     SolverOptions * options)
{

  double error;
  if (options->iparam[SICONOS_VI_IPARAM_ERROR_EVALUATION]==SICONOS_VI_ERROR_EVALUATION_FULL)
    variationalInequality_computeError(problem, z , w, tolerance, options, &error);
  else if (options->iparam[SICONOS_VI_IPARAM_ERROR_EVALUATION]==SICONOS_VI_ERROR_EVALUATION_LIGHT ||
           options->iparam[SICONOS_VI_IPARAM_ERROR_EVALUATION]==SICONOS_VI_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL)
  {
    /* cblas_dcopy(n , x , 1 , xtmp, 1); */
    /* cblas_daxpy(n, -1.0, x_k , 1, xtmp , 1) ; */
    /* error = cblas_dnorm2(n,xtmp,1); */
    double norm_z = cblas_dnrm2(problem->size , z , 1);
    error = norm_z_z_k;
    if (fabs(norm_z) > DBL_EPSILON)
      error /= norm_z;
  }
  else
  {
    return NAN;
    numerics_error("compute_error(VariationalInequality* problem, ...)", "unknown error evaluation strategy");
  }
  return error;
}



void variationalInequality_FixedPointProjection(VariationalInequality* problem, double *x, double *w, int* info, SolverOptions* options)
{

  //verbose=1;
  /* /\* int and double parameters *\/ */
  int* iparam = options->iparam;
  double* dparam = options->dparam;
  /* Number of contacts */
  int n = problem->size;
  /* Maximum number of iterations */
  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];
  /* Tolerance */
  double tolerance = dparam[SICONOS_DPARAM_TOL];

  /*****  Fixed point iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;

  double * xtmp = (double *)calloc(n , sizeof(double));
  double * wtmp = (double *)calloc(n , sizeof(double));

  double rho = 0.0, rho_k =0.0;
  int isVariable = 0;

  if (dparam[SICONOS_VI_EG_DPARAM_RHO] > 0.0)
  {
    rho = dparam[SICONOS_VI_EG_DPARAM_RHO];
    if (verbose > 0)
    {
      printf("----------------------------------- VI - Fixed Point Projection (FPP) - Fixed stepsize with  rho = %14.7e \n", rho);
    }
  }
  else
  {
    /* Variable step in iterations*/
    isVariable = 1;
    rho = -dparam[SICONOS_VI_EG_DPARAM_RHO];
    if (verbose > 0)
    {
      printf("----------------------------------- VI - Fixed Point Projection (FPP) - Variable stepsize with starting rho = %14.7e \n", rho);
    }
  }


  /* Variable for Line_search */
  int success =0;
  double error_k;
  int ls_iter = 0;
  int ls_itermax = 10;
  double tau=dparam[SICONOS_VI_EG_DPARAM_LS_TAU],
    tauinv=dparam[SICONOS_VI_EG_DPARAM_LS_TAUINV],
    L= dparam[SICONOS_VI_EG_DPARAM_LS_L], Lmin = dparam[SICONOS_VI_EG_DPARAM_LS_LMIN];

  DEBUG_PRINTF("tau=%g, tauinv=%g, L= %g, Lmin = %g",dparam[4], dparam[5],  dparam[6], dparam[7] ) ;
  double a1=0.0, a2=0.0;
  double * x_k = NULL;
  double * w_k = NULL;

  if (isVariable)
  {
    x_k = (double *)malloc(n * sizeof(double));
    w_k = (double *)malloc(n * sizeof(double));
  }

  //isVariable=0;
  if (!isVariable)
  {
    /*   double minusrho  = -1.0*rho; */
    while ((iter < itermax) && (hasNotConverged > 0))
    {
      ++iter;

      problem->F(problem,n,x,w);
      cblas_daxpy(n, -1.0, w , 1, x , 1) ;
      cblas_dcopy(n , x , 1 , xtmp, 1);
      problem->ProjectionOnX(problem,xtmp,x);

      /* **** Criterium convergence **** */

      variationalInequality_computeError(problem, x , w, tolerance, options, &error);

      hasNotConverged = determine_convergence(error, &tolerance, iter, options,
                                              problem, x, w, rho);

      if (options->callback)
      {
        options->callback->collectStatsIteration(options->callback->env, n,
                                        x, w,
                                        error, NULL);
      }

    }
  }
  else if (isVariable)
  {
    if (iparam[SICONOS_VI_IPARAM_LINESEARCH_METHOD]==0) /* Armijo rule with Khotbotov ratio (default)   */
    {
      DEBUG_PRINT("Variable step size method with Armijo rule with Khotbotov ratio (default) \n");
      while ((iter < itermax) && (hasNotConverged > 0))
      {
        ++iter;
        /* Store the error */
        error_k = error;

        /* x_k <-- x store the x at the beginning of the iteration */
        cblas_dcopy(n , x , 1 , x_k, 1);
        /* compute w_k =F(x_k) */
        problem->F(problem,n,x,w_k);

        ls_iter = 0 ;
        success =0;
        rho_k=rho / tau;

        while (!success && (ls_iter < ls_itermax))
        {
          /* if (iparam[3] && ls_iter !=0) rho_k = rho_k * tau * min(1.0,a2/(rho_k*a1)); */
          /* else */ rho_k = rho_k * tau ;

          /* x <- x_k  for the std approach*/
          if (iparam[2]==0) cblas_dcopy(n, x_k, 1, x , 1) ;

          /* x <- x - rho_k*  w_k */
          cblas_daxpy(n, -rho_k, w_k , 1, x , 1) ;

          /* xtmp <-  ProjectionOnX(x) */
          problem->ProjectionOnX(problem,x,xtmp);
          problem->F(problem,n,xtmp,w);

          DEBUG_EXPR_WE( for (int i =0; i< 5 ; i++) { printf("xtmp[%i]=%12.8e\t",i,xtmp[i]);
              printf("w[%i]=F[%i]=%12.8e\n",i,i,w[i]);});
          /* velocitytmp <- velocity */
          /* cblas_dcopy(n, w, 1, wtmp , 1) ; */

          /* velocity <- velocity - velocity_k   */
          cblas_daxpy(n, -1.0, w_k , 1, w , 1) ;


          /* a1 =  ||w - w_k|| */
          a1 = cblas_dnrm2(n, w, 1);
          DEBUG_PRINTF("a1 = %12.8e\n", a1);

          /* xtmp <- x */
          cblas_dcopy(n, xtmp, 1,x , 1) ;

          /* xtmp <- x - x_k   */
          cblas_daxpy(n, -1.0, x_k , 1, xtmp , 1) ;

          /* a2 =  || x - x_k || */
          a2 = cblas_dnrm2(n, xtmp, 1) ;
          DEBUG_PRINTF("a2 = %12.8e\n", a2);

          DEBUG_PRINTF("test rho_k*a1 < L * a2 = %e < %e\n", rho_k*a1 , L * a2 ) ;
          success = (rho_k*a1 < L * a2)?1:0;

          /* printf("rho_k = %12.8e\t", rho_k); */
          /* printf("a1 = %12.8e\t", a1); */
          /* printf("a2 = %12.8e\t", a2); */
          /* printf("norm x = %12.8e\t",cblas_dnrm2(n, x, 1) ); */
          /* printf("success = %i\n", success); */

          ls_iter++;
        }

        /* problem->F(problem,x,w); */
        DEBUG_EXPR_WE( for (int i =0; i< 5 ; i++) { printf("x[%i]=%12.8e\t",i,x[i]);
            printf("w[%i]=F[%i]=%12.8e\n",i,i,w[i]);});

        /* **** Criterium convergence **** */
        /* variationalInequality_computeError(problem, x , w, tolerance, options, &error); */

        error = compute_error(problem, x , w, a2,  tolerance, options);

        hasNotConverged = determine_convergence(error, &tolerance, iter, options,
                                                problem, x, w, rho);

        DEBUG_EXPR_WE(
          if ((error < error_k))
          {
            printf("(error < error_k) is satisfied\n");
          };
          );


        /*Update rho*/
        if ((rho_k*a1 < Lmin * a2))
        {
          rho =rho_k*tauinv;
          DEBUG_PRINTF("We compute a new rho_k = %e\n", rho_k);

        }
        else
          rho =rho_k;



      }
    }

    if (iparam[SICONOS_VI_IPARAM_LINESEARCH_METHOD] == 1) /* Armijo rule with Solodov.Tseng ratio */
    {
      DEBUG_PRINT("Variable step size method with Armijo rule with Solodov.Tseng ratio \n");
      while ((iter < itermax) && (hasNotConverged > 0))
      {
        ++iter;


        /* Store the error */
        error_k = error;

        /* x_k <-- x store the x at the beginning of the iteration */
        cblas_dcopy(n , x , 1 , x_k, 1);
        /* compute w_k =F(x_k) */
        problem->F(problem,n,x_k,w_k);

        ls_iter = 0 ;
        success =0;
        rho_k=rho/tau;
        while (!success && (ls_iter < ls_itermax))
        {

          /* if (iparam[3] && ls_iter !=0) rho_k = rho_k * tau * min(1.0,a2*a2/(rho_k*a1)); */
          /* else */ rho_k = rho_k * tau ;

           /* x <- x_k  for the std approach*/
          if (iparam[2]==0) cblas_dcopy(n, x_k, 1, x , 1) ;

          /* x <- x - rho_k*  w_k */
          cblas_daxpy(n, -rho_k, w_k , 1, x , 1) ;

          /* xtmp <-  ProjectionOnX(x) */
          problem->ProjectionOnX(problem,x,xtmp);
          problem->F(problem,n,xtmp,w);

          DEBUG_EXPR_WE( for (int i =0; i< 5 ; i++) { printf("xtmp[%i]=%12.8e\t",i,xtmp[i]);
              printf("w[%i]=F[%i]=%12.8e\n",i,i,w[i]);});
          /* wtmp <- w */
          /* cblas_dcopy(n, w, 1, wtmp , 1) ; */

          /* w <- w - w_k   */
          cblas_daxpy(n, -1.0, w_k , 1, w , 1) ;

          /* xtmp <- x - x_k   */
          cblas_dcopy(n, xtmp, 1,x , 1) ;
          cblas_daxpy(n, -1.0, x_k , 1, xtmp , 1) ;

          /* a1 =  (w - w_k)^T(x - x_k) */
          a1 = cblas_ddot(n, xtmp, 1, w, 1);
          DEBUG_PRINTF("a1 = %12.8e\n", a1);

          /* a2 =  || x - x_k || */
          a2 = cblas_dnrm2(n, xtmp, 1) ;
          DEBUG_PRINTF("a2 = %12.8e\n", a2);

          DEBUG_PRINTF("test rho_k*a1 < L * a2 * a2 = %e < %e\n", rho_k*a1 , L * a2 * a2 ) ;
          success = (rho_k*a1 < L * a2 * a2)?1:0;

          /* printf("rho_k = %12.8e\t", rho_k); */
          /* printf("a1 = %12.8e\t", a1); */
          /* printf("a2 = %12.8e\t", a2); */
          /* printf("norm x = %12.8e\t",cblas_dnrm2(n, x, 1) ); */
          /* printf("success = %i\n", success); */

          ls_iter++;
        }

        /* problem->F(problem,x,w); */
        DEBUG_EXPR_WE( for (int i =0; i< 5 ; i++) { printf("x[%i]=%12.8e\t",i,x[i]);
            printf("w[%i]=F[%i]=%12.8e\n",i,i,w[i]);});

        /* **** Criterium convergence **** */

        error = compute_error(problem, x , w, a2,  tolerance, options);

        hasNotConverged = determine_convergence(error, &tolerance, iter, options,
                                                problem, x, w, rho);
        DEBUG_EXPR_WE(
          if ((error < error_k))
          {
            printf("(error < error_k) is satisfied\n");
          };
          );


        /*Update rho*/
        if ((rho_k*a1 < Lmin * a2*a2))
        {
          rho =rho_k*tauinv;
          DEBUG_PRINTF("We compute a new rho_k = %e \n", rho_k);

        }
        else
          rho = rho_k;



      }
    }

    if (iparam[SICONOS_VI_IPARAM_LINESEARCH_METHOD] == 2) /* Armijo rule with Han.Sun ratio */
    {
      DEBUG_PRINT("Variable step size method with Armijo rule with Han.Sun ratio \n");
      while ((iter < itermax) && (hasNotConverged > 0))
      {
        ++iter;


        /* Store the error */
        error_k = error;

        /* x_k <-- x store the x at the beginning of the iteration */
        cblas_dcopy(n , x , 1 , x_k, 1);
        /* compute w_k =F(x_k) */
        problem->F(problem,n,x_k,w_k);

        ls_iter = 0 ;
        success =0;
        rho_k=rho/tau;
        while (!success && (ls_iter < ls_itermax))
        {
          /* if (iparam[3] && ls_iter !=0) rho_k = rho_k * tau * min(1.0,a1/(rho_k*a2*a2)); */
          /* else */ rho_k = rho_k * tau ;

          /* x <- x_k  for the std approach*/
          if (iparam[2]==0)  cblas_dcopy(n, x_k, 1, x , 1) ;

          /* x <- x - rho_k*  w_k */
          cblas_daxpy(n, -rho_k, w_k , 1, x , 1) ;

          /* xtmp <-  ProjectionOnX(x) */
          problem->ProjectionOnX(problem,x,xtmp);
          problem->F(problem,n,xtmp,w);

          DEBUG_EXPR_WE( for (int i =0; i< 5 ; i++) { printf("xtmp[%i]=%12.8e\t",i,xtmp[i]);
              printf("w[%i]=F[%i]=%12.8e\n",i,i,w[i]);});
          /* wtmp <- w */
          /* cblas_dcopy(n, w, 1, wtmp , 1) ; */

          /* w <- w - w_k   */
          cblas_daxpy(n, -1.0, w_k , 1, w , 1) ;

          /* xtmp <- x */
          cblas_dcopy(n, xtmp, 1,x , 1) ;

          /* xtmp <- x - x_k   */
          cblas_daxpy(n, -1.0, x_k , 1, xtmp , 1) ;

          a1 = cblas_ddot(n, xtmp, 1, w, 1);
          DEBUG_PRINTF("a1 = %12.8e\n", a1);

          a2 = cblas_dnrm2(n, w, 1) ;
          DEBUG_PRINTF("a2 = %12.8e\n", a2);

          DEBUG_PRINTF("test rho_k*a2*a2 < L * a1 = %e < %e\n", rho_k*a2*a2 , L * a1 ) ;
          success = (rho_k*a2*a2 < L * a1)?1:0;

          /* printf("rho_k = %12.8e\t", rho_k); */
          /* printf("a1 = %12.8e\t", a1); */
          /* printf("a2 = %12.8e\t", a2); */
          /* printf("norm x = %12.8e\t",cblas_dnrm2(n, x, 1) ); */
          /* printf("success = %i\n", success); */

          ls_iter++;
        }

        /* problem->F(problem,x,w); */
        DEBUG_EXPR_WE( for (int i =0; i< 5 ; i++) { printf("x[%i]=%12.8e\t",i,x[i]);
            printf("w[%i]=F[%i]=%12.8e\n",i,i,w[i]);});

        /* **** Criterium convergence **** */

        error = compute_error(problem, x , w, a1,  tolerance, options);
        hasNotConverged = determine_convergence(error, &tolerance, iter, options,
                                                problem, x, w, rho);

        DEBUG_EXPR_WE(
          if ((error < error_k))
          {
            printf("(error < error_k) is satisfied\n");
          };
          );


        /*Update rho*/
        if ((rho_k*a2*a2 < Lmin * a1) && (error < error_k))
        {
          rho =rho_k*tauinv;
          DEBUG_PRINTF("We compute a new rho_k = %e \n", rho_k);

        }
        else
          rho = rho_k;



      }
    }
    if (iparam[SICONOS_VI_IPARAM_LINESEARCH_METHOD] > 2)
    {
      fprintf(stderr, "Numerics, VariationalInequality_FixedPointProjection failed. iparam[SICONOS_VI_IPARAM_LINESEARCH_METHOD] > 2 .\n");
      exit(EXIT_FAILURE);
    }
    /* we return the negative value of rho for multiple call to the solver */
    dparam[SICONOS_VI_EG_DPARAM_RHO] = -rho;
  }// end isvariable=1

  *info = hasNotConverged;
  if (isVariable)
  {
    free(x_k);
    free(w_k);
  }

  dparam[SICONOS_DPARAM_TOL] = tolerance;
  dparam[SICONOS_DPARAM_RESIDU] = error;
  iparam[SICONOS_IPARAM_ITER_DONE] = iter;
  free(xtmp);
  free(wtmp);

}


int variationalInequality_FixedPointProjection_setDefaultSolverOptions(SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the FixedPointProjection Solver\n");
  }

  options->solverId = SICONOS_VI_FPP;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 10;
  options->dSize = 10;
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  solver_options_nullify(options);


  options->iparam[SICONOS_IPARAM_MAX_ITER] = 20000;

  options->iparam[SICONOS_VI_IPARAM_LINESEARCH_METHOD] = 0;

  /* options->iparam[SICONOS_VI_IPARAM_ERROR_EVALUATION]=SICONOS_VI_ERROR_EVALUATION_FULL; */
  options->iparam[SICONOS_VI_IPARAM_ERROR_EVALUATION]=SICONOS_VI_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL;
  options->iparam[SICONOS_VI_IPARAM_ERROR_EVALUATION_FREQUENCY]=0;

  options->dparam[SICONOS_DPARAM_TOL] = 1e-3;

  options->dparam[SICONOS_VI_EG_DPARAM_RHO] = -1.0; // rho is variable by default
  options->dparam[SICONOS_VI_EG_DPARAM_LS_TAU] = 2/3.0;  /* tau */
  options->dparam[SICONOS_VI_EG_DPARAM_LS_TAUINV] = 3.0/2.0;  /*tauinv */
  options->dparam[SICONOS_VI_EG_DPARAM_LS_L] = 0.9;  /* L */
  options->dparam[SICONOS_VI_EG_DPARAM_LS_LMIN] = 0.3;  /* Lmin */

  options->internalSolvers = NULL;

  return 0;
}
