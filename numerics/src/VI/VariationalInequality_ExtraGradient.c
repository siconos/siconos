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

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "numerics_verbose.h"

/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h" */

#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"
#endif

void variationalInequality_ExtraGradient(VariationalInequality* problem, double *x, double *w, int* info, SolverOptions* options)
{

  DEBUG_BEGIN("variationalInequality_ExtraGradient(VariationalInequality* problem, ...)\n")
  /* /\* int and double parameters *\/ */
  int* iparam = options->iparam;
  double* dparam = options->dparam;
  /* Number of contacts */
  int n = problem->size;
  /* Maximum number of iterations */
  int itermax = iparam[0];
  /* Tolerance */
  double tolerance = dparam[0];

  DEBUG_EXPR(NV_display(x,n););
  DEBUG_EXPR(NV_display(w,n););
  /*****  Fixed point iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;
  dparam[0] = dparam[2]; // set the tolerance for the local solver


  double * xtmp = (double *)calloc(n,sizeof(double));
  double * wtmp = (double *)calloc(n,sizeof(double));

  double rho = 0.0, rho_k =0.0;
  int isVariable = 0;

  if (dparam[3] > 0.0)
  {
    rho = dparam[3];
    if (verbose > 0)
    {
      printf("----------------------------------- VI - Extra Gradient (EG) - Fixed stepsize with  rho = %14.7e \n", rho);
    }
  }
  else
  {
    /* Variable step in iterations*/
    isVariable = 1;
    rho = -dparam[3];
    if (verbose > 0)
    {
      printf("----------------------------------- VI - Extra Gradient (EG) - Variable stepsize with starting rho = %14.7e \n", rho);
    }

  }

  /* Variable for Line_search */
  int success =0;
  double error_k, light_error_sum =0.0;
  int ls_iter = 0;
  int ls_itermax = 10;
  double tau=dparam[4], tauinv=dparam[5], L= dparam[6], Lmin = dparam[7];
  double a1=0.0, a2=0.0;
  double * x_k =0;
  double * w_k =0;

  if (isVariable)
  {
    x_k = (double *)malloc(n * sizeof(double));
    w_k = (double *)malloc(n * sizeof(double));
  }
  /* memcpy(x,x_k,n * sizeof(double)); */
  /* memcpy(w,w_k,n * sizeof(double)); */
  //isVariable=0;
  if (!isVariable)
  {
    /*   double minusrho  = -1.0*rho; */
    while ((iter < itermax) && (hasNotConverged > 0))
    {
      ++iter;

      /* xtmp <- x  */
      cblas_dcopy(n , x , 1 , xtmp, 1);
      
     
      /* wtmp <- F(xtmp) */
      problem->F(problem, n, xtmp,wtmp);
      
      /* xtmp <- xtmp - F(xtmp) */
      cblas_daxpy(n, -1.0, wtmp , 1, xtmp , 1) ;

      /* wtmp <-  ProjectionOnX(xtmp) */
      problem->ProjectionOnX(problem,xtmp,wtmp);

      /* x <- x - wtmp */
      cblas_daxpy(n, -1.0, wtmp , 1, x , 1) ;

      /* x <-  ProjectionOnX(x) */
      cblas_dcopy(n , xtmp , 1 , x, 1);

      problem->ProjectionOnX(problem,xtmp,x);


      /* problem->F(problem,x,w); */
      /* cblas_daxpy(n, -1.0, w , 1, x , 1) ; */
      /* cblas_dcopy(n , x , 1 , xtmp, 1); */
      /* problem->ProjectionOnX(problem,xtmp,x); */

      /* **** Criterium convergence **** */
      if (options->iparam[SICONOS_VI_ERROR_EVALUATION] == SICONOS_VI_ERROR_EVALUATION_FULL )
      {
        variationalInequality_computeError(problem, x , w, tolerance, options, &error);
      }
      else if (options->iparam[SICONOS_VI_ERROR_EVALUATION] == SICONOS_VI_ERROR_EVALUATION_LIGHT )
      {
        cblas_dcopy(n, xtmp, 1,x , 1) ;
        cblas_daxpy(n, -1.0, x_k , 1, xtmp , 1) ;
        light_error_sum = cblas_dnrm2(n,xtmp,1);
        double norm_x= cblas_dnrm2(n,x,1);
        if (fabs(norm_x) > DBL_EPSILON)
          light_error_sum /= norm_x;
        error=light_error_sum;
      }

      if (options->callback)
      {
        options->callback->collectStatsIteration(options->callback->env, n,
                                        x, w,
                                        error, NULL);
      }

      if (verbose > 0)
      {
        printf("----------------------------------- VI - Extra Gradient (EG) - Iteration %i rho = %14.7e \tError = %14.7e\n", iter, rho, error);
      }
      if (error < tolerance) hasNotConverged = 0;
      *info = hasNotConverged;
    }
  }

  if (isVariable)
  {
    if (iparam[1]==0)/* Armijo rule with Khotbotov ratio (default)   */
    {
      while ((iter < itermax) && (hasNotConverged > 0))
      {
        ++iter;


        /* Store the error */
        error_k = error;

        /* x_k <-- x store the x at the beginning of the iteration */
        cblas_dcopy(n , x , 1 , x_k, 1);
        DEBUG_EXPR(NV_display(x_k,n););
        problem->F(problem, n, x, w_k);

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

        /* velocitytmp <- q  */
        /* cblas_dcopy(n , q , 1 , velocitytmp, 1); */
        /* NM_gemv(alpha, M, reaction, beta, velocitytmp); */

        problem->F(problem, n, x,wtmp);

        /* x <- x - rho_k*  wtmp */
        cblas_daxpy(n, -rho_k, wtmp , 1, x , 1) ;

        /* wtmp <-  ProjectionOnX(xtmp) */
        cblas_dcopy(n , x , 1 , xtmp, 1);
        problem->ProjectionOnX(problem,xtmp,x);

        DEBUG_EXPR(NV_display(x,n););
        DEBUG_EXPR(NV_display(w,n););


        /* **** Criterium convergence **** */
        if (options->iparam[SICONOS_VI_ERROR_EVALUATION] == SICONOS_VI_ERROR_EVALUATION_FULL )
        {
          variationalInequality_computeError(problem, x , w, tolerance, options, &error);
        }
        else if (options->iparam[SICONOS_VI_ERROR_EVALUATION] == SICONOS_VI_ERROR_EVALUATION_LIGHT )
        {
          cblas_dcopy(n, xtmp, 1,x , 1) ;
          cblas_daxpy(n, -1.0, x_k , 1, xtmp , 1) ;
          light_error_sum = cblas_dnrm2(n,xtmp,1);
          double norm_x= cblas_dnrm2(n,x,1);
          if (fabs(norm_x) > DBL_EPSILON)
            light_error_sum /= norm_x;
          error=light_error_sum;
        }

        DEBUG_PRINTF("error = %12.8e\t error_k = %12.8e\n",error,error_k);
        /*Update rho*/
        if ((rho_k*a1 < Lmin * a2) && (error < error_k))
        {
          rho =rho_k*tauinv;
        }
        else
          rho =rho_k;


        if (verbose > 0)
        {
          printf("----------------------------------- VI - Extra Gradient (EG) - Iteration %i rho = %14.7e \tError = %14.7e\n", iter, rho, error);
        }
        if (error < tolerance) hasNotConverged = 0;
        *info = hasNotConverged;
      }
    }// end iparam[1]==0

    if (iparam[1]==1) /* Armijo rule with Solodov.Tseng ratio */
    {
      while ((iter < itermax) && (hasNotConverged > 0))
      {
        ++iter;


        /* Store the error */
        error_k = error;

        /* x_k <-- x store the x at the beginning of the iteration */
        cblas_dcopy(n , x , 1 , x_k, 1);

        problem->F(problem, n, x, w_k);

        ls_iter = 0 ;
        success =0;
        rho_k=rho / tau;

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



        /* velocitytmp <- q  */
        /* cblas_dcopy(n , q , 1 , velocitytmp, 1); */
        /* NM_gemv(alpha, M, reaction, beta, velocitytmp); */

        problem->F(problem, n, x,wtmp);

        /* x <- x - rho_k*  wtmp */
        cblas_daxpy(n, -rho_k, wtmp , 1, x , 1) ;

        /* wtmp <-  ProjectionOnX(xtmp) */
        cblas_dcopy(n , x , 1 , xtmp, 1);
        problem->ProjectionOnX(problem,xtmp,x);
        DEBUG_EXPR_WE( for (int i =0; i< 5 ; i++)
                       {
                         printf("x[%i]=%12.8e\t",i,x[i]);    printf("w[%i]=F[%i]=%12.8e\n",i,i,w[i]);
                       }
          );



        /* **** Criterium convergence **** */
        if (options->iparam[SICONOS_VI_ERROR_EVALUATION] == SICONOS_VI_ERROR_EVALUATION_FULL )
        {
          variationalInequality_computeError(problem, x , w, tolerance, options, &error);
        }
        else if (options->iparam[SICONOS_VI_ERROR_EVALUATION] == SICONOS_VI_ERROR_EVALUATION_LIGHT )
        {
          cblas_dcopy(n, xtmp, 1,x , 1) ;
          cblas_daxpy(n, -1.0, x_k , 1, xtmp , 1) ;
          light_error_sum = cblas_dnrm2(n,xtmp,1);
          double norm_x= cblas_dnrm2(n,x,1);
          if (fabs(norm_x) > DBL_EPSILON)
            light_error_sum /= norm_x;
          error=light_error_sum;
        }
        DEBUG_PRINTF("error = %12.8e\t error_k = %12.8e\n",error,error_k);
        /*Update rho*/
        if ((rho_k*a1 < Lmin * a2*a2) && (error < error_k))
        {
          rho =rho_k*tauinv;
        }
        else
          rho =rho_k;


        if (verbose > 0)
        {
          printf("----------------------------------- VI - Extra Gradient (EG) - Iteration %i rho = %14.7e \tError = %14.7e\n", iter, rho, error);
        }
        if (error < tolerance) hasNotConverged = 0;
        *info = hasNotConverged;
      }
    }// end iparam[1]==1



    free(x_k);
    free(w_k);
  }


  if (verbose > 0)
  {
    printf("----------------------------------- VI - Extra Gradient (EG) - #Iteration %i Final Residual = %14.7e\n", iter, error);
  }
  dparam[0] = tolerance;
  dparam[1] = error;
  iparam[7] = iter;
  free(xtmp);
  free(wtmp);
  DEBUG_END("variationalInequality_ExtraGradient(VariationalInequality* problem, ...)\n")

}


int variationalInequality_ExtraGradient_setDefaultSolverOptions(SolverOptions* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the ExtraGradient Solver\n");
  }

  options->solverId = SICONOS_VI_EG;
  options->numberOfInternalSolvers = 0;
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
  options->iparam[0] = 2000000;
  options->dparam[0] = 1e-3;
  options->dparam[3] = 1e-3;
  options->dparam[3] = -1.0; // rho is variable by default
  options->dparam[4] = 2/3.0;  /* tau */
  options->dparam[5] = 3.0/2.0;  /*tauinv */
  options->dparam[6] = 0.9;  /* L */
  options->dparam[7] = 0.3;  /* Lmin */

  options->internalSolvers = NULL;

  return 0;
}
