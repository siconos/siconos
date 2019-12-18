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
#include <assert.h>                              // for assert
#include <math.h>                                // for pow, NAN
#include <stdio.h>                               // for printf, NULL
#include <stdlib.h>                              // for calloc, free, malloc
#include "NumericsFwd.h"                         // for SolverOptions, Varia...
#include "SolverOptions.h"                       // for SolverOptions, solve...
#include "VI_cst.h"                              // for SICONOS_VI_HP
#include "VariationalInequality.h"               // for VariationalInequality
#include "VariationalInequality_Solvers.h"       // for variationalInequalit...
#include "VariationalInequality_computeError.h"  // for variationalInequalit...
#include "debug.h"                               // for DEBUG_PRINTF, DEBUG_...
#include "numerics_verbose.h"                    // for verbose
#include "SiconosBlas.h"                               // for cblas_dcopy, cblas_d...

void variationalInequality_HyperplaneProjection(VariationalInequality* problem, double *x, double *w, int* info, SolverOptions* options)
{
  /* /\* int and double parameters *\/ */
  int* iparam = options->iparam;
  double* dparam = options->dparam;
  /* Number of contacts */
  int n = problem->size;
  /* Maximum number of iterations */
  int itermax = iparam[SICONOS_IPARAM_MAX_ITER];
  /* Maximum number of iterations in Line--search */
  int lsitermax = iparam[SICONOS_VI_IPARAM_LS_MAX_ITER];
  assert(lsitermax >0);
  /* Tolerance */
  double tolerance = dparam[SICONOS_DPARAM_TOL];


  /*****  Fixed point iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;

  double * xtmp = (double *)calloc(n, sizeof(double));
  double * wtmp = (double *)calloc(n, sizeof(double));
  double * xtmp2 = (double *)calloc(n, sizeof(double));
  double * xtmp3 = (double *)calloc(n, sizeof(double));

  int isVariable = 0;

  double tau = 1.0;
  double sigma = 0.99;

  if(dparam[SICONOS_VI_DPARAM_LS_TAU] > 0.0)
  {
    tau = dparam[SICONOS_VI_DPARAM_LS_TAU];
  }
  else
  {
    printf("Hyperplane Projection method. tau <=0  is not well defined\n");
    printf("Hyperplane Projection method. tau is set to 1.0\n");
  }

  if(dparam[SICONOS_VI_DPARAM_SIGMA] > 0.0 && dparam[SICONOS_VI_DPARAM_SIGMA] < 1.0)
  {
    sigma = dparam[SICONOS_VI_DPARAM_SIGMA];
  }
  else
  {
    printf("Hyperplane Projection method. 0<sigma <1  is not well defined\n");
    printf("Hyperplane Projection method. sigma is set to %6.4e\n", sigma);
  }


  isVariable=0;


  if(!isVariable)
  {
    /*   double minusrho  = -1.0*rho; */
    while((iter < itermax) && (hasNotConverged > 0))
    {
      ++iter;
      /** xtmp <-- x (x_k) */
      cblas_dcopy(n, x, 1, xtmp, 1);

      /* xtmp (y_k)= P_X(x_k-tau F(x_k)) */
      problem->F(problem, n, xtmp, wtmp);
      cblas_daxpy(n, -tau, wtmp, 1, xtmp, 1) ;
      cblas_dcopy(n, xtmp, 1, xtmp2, 1);
      problem->ProjectionOnX(problem, xtmp2,xtmp);

      // Armijo line search

      int stopingcriteria = 1;
      int ls_iter = -1;
      double alpha = 1.0;
      double lhs = NAN;
      double rhs;
      // xtmp3 = z_k-y_k
      cblas_dcopy(n, x, 1, xtmp3, 1);
      cblas_daxpy(n, -1.0, xtmp, 1, xtmp3, 1);
      rhs = cblas_dnrm2(n,xtmp3, 1);
      rhs = sigma / tau * rhs * rhs;
      DEBUG_EXPR_WE(
        for(int i =0; i< n ; i++)
    {
      printf("(y_k) xtmp[%i]=%6.4e\t",i,xtmp[i]);
        printf("(x_k-y_k) xtmp3[%i]=%6.4e\n",i,xtmp3[i]);
      }
      );
      while(stopingcriteria && (ls_iter < lsitermax))
      {
        ls_iter++ ;
        /* xtmp2 = alpha * y_k + (1-alpha) x_k */
        alpha = 1.0 / (pow(2.0, ls_iter));
        DEBUG_PRINTF("alpha = %6.4e\n", alpha);
        cblas_dcopy(n,xtmp, 1, xtmp2, 1);
        cblas_dscal(n, alpha, xtmp2, 1);
        cblas_daxpy(n, 1.0-alpha, x, 1, xtmp2, 1);



        /* wtmp =  */

        problem->F(problem, n, xtmp2,wtmp);
        DEBUG_EXPR_WE(
          for(int i =0; i< n ; i++)
      {
        printf("(z_k) xtmp2[%i]=%6.4e\n",i,xtmp2[i]);
          printf("F(z_k) wtmp[%i]=%6.4e\n",i,wtmp[i]);
        }
        );
        lhs = cblas_ddot(n, wtmp, 1, xtmp3, 1);

        if(lhs >= rhs)  stopingcriteria = 0;

        DEBUG_PRINTF("ls_iter= %i, lsitermax =%i, stopingcriteria  %i\n",ls_iter,lsitermax,stopingcriteria);
        DEBUG_PRINTF("Number of iteration in Armijo line search = %i\t, lhs = %6.4e\t, rhs = %6.4e\t, alpha = %6.4e\t, sigma = %6.4e\t, tau = %6.4e\n", ls_iter, lhs, rhs, alpha, sigma, tau);

      }



      cblas_dcopy(n, x, 1, xtmp3, 1);
      cblas_daxpy(n, -1.0, xtmp2, 1, xtmp3, 1);
      DEBUG_PRINTF("norm(x-x_k) = %6.4e\n",cblas_dnrm2(n, xtmp3, 1));
      lhs=cblas_ddot(n, wtmp, 1, xtmp3, 1);

      double nonorm = cblas_dnrm2(n, wtmp, 1);
      double rhoequiv = lhs / (nonorm * nonorm);

      /* rhoequiv =1.0;  */
      DEBUG_PRINTF("nonorm = %6.4e\n", nonorm);
      DEBUG_PRINTF("lhs = %6.4e\n", lhs);
      DEBUG_PRINTF("rho equiv = %6.4e\n", rhoequiv);

      cblas_daxpy(n, -rhoequiv, wtmp, 1, x, 1);

      cblas_dcopy(n, x, 1, xtmp, 1);
      problem->ProjectionOnX(problem, xtmp,x);



      /* **** Criterium convergence **** */
      variationalInequality_computeError(problem, x, w, tolerance, options, &error);
      DEBUG_EXPR_WE(
        for(int i =0; i< n ; i++)
    {
      printf("x[%i]=%6.4e\t",i,x[i]);
        printf("w[%i]=F[%i]=%6.4e\n",i,i,w[i]);
      }
      );
      if(options->callback)
      {
        options->callback->collectStatsIteration(options->callback->env, n,
            x, w,
            error, NULL);
      }

      if(verbose > 0)
      {
        printf("--------------- VI - Hyperplane Projection (HP) - Iteration %i tau = %14.7e \t rhoequiv = %14.7e \tError = %14.7e\n", iter, tau, rhoequiv, error);
      }
      if(error < tolerance) hasNotConverged = 0;
      *info = hasNotConverged;
    }
  }



  if(verbose > 0)
  {
    printf("--------------- VI - Hyperplane Projection (HP) - #Iteration %i Final Residual = %14.7e\n", iter, error);
  }
  dparam[SICONOS_DPARAM_RESIDU] = error;
  iparam[SICONOS_IPARAM_ITER_DONE] = iter;
  free(xtmp);
  free(xtmp2);
  free(xtmp3);
  free(wtmp);

}


void variationalInequality_HyperplaneProjection_set_default(SolverOptions* options)
{
  options->iparam[SICONOS_VI_IPARAM_LS_MAX_ITER] = 100;
  options->dparam[SICONOS_VI_DPARAM_LS_TAU] = 1.0; /*tau */
  options->dparam[SICONOS_VI_DPARAM_SIGMA] = 0.8; /* sigma */
}
