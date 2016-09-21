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

/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "debug.h"
#include "numerics_verbose.h"
#include <assert.h>

void variationalInequality_HyperplaneProjection(VariationalInequality* problem, double *x, double *w, int* info, SolverOptions* options)
{
  /* /\* int and double parameters *\/ */
  int* iparam = options->iparam;
  double* dparam = options->dparam;
  /* Number of contacts */
  int n = problem->size;
  /* Maximum number of iterations */
  int itermax = iparam[0];
  /* Maximum number of iterations in Line--search */
  int lsitermax = iparam[1];
  assert(lsitermax >0);
  /* Tolerance */
  double tolerance = dparam[0];


  /*****  Fixed point iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;
  dparam[0] = dparam[2]; // set the tolerance for the local solver


  double * xtmp = (double *)malloc(n * sizeof(double));
  double * wtmp = (double *)malloc(n * sizeof(double));
  double * xtmp2 = (double *)malloc(n * sizeof(double));
  double * xtmp3 = (double *)malloc(n * sizeof(double));

  int isVariable = 0;

  double tau = 1.0;
  double sigma = 0.99;

  if (dparam[3] > 0.0)
  {
    tau = dparam[3];
  }
  else
  {
    printf("Hyperplane Projection method. tau <=0  is not well defined\n");
    printf("Hyperplane Projection method. tau is set to 1.0\n");
  }

  if (dparam[4] > 0.0 && dparam[4] < 1.0)
  {
    sigma = dparam[4];
  }
  else
  {
    printf("Hyperplane Projection method. 0<sigma <1  is not well defined\n");
    printf("Hyperplane Projection method. sigma is set to %6.4e\n", sigma);
  }


  isVariable=0;


  if (!isVariable)
  {
    /*   double minusrho  = -1.0*rho; */
    while ((iter < itermax) && (hasNotConverged > 0))
    {
      ++iter;
      /** xtmp <-- x (x_k) */
      cblas_dcopy(n , x , 1 , xtmp, 1);



      
      


      

      /* xtmp (y_k)= P_X(x_k-tau F(x_k)) */
      problem->F(problem, n, xtmp, wtmp);
      cblas_daxpy(n, -tau, wtmp , 1, xtmp , 1) ;
      cblas_dcopy(n , xtmp, 1 , xtmp2, 1);
      problem->ProjectionOnX(problem, xtmp2,xtmp);

      // Armijo line search

      int stopingcriteria = 1;
      int ls_iter = -1;
      double alpha = 1.0;
      double lhs = NAN;
      double rhs;
      // xtmp3 = z_k-y_k
      cblas_dcopy(n , x , 1 , xtmp3, 1);
      cblas_daxpy(n, -1.0, xtmp, 1, xtmp3, 1);
      rhs = cblas_dnrm2(n,xtmp3, 1);
      rhs = sigma / tau * rhs * rhs;
      DEBUG_EXPR_WE(
        for (int i =0; i< n ; i++)
        {
          printf("(y_k) xtmp[%i]=%6.4e\t",i,xtmp[i]);    printf("(x_k-y_k) xtmp3[%i]=%6.4e\n",i,xtmp3[i]);
        }
        );
      while (stopingcriteria && (ls_iter < lsitermax))
      {
        ls_iter++ ;
        /* xtmp2 = alpha * y_k + (1-alpha) x_k */
        alpha = 1.0 / (pow(2.0, ls_iter));
        DEBUG_PRINTF("alpha = %6.4e\n", alpha);
        cblas_dcopy(n ,xtmp , 1 , xtmp2, 1);
        cblas_dscal(n , alpha, xtmp2, 1);
        cblas_daxpy(n, 1.0-alpha, x, 1, xtmp2, 1);



        /* wtmp =  */

        problem->F(problem, n, xtmp2,wtmp);
        DEBUG_EXPR_WE(
          for (int i =0; i< n ; i++)
          {
            printf("(z_k) xtmp2[%i]=%6.4e\n",i,xtmp2[i]);    printf("F(z_k) wtmp[%i]=%6.4e\n",i,wtmp[i]);
          }
          );
        lhs = cblas_ddot(n, wtmp, 1, xtmp3, 1);

        if (lhs >= rhs)  stopingcriteria = 0;

        DEBUG_PRINTF("ls_iter= %i, lsitermax =%i, stopingcriteria  %i\n",ls_iter,lsitermax,stopingcriteria);
        DEBUG_PRINTF("Number of iteration in Armijo line search = %i\t, lhs = %6.4e\t, rhs = %6.4e\t, alpha = %6.4e\t, sigma = %6.4e\t, tau = %6.4e\n", ls_iter, lhs, rhs, alpha, sigma, tau);

      }



      cblas_dcopy(n , x , 1 , xtmp3, 1);
      cblas_daxpy(n, -1.0, xtmp2, 1, xtmp3, 1);
      DEBUG_PRINTF("norm(x-x_k) = %6.4e\n",cblas_dnrm2(n, xtmp3, 1) );
      lhs=cblas_ddot(n, wtmp, 1, xtmp3, 1);

      double nonorm = cblas_dnrm2(n, wtmp, 1);
      double rhoequiv = lhs / (nonorm * nonorm);

      /* rhoequiv =1.0;  */
      DEBUG_PRINTF("nonorm = %6.4e\n", nonorm);
      DEBUG_PRINTF("lhs = %6.4e\n", lhs);
      DEBUG_PRINTF("rho equiv = %6.4e\n", rhoequiv);

      cblas_daxpy(n, -rhoequiv, wtmp, 1, x  , 1);

      cblas_dcopy(n , x, 1 , xtmp, 1);
      problem->ProjectionOnX(problem, xtmp,x);



      /* **** Criterium convergence **** */
      variationalInequality_computeError(problem, x , w, tolerance, options, &error);
      DEBUG_EXPR_WE(
         for (int i =0; i< n ; i++)
         {
           printf("x[%i]=%6.4e\t",i,x[i]);    printf("w[%i]=F[%i]=%6.4e\n",i,i,w[i]);
         }
        );
      if (options->callback)
      {
        options->callback->collectStatsIteration(options->callback->env, n,
                                        x, w,
                                        error, NULL);
      }

      if (verbose > 0)
      {
        printf("----------------------------------- VI - Hyperplane Projection (HP) - Iteration %i tau = %14.7e \t rhoequiv = %14.7e \tError = %14.7e\n", iter, tau, rhoequiv, error);
      }
      if (error < tolerance) hasNotConverged = 0;
      *info = hasNotConverged;
    }
  }



  if (verbose > 0)
  {
    printf("----------------------------------- VI - Hyperplane Projection (HP) - #Iteration %i Final Residual = %14.7e\n", iter, error);
  }
  dparam[0] = tolerance;
  dparam[1] = error;
  iparam[7] = iter;
  free(xtmp);
  free(xtmp2);
  free(xtmp3);
  free(wtmp);

}


int variationalInequality_HyperplaneProjection_setDefaultSolverOptions(SolverOptions* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the HyperplaneProjection Solver\n");
  }

  options->solverId = SICONOS_VI_HP;
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
  options->iparam[0] = 20000;
  options->iparam[1] = 100;

  options->dparam[0] = 1e-3;
  options->dparam[3] = 1.0; /*tau */
  options->dparam[4] = 0.8; /* sigma */

  options->internalSolvers = NULL;

  return 0;
}
