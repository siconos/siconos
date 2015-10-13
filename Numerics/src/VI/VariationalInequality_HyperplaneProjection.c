/* Siconos-Numerics, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
//#include "projectionOnCone.h"
#include "VariationalInequality_Solvers.h"
#include "VariationalInequality_computeError.h"
#include "SiconosBlas.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define VERBOSE_DEBUG
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
    printf("Hyperplane Projection method. rho is set to 1.0\n");

  }
  if (dparam[4] > 0.0 && dparam[4] < 1.0)
  {
    sigma = dparam[4];
  }
  else
  {
    printf("Hyperplane Projection method. 0<sigma <1  is not well defined\n");
    printf("Hyperplane Projection method. sigma is set to 0.99\n");
  }


  isVariable=0;
  if (!isVariable)
  {
    /*   double minusrho  = -1.0*rho; */
    while ((iter < itermax) && (hasNotConverged > 0))
    {
      ++iter;
      /** xtmp <-- x */
      cblas_dcopy(n , x , 1 , xtmp, 1);   

      problem->F(problem, n, xtmp, wtmp);

      double rho = 1;
      cblas_daxpy(n, -rho, wtmp , 1, xtmp , 1) ;
      

      cblas_dcopy(n , xtmp, 1 , xtmp2, 1);   
      problem->ProjectionOnX(problem, xtmp2,xtmp);

      // Armijo line search
      
      int stopingcriteria = 1;
      int i = -1;
      double alpha = 1.0;
      double lhs = NAN;
      double rhs;
      // z_k-y_k
      cblas_dcopy(n , x , 1 , xtmp3, 1);

      cblas_daxpy(n, -1.0, xtmp, 1, xtmp3, 1);

      while (stopingcriteria && (i < lsitermax))
      {
        i++ ;
        cblas_dcopy(n ,xtmp , 1 , xtmp2, 1);

        alpha = 1.0 / (pow(2.0, i));

#ifdef VERBOSE_DEBUG
        printf("alpha = %f\n", alpha);
#endif
        cblas_dscal(n , alpha, xtmp2, 1);
        alpha  = 1.0 - alpha;
        
        cblas_daxpy(n, alpha, x, 1, xtmp2, 1);
        problem->F(problem, n, xtmp2,wtmp);

        lhs = cblas_ddot(n, wtmp, 1, xtmp3, 1);
        rhs = cblas_dnrm2(n,xtmp3, 1);

        rhs = sigma / rho * rhs * rhs;
        if (lhs >= rhs)  stopingcriteria = 0;

#ifdef VERBOSE_DEBUG
        printf("Number of iteration in Armijo line search = %i\n", i);
        printf("lhs = %f\n", lhs);
        printf("rhs = %f\n", rhs);
        printf("alpha = %f\n", alpha);
        printf("sigma = %f\n", sigma);
        printf("rho = %f\n", rho);
#endif
      }
      double nonorm = cblas_dnrm2(n, wtmp, 1);
      double rhoequiv = lhs / (nonorm * nonorm);
#ifdef VERBOSE_DEBUG
      printf("rho equiv = %f\n", rhoequiv);
#endif
      cblas_daxpy(n, -rhoequiv, wtmp, 1, x  , 1);
      
      cblas_dcopy(n , x, 1 , xtmp, 1);        
      problem->ProjectionOnX(problem, xtmp,x);
      
      /* **** Criterium convergence **** */
      variationalInequality_computeError(problem, x , w, tolerance, options, &error);
  
      if (options->callback)
      {
        options->callback->collectStatsIteration(options->callback->env, n,
                                        x, w,
                                        error, NULL);
      }

      if (verbose > 0)
      {
        printf("----------------------------------- VI - Hyperplane Projection (HP) - Iteration %i rho = %14.7e \t rhoequiv = %14.7e \tError = %14.7e\n", iter, rho, rhoequiv, error);
      }
      if (error < tolerance) hasNotConverged = 0;
      *info = hasNotConverged;
    }
  }



  if (verbose > 0)
  {
    printf("----------------------------------- VI - Hyperplane Projection (HP) - #Iteration %i Final Error = %14.7e\n", iter, error);
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
  null_SolverOptions(options);
  for (i = 0; i < 8; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->iparam[0] = 20000;
  options->dparam[0] = 1e-3;
  options->dparam[3] = 1e-3;
  options->dparam[3] = -1.0; // rho is variable by default
  
  options->internalSolvers = NULL;

  return 0;
}
