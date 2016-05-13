/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.

 * Copyright 2016 INRIA.

 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at

 * http://www.apache.org/licenses/LICENSE-2.0

 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
#include "projectionOnCone.h"
#include "fc3d_Solvers.h"
#include "fc3d_compute_error.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SiconosBlas.h"
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"


void fc3d_fixedPointProjection(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options)
{
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;
  /* Number of contacts */
  int nc = problem->numberOfContacts;
  double* q = problem->q;
  NumericsMatrix* M = problem->M;
  double* mu = problem->mu;
  /* Dimension of the problem */
  int n = 3 * nc;
  /* Maximum number of iterations */
  int itermax = iparam[0];
  /* Tolerance */
  double tolerance = dparam[0];

  /*****  Fixed point iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;
  int contact; /* Number of the current row of blocks in M */
  int nLocal = 3;
  dparam[0] = dparam[2]; // set the tolerance for the local solver
  double * velocitytmp = (double *)malloc(n * sizeof(double));


  double rho = 0.0, rho_k =0.0;
  int isVariable = 0;
  double rhomax = 0.0;

  if (dparam[3] < 0.0)
  {
    rho = -dparam[3];
  }
  /* Variable step in fixed*/
  isVariable = 1;
  rhomax = dparam[3];
  rho = rhomax;
  if (verbose > 0)
  {
    printf("----------------------------------- FC3D -  Fixed Point Projection (FPP) - Variable stepsize with starting rho = %14.7e \n", rho);
  }



  double alpha = 1.0;
  double beta = 1.0;

  /* Variable for Line_search */
  int success = 0;
  double error_k;
  int ls_iter = 0;
  int ls_itermax = 10;
  double tau=0.6, L= 0.9, Lmin =0.3, taumin=0.7;



  double a1=0.0, a2=0.0;
  double * reaction_k = 0;
  double * velocity_k = 0;
  double * reactiontmp = 0;
  if (isVariable)
  {
    reaction_k = (double *)malloc(n * sizeof(double));
    velocity_k = (double *)malloc(n * sizeof(double));
    reactiontmp = (double *)malloc(n * sizeof(double));
  }



  // Compute new rho if variable
  if (isVariable)
  {
    while ((iter < itermax) && (hasNotConverged > 0))
    {
      ++iter;

      /* Store the error */
      error_k = error;
      /* store the reaction at the beginning of the iteration */
      cblas_dcopy(n , reaction , 1 , reaction_k, 1);

      /* velocity_k <- q  */
      cblas_dcopy(n , q , 1 , velocity_k, 1);

      /* velocity_k <- q + M * reaction_k  */
      beta = 1.0;
      prodNumericsMatrix(n, n, alpha, M, reaction_k, beta, velocity_k);

      ls_iter = 0 ;
      success =0;

      while (!success && (ls_iter < ls_itermax))
      {

        rho_k = rho * pow(tau,ls_iter);



        /* projection for each contact */
        for (contact = 0 ; contact < nc ; ++contact)
        {
          int pos = contact * nLocal;
          double  normUT = sqrt(velocity_k[pos + 1] * velocity_k[pos + 1] + velocity_k[pos + 2] * velocity_k[pos + 2]);
          /* reaction[pos] = reaction_k[pos] -  rho_k * (velocity_k[pos] + mu[contact] * normUT); */
          /* reaction[pos + 1] = reaction_k[pos+1] - rho_k * velocity_k[pos + 1]; */
          /* reaction[pos + 2] = reaction_k[pos+2] - rho_k * velocity_k[pos + 2]; */

          reaction[pos] -= rho_k * (velocity_k[pos] + mu[contact] * normUT);
          reaction[pos + 1] -= rho_k * velocity_k[pos + 1];
          reaction[pos + 2] -= rho_k * velocity_k[pos + 2];

          projectionOnCone(&reaction[pos], mu[contact]);
        }


        /* velocity <- q + M * reaction  */
        cblas_dcopy(n , q , 1 , velocity, 1);
        beta = 1.0;
        prodNumericsMatrix(n, n, alpha, M, reaction, beta, velocity);

        DEBUG_EXPR_WE( for (int i =0; i< 5 ; i++)
                       {
                         printf("reaction[%i]=%12.8e\t",i,reaction[i]);    printf("velocity[%i]=F[%i]=%12.8e\n",i,i,velocity[i]);
                       }
          );


        /* velocitytmp <- velocity */
        cblas_dcopy(n, velocity, 1, velocitytmp , 1) ;

        for (contact = 0 ; contact < nc ; ++contact)
        {
          int pos = contact * nLocal;
          double  normUT = sqrt(velocitytmp[pos + 1] * velocitytmp[pos + 1]
                                + velocitytmp[pos + 2] * velocitytmp[pos + 2]);
          double  normUT_k = sqrt(velocity_k[pos + 1] * velocity_k[pos + 1] + velocity_k[pos + 2] * velocity_k[pos + 2]);
          velocitytmp[pos] += mu[contact] * (normUT -normUT_k)  ;
        }

        /* velocitytmp <- velocity - velocity_k   */
        cblas_daxpy(n, -1.0, velocity_k , 1, velocitytmp , 1) ;

        a1 = cblas_dnrm2(n, velocitytmp, 1);
        DEBUG_PRINTF("a1 = %12.8e\n", a1);

        /* reactiontmp <- reaction */
        cblas_dcopy(n, reaction, 1, reactiontmp , 1) ;

        /* reactiontmp <- reaction - reaction_k   */
        cblas_daxpy(n, -1.0, reaction_k , 1, reactiontmp , 1) ;

        a2 = cblas_dnrm2(n, reactiontmp, 1) ;
        DEBUG_PRINTF("a2 = %12.8e\n", a2);

        success = (rho_k*a1 <= L * a2)?1:0;
        /* printf("rho_k = %12.8e\t", rho_k); */
        /* printf("a1 = %12.8e\t", a1); */
        /* printf("a2 = %12.8e\t", a2); */
        /* printf("norm reaction = %12.8e\t",cblas_dnrm2(n, reaction, 1) ); */
        /* printf("success = %i\n", success); */

        ls_iter++;
      }

      DEBUG_EXPR_WE( for (int i =0; i< 5 ; i++)
                     {
                       printf("reaction[%i]=%12.8e\t",i,reaction[i]);    printf("velocity[%i]=F[%i]=%12.8e\n",i,i,velocity[i]);
                     };
        );


      /* **** Criterium convergence **** */
      fc3d_compute_error(problem, reaction , velocity, tolerance, options, &error);
      DEBUG_EXPR_WE(
        if ((error < error_k))
        {
          printf("(error < error_k) is satisfied\n");
        }
        );
      /*Update rho*/
      if ((rho_k*a1 < Lmin * a2) && (error < error_k))
      {
        rho =rho_k/taumin;
        DEBUG_PRINTF("We compute a new rho_k = \n", rho_k);

      }
      else
        rho =rho_k;

    if (verbose > 0)
      printf("----------------------------------- FC3D -  Fixed Point Projection (FPP) - Iteration %i rho = %14.7e \tError = %14.7e\n", iter, rho, error);

    if (error < tolerance) hasNotConverged = 0;
    *info = hasNotConverged;
  }

  }



  if (verbose > 0)
    printf("----------------------------------- FC3D - Fixed Point Projection (FPP) - #Iteration %i Final Residual = %14.7e\n", iter, error);
  iparam[7] = iter;
  dparam[0] = tolerance;
  dparam[1] = error;
  free(velocitytmp);
  free(reaction_k);
  free(velocity_k);
  free(reactiontmp);

}


int fc3d_fixedPointProjection_setDefaultSolverOptions(SolverOptions* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the FPP Solver\n");
  }

  /*strcpy(options->solverName,"DSFP");*/
  options->solverId = SICONOS_FRICTION_3D_FPP;
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
  options->dparam[3] = 1.0; /* Default value for rho (line search activated)*/

  options->internalSolvers = NULL;

  return 0;
}
