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
#include "projectionOnCone.h"
#include "fc3d_Solvers.h"
#include "fc3d_compute_error.h"
#include "SiconosBlas.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "numerics_verbose.h"

#define VERBOSE_DEBUG

void fc3d_HyperplaneProjection(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options)
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
  /* Maximum number of iterations in Line--search */
  int lsitermax = iparam[1];
  /* Tolerance */
  double tolerance = dparam[0];
  double norm_q = cblas_dnrm2(nc*3 , problem->q , 1);
 




  /*****  Fixed point iterations *****/
  int iter = 0; /* Current iteration number */
  double error = 1.; /* Current error */
  int hasNotConverged = 1;
  int contact; /* Number of the current row of blocks in M */
  int nLocal = 3;
  dparam[0] = dparam[2]; // set the tolerance for the local solver
  double * velocitytmp = (double *)calloc(n, sizeof(double));
  double * reactiontmp = (double *)calloc(n, sizeof(double));
  double * reactiontmp2 = (double *)calloc(n, sizeof(double));
  double * reactiontmp3 = (double *)calloc(n, sizeof(double));

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
 
  /*   double minusrho  = -1.0*rho; */
  while ((iter < itermax) && (hasNotConverged > 0))
  {
    ++iter;

    cblas_dcopy(n , q , 1 , velocitytmp, 1);
    cblas_dcopy(n , reaction , 1 , reactiontmp, 1);

    prodNumericsMatrix(n, n, 1.0, M, reactiontmp, 1.0, velocitytmp);


    // projection for each contact

    double rho = 1;

    for (contact = 0 ; contact < nc ; ++contact)
    {
      int pos = contact * nLocal;
      double  normUT = sqrt(velocitytmp[pos + 1] * velocitytmp[pos + 1] + velocitytmp[pos + 2] * velocitytmp[pos + 2]);
      reactiontmp[pos] -= rho * (velocitytmp[pos] + mu[contact] * normUT);
      reactiontmp[pos + 1] -= rho * velocitytmp[pos + 1];
      reactiontmp[pos + 2] -= rho * velocitytmp[pos + 2];
      projectionOnCone(&reactiontmp[pos], mu[contact]);
    }

    // Armijo line search

    int stopingcriteria = 1;
    int i = -1;
    double alpha ;
    double lhs = NAN;
    double rhs;
    // z_k-y_k
    cblas_dcopy(n , reaction , 1 , reactiontmp3, 1);
    cblas_daxpy(n, -1.0, reactiontmp, 1, reactiontmp3, 1);


    while (stopingcriteria && (i < lsitermax))
    {
      i++ ;
      cblas_dcopy(n , reactiontmp , 1 , reactiontmp2, 1);
      alpha = 1.0 / (pow(2.0, i));
#ifdef VERBOSE_DEBUG
      printf("alpha = %f\n", alpha);
#endif
      cblas_dscal(n , alpha, reactiontmp2, 1);
      alpha  = 1.0 - alpha;

      cblas_daxpy(n, alpha, reaction, 1, reactiontmp2, 1);

      cblas_dcopy(n , q , 1 , velocitytmp, 1);

      prodNumericsMatrix(n, n, 1.0, M, reactiontmp2, 1.0, velocitytmp);



      /* #ifdef VERBOSE_DEBUG */
      /*     for (contact = 0 ; contact < nc ; ++contact) */
      /*     { */
      /*       for(int kk=0; kk<3;kk++) printf("reactiontmp2[%i]=%12.8e\t",contact*nLocal+kk,  reactiontmp2[contact*nLocal+kk]); */
      /*       printf("\n"); */
      /*     } */
      /* #endif   */
      lhs = cblas_ddot(n, velocitytmp, 1, reactiontmp3, 1);
      rhs = cblas_dnrm2(n, reactiontmp3, 1);
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

    double nonorm = cblas_dnrm2(n, velocitytmp, 1);
    double rhoequiv = lhs / (nonorm * nonorm);
#ifdef VERBOSE_DEBUG
    printf("rho equiv = %f\n", rhoequiv);
#endif
    cblas_daxpy(n, -rhoequiv, velocitytmp, 1, reaction  , 1);


    // projection for each contact
    for (contact = 0 ; contact < nc ; ++contact)
    {
      int pos = contact * nLocal;
      projectionOnCone(&reaction[pos], mu[contact]);
    }

    /* **** Criterium convergence **** */
    fc3d_compute_error(problem, reaction , velocity, tolerance, options, norm_q, &error);

    if (options->callback)
    {
      options->callback->collectStatsIteration(options->callback->env, nc * 3, 
                                      reaction, velocity, 
                                      error, NULL);
    }

    if (verbose > 0)
      printf("----------------------------------- FC3D - Hyperplane Projection (HP) - Iteration %i rho = %14.7e \t rhoequiv = %14.7e \tError = %14.7e\n", iter, rho, rhoequiv, error);

    if (error < tolerance) hasNotConverged = 0;
    *info = hasNotConverged;
  }
  if (verbose > 0)
    printf("----------------------------------- FC3D - Hyperplane Projection (HP) - #Iteration %i Final Residual = %14.7e\n", iter, error);
  dparam[0] = tolerance;
  dparam[1] = error;
  iparam[7] = iter;
  free(velocitytmp);
  free(reactiontmp);
  free(reactiontmp2);
  free(reactiontmp3);

}


int fc3d_HyperplaneProjection_setDefaultSolverOptions(SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the HyperplaneProjection Solver\n");
  }

  /*strcpy(options->solverName,"DSFP");*/
  options->solverId = SICONOS_FRICTION_3D_HP;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 8;
  options->dSize = 8;
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  options->dWork = NULL;
  solver_options_nullify(options);
  options->iparam[0] = 2000000;
  options->iparam[1] = 50;
  options->dparam[0] = 1e-3;
  options->dparam[3] = 1.0;
  options->dparam[4] = 0.99;

  options->internalSolvers = NULL;

  return 0;
}
