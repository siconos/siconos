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
#include "projectionOnCone.h"
#include "FrictionContact3D_Solvers.h"
#include "FrictionContact3D_compute_error.h"
#include "LA.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void frictionContact3D_DeSaxceFixedPoint(FrictionContactProblem* problem, double *reaction, double *velocity, int* info, SolverOptions* options)
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

  double rho = 0.0;
  int isVariable = 0;
  double rhomax = 0.0;
  if (dparam[3] > 0.0)
  {
    rho = dparam[3];
  }
  else
  {
    /* Variable step in fixed*/
    isVariable = 1;
    printf("Variable step in fixed point iterations\n");
    rhomax = -dparam[3];
    rho = rhomax;
  }
  double * work1tmp = 0;
  double * work2tmp = 0;
  double * direction = 0;
  if (isVariable)
  {
    work1tmp = (double *)malloc(n * sizeof(double));
    work2tmp = (double *)malloc(n * sizeof(double));
    direction = (double *)malloc(n * sizeof(double));
  }
  double alpha = 1.0;
  double beta = 1.0;

  /*   double minusrho  = -1.0*rho; */
  while ((iter < itermax) && (hasNotConverged > 0))
  {
    ++iter;

    DCOPY(n , q , 1 , velocitytmp, 1);
    if (isVariable)   DCOPY(n , reaction , 1 , work1tmp, 1);
    beta = 1.0;
    prodNumericsMatrix(n, n, alpha, M, reaction, beta, velocitytmp);
    // projection for each contact
    for (contact = 0 ; contact < nc ; ++contact)
    {
      int pos = contact * nLocal;
      double  normUT = sqrt(velocitytmp[pos + 1] * velocitytmp[pos + 1] + velocitytmp[pos + 2] * velocitytmp[pos + 2]);
      reaction[pos] -= rho * (velocitytmp[pos] + mu[contact] * normUT);
      reaction[pos + 1] -= rho * velocitytmp[pos + 1];
      reaction[pos + 2] -= rho * velocitytmp[pos + 2];
      projectionOnCone(&reaction[pos], mu[contact]);
    }

    // Compute new rho if variable
    if (isVariable)
    {
      DCOPY(n , work1tmp , 1 , direction , 1);
      DSCAL(n, -1.0, direction, 1);
      DAXPY(n, 1.0, reaction, 1, direction , 1) ;  // warning compyte -d and not d
      beta = 0.0;
      double alpha1 = DNRM2(n, direction, 1);
      /*    prodNumericsMatrix(n,n, alpha, M, direction, beta,work1tmp ); */
      /*    double alpha2=DDOT(n,direction,1,work1tmp,1); */
      prodNumericsMatrix(n, n, alpha, M, work2tmp, beta, work2tmp);
      double alpha3 = DDOT(n, direction, 1, work2tmp, 1);

      /*    if (alpha2 < 1e-7*alpha1) */
      /*        rho= rhomax; */
      /*    else */
      /*        rho = alpha1/alpha2; */

      if (alpha3 < 1e-7 * alpha1)
        rho = rhomax;
      else
        rho = alpha1 / alpha3;
      /*        printf("alpha1= %14.7e \t, alpha2=%14.7e\t  alpha3=%14.7e\t,rho = %14.7e\n",alpha1,alpha2,alpha3,rho); */
    }


    /* **** Criterium convergence **** */
    FrictionContact3D_compute_error(problem, reaction , velocity, tolerance, options, &error);

    if (verbose > 0)
      printf("----------------------------------- FC3D - DeSaxce Fixed Point (DSFP) - Iteration %i rho = %14.7e \tError = %14.7e\n", iter, rho, error);

    if (error < tolerance) hasNotConverged = 0;
    *info = hasNotConverged;
  }
  printf("----------------------------------- FC3D - DeSaxce Fixed point (DSFP) - #Iteration %i Final Error = %14.7e\n", iter, error);
  dparam[0] = tolerance;
  dparam[1] = error;
  free(velocitytmp);
  if (isVariable)
  {
    free(work1tmp);
    free(work2tmp);
    free(direction);
  }
}


int frictionContact3D_DeSaxceFixedPoint_setDefaultSolverOptions(SolverOptions* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the DSFP Solver\n");
  }

  /*strcpy(options->solverName,"DSFP");*/
  options->solverId = SICONOS_FRICTION_3D_DSFP;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 5;
  options->dSize = 5;
  options->iparam = (int *)malloc(options->iSize * sizeof(int));
  options->dparam = (double *)malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  options->iWork = NULL;
  for (i = 0; i < 5; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->iparam[0] = 20000;
  options->dparam[0] = 1e-3;
  options->dparam[3] = 1e-3;

  options->internalSolvers = NULL;

  return 0;
}
