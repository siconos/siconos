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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "LA.h"
#include "LCP_Solvers.h"

void lcp_newton_FB(LinearComplementarityProblem* problem, double *z, double *w, int *info , SolverOptions* options)
{
  /* matrix M/vector q of the lcp */
  double * M = problem->M->matrix0;

  double * q = problem->q;

  /* size of the LCP */
  int n = problem->size;

  int i, j, iter;
  int m, k;

  int incx, incy;
  double err, a1, b1;
  double alpha, normi;
  int infoDGESV;

  int *ipiv;
  double *beta, *mbeta;
  double *JacPhi, *JacPhi_copy, *Phi;
  int itermax = options->iparam[0];
  double tol = options->dparam[0];


  printf("The Algorithm lcp_newton_FB is not reliable yet, report to siconos.gforge.inria.fr if you need it soon \n");
  return ;


  incx = 1;
  incy = 1;
  /*output*/

  options->iparam[1] = 0;
  options->dparam[1] = 0.0;

  for (i = 0; i < n; i++) z[i] = 0.0;


  // Creation of the gradient of the function H

  JacPhi = (double *)malloc(n * n * sizeof(double));
  JacPhi_copy   = (double *)malloc(n * n * sizeof(double));

  for (j = 0; j < n; j++)
    for (i = 0; i < n; i++) JacPhi[j * n + i] = 0.0;


  // Creation of the RHS Phi,
  Phi = (double *)malloc(n * sizeof(double));
  for (i = 0; i < n; i++) Phi[i] = 0.0;

  // FP- m was not set, I guess m=n ??
  m = n;

  beta = (double *)malloc(m * sizeof(double));
  mbeta = (double *)malloc(m * sizeof(double));

  iter = 0;
  err  = 1.;


  // Newton Iteration
  while ((iter < itermax) && (err > tol))
  {
    ++iter;

    // Construction of the directional derivatives of Phi, JacPhi
    // q --> w
    DCOPY(n , q , incx , w , incy);
    // Mz+q --> w
    a1 = 1.;
    b1 = 1.;
    DGEMV(LA_TRANS , n , n , a1 , M , n , z , incx , b1 , w , incy);
    for (i = 0; i < n; i++) printf("z[%i]=%e", i, z[i]);
    printf("\n");
    for (i = 0; i < n; i++) printf("w[%i]=%e", i, w[i]);
    printf("\n");



    for (i = 0; i < n; i++)
    {
      if ((z[i] == 0) && (w[i] == 0))
      {
        beta[i] = 1.0;
      }
      else
      {
        beta[i] = 0.0;
      }

    }
    for (i = 0; i < n; i++) printf("beta[%i]=%e", i, beta[i]);
    printf("\n");

    // M^T.beta --> mbeta
    a1 = 1.;
    b1 = 0.0;
    DGEMV(LA_NOTRANS , n , n , a1 , M , n , beta , incx , b1 , mbeta  , incy);
    for (i = 0; i < n; i++) printf("mbeta[%i]=%e", i, mbeta[i]);
    printf("\n");



    for (i = 0; i < n; i++)
    {

      if ((z[i] == 0) && (w[i] == 0))
      {
        normi = sqrt(beta[i] * beta[i] + mbeta[i] * mbeta[i]);
        for (j = 0; j < n; j++)
        {
          JacPhi[j * n + i] = (mbeta[i] / normi - 1.0) * M[j * n + i];
        }
        JacPhi[i * n + i] += (beta[i] / normi - 1.0);

      }
      else
      {
        normi = (z[i] * z[i] + w[i] * w[i]);
        printf("normi=%e", normi);
        for (j = 0; j < n; j++)
        {
          JacPhi[j * n + i] = (w[i] / normi - 1.0) * M[j * n + i];

        }
        JacPhi[i * n + i] += (z[i] / normi - 1.0);


      }

    }
    for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++) printf("JacPhi[%i][%i]=%e", i, j, JacPhi[j * n + i]);
      printf("\n");
    }

    // Computation of the value Phi
    for (i = 0; i < n; i++)
    {
      Phi[i] = sqrt(z[i] * z[i] + w[i] * w[i]) - (z[i] + w[i]);

    }




    // Computation of the element of the subgradient.

    DCOPY(n , JacPhi , incx , JacPhi_copy , incy);
    k = 1;
    DGESV(m, k, JacPhi_copy, m, ipiv, beta, m, &infoDGESV);

    if (infoDGESV)
    {
      if (verbose > 0)
      {
        printf("Problem in DGESV\n");
      }
      options->iparam[1] = iter;
      options->dparam[1] = err;
      *info = 2;

      return ;

    }


    // iteration
    alpha = -1.0;
    DAXPY(n , alpha , beta , incx , z , incy);     //  z-beta --> z


    // Construction of the RHS for the next iterate and for the error evaluation
    // q --> w
    DCOPY(n , q , incx , w , incy);
    // Mz+q --> w
    a1 = 1.;
    b1 = 1.;
    DGEMV(LA_TRANS , n , n , a1 , M , n , z , incx , b1 , w , incy);

    for (i = 0; i < n; i++)
    {
      Phi[i] = sqrt(z[i] * z[i] + w[i] * w[i]) - (z[i] + w[i]);

    }


    // Error Evaluation



    err = DNRM2(n , Phi , incx);
    err = 1 / 2 * err * err;

  }

  options->iparam[1] = iter;
  options->dparam[1] = err;

  if (verbose > 0)
  {
    if (err > tol)
    {
      printf(" No convergence of NLGS after %d iterations\n" , iter);
      printf(" The residue is : %g \n", err);
      *info = 1;
    }
    else
    {
      printf(" Convergence of NLGS after %d iterations\n" , iter);
      printf(" The residue is : %g \n", err);
      *info = 0;
    }
  }
  else
  {
    if (err > tol) *info = 1;
    else *info = 0;
  }

  free(JacPhi);
  free(JacPhi_copy);
  free(Phi);
  free(beta);
  free(mbeta);



}
int linearComplementarity_newton_FB_setDefaultSolverOptions(SolverOptions* options)
{
  int i;
  if (verbose > 0)
  {
    printf("Set the Default SolverOptions for the NewtonFB Solver\n");
  }



  options->solverId = SICONOS_LCP_NEWTONFB;
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
  options->iparam[0] = 1000;
  options->dparam[0] = 1e-6;



  return 0;
}

