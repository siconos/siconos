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
#include "LCP_Solvers.h"
#include "SiconosLapack.h"

#define DEBUG_STDOUT
#define DEBUG_MESSAGES
#include "debug.h"

void lcp_newton_FB(LinearComplementarityProblem* problem, double *z, double *w, int *info , SolverOptions* options)
{
  /* matrix M/vector q of the lcp */
  assert(problem);
  assert(problem->M);
  double * M = problem->M->matrix0;
  double * q = problem->q;

  assert(M);
  assert(q);
  /* size of the LCP */
  unsigned int n = problem->size;
  assert(n>0);

  unsigned int iter, done;
  unsigned int nn = n*n;

  int incx, incy;
  double theta, preRHS, gamma, sigma, tau, threshold, p, rho;
  double theta_iter, normi, err;
  int infoDGESV;

  int *ipiv;
  double *workV1, *workV2;
  double *H, *JacThetaFB, *F_FB;
  int itermax = options->iparam[0];
  double tol = options->dparam[0];


  incx = 1;
  incy = 1;

  // gamma in (0, 1) or (0, 0.5) ??
  // inconsistency between Facchinei - Pang and A Theoretical and Numerical
  // Comparison of Some Semismooth Algorithm for Complementarity Problems
  // The following values are taken from the latter paper.
  p = 2.1;
  sigma = 0.9; // sigma is the gamma' in VFBLSA (Facchinei - Pang)
  gamma = 1e-4;
  rho = 1e-8;

  /*output*/
  options->iparam[1] = 0;
  options->dparam[1] = 0.0;

  // Maybe there is a better way to initialize
  for (unsigned int i = 0; i < n; i++) z[i] = 0.0;


  F_FB = (double *)malloc(n * sizeof(double));
  H = (double *)malloc(nn * sizeof(double));
  JacThetaFB = (double *)malloc(nn * sizeof(double));

  ipiv = (int *)malloc(n * sizeof(int));

  workV1 = (double *)malloc(n * sizeof(double));
  workV2 = (double *)malloc(n * sizeof(double));



  iter = 0;

  // Construction of the RHS for the next iterate and for the thetaor evaluation
  // q --> w
  cblas_dcopy(n , q , incx , w , incy);
  // Mz+q --> w
  cblas_dgemv(CblasColMajor,CblasNoTrans, n , n , 1.0 , M , n , z , incx , 1.0 , w , incy);

  for (unsigned int i = 0; i < n; i++)
  {
    F_FB[i] = sqrt(z[i] * z[i] + w[i] * w[i]) - (z[i] + w[i]);
  }


  // Merit Evaluation
  theta = cblas_dnrm2(n , F_FB , incx);
  theta = 0.5 * theta * theta;

  lcp_compute_error(problem, z, w, tol, &err);

  // Newton Iteration
  while ((iter < itermax) && (err > tol))
  {
    ++iter;

    // Computation of the new value w = F(z) = Mz + q
    // q --> w
    cblas_dcopy(n , q , incx , w , incy);
    // Mz+q --> w
    cblas_dgemv(CblasColMajor,CblasNoTrans, n , n , 1.0, M , n , z , incx , 1.0, w , incy);
    DEBUG_PRINT("z ");
    DEBUG_EXPR_WE(for (unsigned int i = 0; i < n; ++i)
      { DEBUG_PRINTF("% 2.2e ", z[i]) }
      DEBUG_PRINT("\n"));

    DEBUG_PRINT("w ");
    DEBUG_EXPR_WE(for (unsigned int i = 0; i < n; ++i)
      { DEBUG_PRINTF("% 2.2e ", w[i]) }
      DEBUG_PRINT("\n"));


    // Construction of the directional derivatives of F_FB, JacThetaFB

    // constructing the set beta
    for (unsigned int i = 0; i < n; ++i)
    {
      if ((z[i] == 0.0) && (w[i] == 0.0))
      {
        workV1[i] = 1.0;
      }
      else
      {
        workV1[i] = 0.0;
      }
    }

    // Needed for H
    // M^T.workV1 --> workV2
    cblas_dgemv(CblasColMajor,CblasTrans, n , n , 1.0 , M , n , workV1 , incx , 0.0 , workV2  , incy);

    // Construction of H
    for (unsigned int i = 0; i < n; ++i)
    {
      if (workV1[i] != 0.0) // i in beta
      {
        normi = sqrt(workV1[i] * workV1[i] + workV2[i] * workV2[i]);
        for (unsigned int j = 0; j < n; j++)
        {
          H[j * n + i] = (workV2[i] / normi - 1.0) * M[j * n + i];
        }
        H[i * n + i] += (workV1[i] / normi - 1.0);

      }
      else // i not in beta
      {
        normi = sqrt(z[i] * z[i] + w[i] * w[i]);
        for (unsigned int j = 0; j < n; j++)
        {
          H[j * n + i] = (w[i] / normi - 1.0) * M[j * n + i];

        }
        H[i * n + i] += (z[i] / normi - 1.0);
      }

    }
    DEBUG_PRINT("Directional derivative matrix H\n");
    DEBUG_EXPR_WE(for (unsigned int i = 0; i < n; ++i)
      { for(unsigned int j = 0 ; j < n; ++j)
      { DEBUG_PRINTF("% 2.2e ", H[j * n + i]) }
      DEBUG_PRINT("\n")});

    // Computation of the value F_FB
    for (unsigned int i = 0; i < n; i++)
    {
      F_FB[i] = sqrt(z[i] * z[i] + w[i] * w[i]) - (z[i] + w[i]);
    }
    DEBUG_PRINT("F_FB ");
    DEBUG_EXPR_WE(for (unsigned int i = 0; i < n; ++i)
      { DEBUG_PRINTF("% 2.2e ", F_FB[i]) }
      DEBUG_PRINT("\n"));


    // Computation of JacThetaFB
    // JacThetaFB = H^T * F_FB
    cblas_dgemv(CblasColMajor,CblasTrans, n, n, 1.0, H, n, F_FB, incx, 0.0, JacThetaFB, incy);

    // Find direction by solving H * d = -F_FB
    cblas_dcopy(n, F_FB, incx, workV1, incy);
    cblas_dscal(n, -1.0, workV1, incx);
    DGESV(n, 1, H, n, ipiv, workV1, n, &infoDGESV);

    // xhub :: we should be able to continue even if DGESV fails!
    if (infoDGESV)
    {
      if (verbose > 0)
      {
        printf("Problem in DGESV\n");
      }
      options->iparam[1] = iter;
      options->dparam[1] = theta;
      *info = 2;

      goto lcp_newton_FB_free;
    }

   // workV1 contains the direction d
   cblas_dcopy(n , q , incx , w , incy);
   cblas_dcopy(n, z, incx, workV2, incy);
   cblas_daxpy(n, 1.0, workV1, incx, workV2, incy);     //  z + d --> z

   // compute new F_FB value and also the new merit
   cblas_dgemv(CblasColMajor,CblasNoTrans, n , n , 1.0 , M , n , workV2 , incx , 1.0 , w , incy);
   for (unsigned int i = 0; i < n; ++i)
   {
     F_FB[i] = sqrt(workV2[i] * workV2[i] + w[i] * w[i]) - (workV2[i] + w[i]);
   }

   theta_iter = cblas_dnrm2(n, F_FB, incx);
   theta_iter = 0.5 * theta_iter * theta_iter;

   tau = 1.0;
   if (theta_iter > sigma * theta) // Newton direction not good enough
   {
     if (verbose > 1)
       printf("lcp_newton_FB :: Newton direction not acceptable %g %g\n", theta_iter, theta);

    // Computations for the linear search
    // preRHS = gamma * <JacThetaFB, d>
    // TODO: find a better name for this variable
    preRHS = cblas_ddot(n, JacThetaFB, incx, workV1, incy);

    threshold = -rho*pow(cblas_dnrm2(n, workV1, incx), p);
    if (preRHS > threshold)
    {
      if (verbose > 1)
        printf("lcp_newton_FB :: direction not acceptable %g %g\n", preRHS, threshold);
      cblas_dcopy(n, JacThetaFB, incx, workV1, incy);
      cblas_dscal(n, -1.0, workV1, incx);
      preRHS = cblas_ddot(n, JacThetaFB, incx, workV1, incy);
    }
    preRHS *= gamma;
    done = 0;
    tau = 2.0; // XXX xhub :: this has to be checked

    // Line search
    // xhub :: if we reuse of tau, we may be hit hard by numerical issues
    // TODO: at the end reset tau to a good value of 2**-i
    do
    {
      // workV1 contains the direction d
      cblas_dcopy(n , q , incx , w , incy);
      cblas_dcopy(n, z, incx, workV2, incy);
      cblas_daxpy(n, tau, workV1, incx, workV2, incy);     //  z + tau*d --> z

      // compute new F_FB
      cblas_dgemv(CblasColMajor,CblasNoTrans, n , n , 1.0 , M , n , workV2 , incx , 1.0 , w , incy);
      for (unsigned int i = 0; i < n; ++i)
      {
        F_FB[i] = sqrt(workV2[i] * workV2[i] + w[i] * w[i]) - (workV2[i] + w[i]);
      }

      theta_iter = cblas_dnrm2(n, F_FB, incx);
      theta_iter = 0.5 * theta_iter * theta_iter;

      if (theta_iter <= theta + tau*preRHS)
      {
        done = 1;
        if (verbose > 1)
          printf("lcp_newton_FB :: tau %g\n", tau);
      }
      else
      {
        // tau too large, decrease it
        tau /= 2.0;
      }

// xhub :: tentative reuse of the previous tau
//       if (theta_iter <= theta + tau*preRHS)
//       {
//         if (alreadyFeasible == 0) // it is the first feasible tau => win
//           done = 1;
//         else if (tau < 1.9)
//           tau *= 2.0;
//         alreadyFeasible = 1;
//       }
//       else
//       {
//         // tau too large, decrease it
//         tau /= 2.0;
//         if (alreadyFeasible == 1) // first time tau is too large
//           done = 1;
//         alreadyFeasible = 0;
//       }
// 
    }
    while(!done);

   }

   // update
    cblas_daxpy(n , tau, workV1 , incx , z , incy);     //  z + tau*d --> z

    // Construction of the RHS for the next iterate
    // q --> w
    cblas_dcopy(n , q , incx , w , incy);
    // Mz+q --> w
    cblas_dgemv(CblasColMajor,CblasNoTrans, n , n , 1.0 , M , n , z , incx , 1.0 , w , incy);

    for (unsigned int i = 0; i < n; i++)
    {
      F_FB[i] = sqrt(z[i] * z[i] + w[i] * w[i]) - (z[i] + w[i]);
    }

    // Merit Evaluation
    theta = cblas_dnrm2(n , F_FB , incx);
    theta = 0.5 * theta * theta;

    // Error Evaluation
    lcp_compute_error(problem, z, w, tol, &err);
  }

  options->iparam[1] = iter;
  options->dparam[1] = theta;

  if (verbose > 0)
  {
    if (theta > tol)
    {
      printf(" No convergence of Newton FB  after %d iterations\n" , iter);
      printf(" The residue is : %g \n", theta);
      *info = 1;
    }
    else
    {
      printf(" Convergence of Newton FB after %d iterations\n" , iter);
      printf(" The residue is : %g \n", theta);
      *info = 0;
    }
  }
  else
  {
    if (theta > tol) *info = 1;
    else *info = 0;
  }

lcp_newton_FB_free:

  free(H);
  free(JacThetaFB);
  free(F_FB);
  free(workV1);
  free(workV2);
  free(ipiv);
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
  options->iWork = NULL;   options->callback = NULL; options->numericsOptions = NULL;
  for (i = 0; i < 5; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->iparam[0] = 1000;
  options->dparam[0] = 1e-6;



  return 0;
}

