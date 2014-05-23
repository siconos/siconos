/* Siconos-Numerics, Copyright INRIA 2005-2014.
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
#include <float.h>

#include "SiconosLapack.h"
#include "Newton_Methods.h"
#include "LineSearch.h"

#define FBLSA_DWORK_KEEP 2
#define FBLSA_NONMONOTONE_LS 3
#define FBLSA_NONMONOTONE_LS_M 4


//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"

void newton_FBLSA(unsigned int n, double *z, double *F, int *info, void* data, SolverOptions* options, functions_FBLSA* functions)
{
  /* size of the CP */
  assert(n>0);

  unsigned int iter;
  unsigned int nn = n*n;

  int incx, incy;
  double theta, preRHS, gamma, sigma, tau, threshold, p, rho;
  double theta_iter = 0.0;
  double err;

  int *ipiv;
  double *workV1, *workV2;
  double *H, *JacThetaF_merit, *F_merit;
  int itermax = options->iparam[0];
  double tol = options->dparam[0];


  incx = 1;
  incy = 1;

  // gamma in (0, 1) or (0, 0.5) ??
  // inconsistency between Facchinei--Pang and "A Theoretical and Numerical
  // Comparison of Some Semismooth Algorithm for Complementarity Problems"
  // The following values are taken from the latter.
  p = 2.1;
  sigma = 0.9; // sigma is the gamma' in VFBLSA (Facchinei--Pang)
  gamma = 1e-4;
  rho = 1e-8;

  /*output*/
  options->iparam[1] = 0;
  options->dparam[1] = 0.0;

  // Maybe there is a better way to initialize
//  for (unsigned int i = 0; i < n; i++) z[i] = 0.0;

  // if we keep the work space across calls
  if (options->iparam[FBLSA_DWORK_KEEP] && (options->dWork == NULL))
  {
    options->dWork = (double *)malloc((3*n + 2*nn) * sizeof(double));
    options->iWork = (int *)malloc(n * sizeof(int));
  }
  if (options->dWork != NULL)
  {
    F_merit = options->dWork;
    workV1 = F_merit + n;
    workV2 = workV1 + n;
    H = workV2 + n;
    JacThetaF_merit = H + nn;

    assert(options->iWork);
    ipiv = options->iWork;
  }
  else
  {
    F_merit = (double *)calloc(n,  sizeof(double));
    H = (double *)calloc(nn, sizeof(double));
    JacThetaF_merit = (double *)calloc(nn, sizeof(double));
    workV1 = (double *)calloc(n, sizeof(double));
    workV2 = (double *)calloc(n, sizeof(double));

    ipiv = (int *)calloc(n, sizeof(int));
  }

  newton_stats stats_iteration;
  if (options->callback)
  {
    stats_iteration.id = NEWTON_STATS_ITERATION;
  }

  linesearch_data ls_data;
  ls_data.compute_F = functions->compute_F;
  ls_data.compute_F_merit = functions->compute_F_merit;
  ls_data.z = z;
  ls_data.zc = workV2;
  ls_data.F = F;
  ls_data.F_merit = F_merit;
  ls_data.desc_dir = workV1;
  ls_data.data = data;

  if (options->iparam[FBLSA_NONMONOTONE_LS] > 0)
  {
    ls_data.nonmonotone = options->iparam[FBLSA_NONMONOTONE_LS];
    ls_data.M = options->iparam[FBLSA_NONMONOTONE_LS_M];
    ls_data.m = 0;
    ls_data.previous_theta = calloc(ls_data.M, sizeof(double));
  }
  else
  {
    ls_data.nonmonotone = 0;
    ls_data.M = 0;
    ls_data.m = 0;
    ls_data.previous_theta = NULL;
  }

  // if error based on the norm of JacThetaF_merit, do something not too stupid
  // here
  JacThetaF_merit[0] = 1e20;

  iter = 0;

  functions->compute_F(data, z, F);
  functions->compute_F_merit(data, z, F, F_merit);

  // Merit Evaluation
  theta = cblas_dnrm2(n , F_merit , incx);
  theta = 0.5 * theta * theta;

  functions->compute_error(data, z, F, JacThetaF_merit, tol, &err);

  // Newton Iteration
  while ((iter < itermax) && (err > tol))
  {
    ++iter;

    functions->compute_F(data, z, F);

    DEBUG_PRINT("z ");
    DEBUG_EXPR_WE(for (unsigned int i = 0; i < n; ++i)
        { DEBUG_PRINTF("% 2.2e ", z[i]) }
        DEBUG_PRINT("\n"));

    DEBUG_PRINT("F ");
    DEBUG_EXPR_WE(for (unsigned int i = 0; i < n; ++i)
        { DEBUG_PRINTF("% 2.2e ", F[i]) }
        DEBUG_PRINT("\n"));


    // Construction of the directional derivatives of F_merit, JacThetaF_merit

    // Construction of H and F_desc
    if (functions->compute_RHS_desc) // different merit function for the descent calc.(usually min)
    {
      functions->compute_H_desc(data, z, F, workV1, workV2, H);
      functions->compute_RHS_desc(data, z, F, F_merit);
    }
    else
    {
      functions->compute_H(data, z, F, workV1, workV2, H);
      functions->compute_F_merit(data, z, F, F_merit);
      cblas_dgemv(CblasColMajor,CblasTrans, n, n, 1.0, H, n, F_merit, incx, 0.0, JacThetaF_merit, incy);
    }

    DEBUG_PRINT("Directional derivative matrix H\n");
    DEBUG_EXPR_WE(for (unsigned int i = 0; i < n; ++i)
        { for(unsigned int j = 0 ; j < n; ++j)
        { DEBUG_PRINTF("% 2.2e ", H[j * n + i]) }
        DEBUG_PRINT("\n")});

    DEBUG_PRINT("F_desc ");
    DEBUG_EXPR_WE(for (unsigned int i = 0; i < n; ++i)
        { DEBUG_PRINTF("% 2.2e ", F_merit[i]) }
        DEBUG_PRINT("\n"));


    // Find direction by solving H * d = -F_desc
    cblas_dcopy(n, F_merit, incx, workV1, incy);
    cblas_dscal(n, -1.0, workV1, incx);
    int infoDGESV = 0;
    DGESV(n, 1, H, n, ipiv, workV1, n, &infoDGESV);

    DEBUG_PRINT("d_k ");
    DEBUG_EXPR_WE(for (unsigned int i = 0; i < n; ++i)
        { DEBUG_PRINTF("% 2.10e ", workV1[i]) }
        DEBUG_PRINT("\n"));

    // Computation of JacThetaF_merit
    // JacThetaF_merit = H^T * F_merit
    if (functions->compute_RHS_desc) // need to compute H and F_merit for the merit
    {
      // /!\ maide! workV1 cannot be used since it contains the descent
      // direction ...
      functions->compute_H(data, z, F, F_merit, workV2, H);
      functions->compute_F_merit(data, z, F, F_merit);
      cblas_dgemv(CblasColMajor, CblasTrans, n, n, 1.0, H, n, F_merit, incx, 0.0, JacThetaF_merit, incy);
    }


    DEBUG_PRINT("JacThetaF_merit ");
    DEBUG_EXPR_WE(for (unsigned int i = 0; i < n; ++i)
        { DEBUG_PRINTF("% 2.2e ", JacThetaF_merit[i]) }
        DEBUG_PRINT("\n"));


    // xhub :: we should be able to continue even if DGESV fails!
    if (infoDGESV)
    {
      if (functions->compute_RHS_desc) // we are safe here
      {
        DEBUG_PRINT("functions->compute_RHS_desc : no min descent direction ! searching for FB descent direction\n");
        cblas_dcopy(n, F_merit, incx, workV1, incy);
        cblas_dscal(n, -1.0, workV1, incx);
        infoDGESV = 0;
        DGESV(n, 1, H, n, ipiv, workV1, n, &infoDGESV);
      }
      else
      {
        if (verbose > 0)
        {
          printf("Problem in DGESV\n");
        }
        options->iparam[1] = iter;
        options->dparam[1] = theta;
        *info = 2;

        goto newton_FBLSA_free;
      }
    }

    if (infoDGESV == 0)
    {
      // workV1 contains the direction d
      cblas_dcopy(n, z, incx, workV2, incy);
      cblas_daxpy(n, 1.0, workV1, incx, workV2, incy);     //  z + d --> z

      // compute new F_merit value and also the new merit
      functions->compute_F(data, workV2, F);
      functions->compute_F_merit(data, workV2, F, F_merit);


      theta_iter = cblas_dnrm2(n, F_merit, incx);
      theta_iter = 0.5 * theta_iter * theta_iter;
    }
    else
    {
      cblas_dcopy(n, JacThetaF_merit, incx, workV1, incy);
      cblas_dscal(n, -1.0, workV1, incx);
    }

    tau = 1.0;
    if ((theta_iter > sigma * theta) || (infoDGESV > 0 && functions->compute_RHS_desc)) // Newton direction not good enough or min part failed
    {
      if (verbose > 1)
        printf("newton_FBLSA :: pure Newton direction not acceptable theta_iter = %g > %g = theta\n", theta_iter, theta);

      // Computations for the linear search
      // preRHS = gamma * <JacThetaF_merit, d>
      // TODO: find a better name for this variable
      preRHS = cblas_ddot(n, JacThetaF_merit, incx, workV1, incy);

      // we should not compute this if min descent search has failed
      threshold = -rho*pow(cblas_dnrm2(n, workV1, incx), p);
      if (preRHS > threshold)
      {
        if (verbose > 1)
          printf("newton_FBLSA :: direction not acceptable %g > %g\n", preRHS, threshold);
        cblas_dcopy(n, JacThetaF_merit, incx, workV1, incy);
        cblas_dscal(n, -1.0, workV1, incx);
        preRHS = cblas_ddot(n, JacThetaF_merit, incx, workV1, incy);
        DEBUG_PRINT("steepest descent ! d_k ");
        DEBUG_EXPR_WE(for (unsigned int i = 0; i < n; ++i)
            { DEBUG_PRINTF("% 2.2e ", workV1[i]) }
            DEBUG_PRINT("\n"));

      }
      preRHS *= gamma;

      // Line search
      tau = linesearch_Armijo2(n, theta, preRHS, &ls_data);
    }

    cblas_daxpy(n , tau, workV1 , incx , z , incy);     //  z + tau*d --> z

    // Construction of the RHS for the next iterate
    functions->compute_F(data, z, F);
    functions->compute_F_merit(data, z, F, F_merit);

    // Merit Evaluation
    theta = cblas_dnrm2(n , F_merit , incx);
    theta = 0.5 * theta * theta;

    // Error Evaluation
    functions->compute_error(data, z, F, JacThetaF_merit, tol, &err);

    if (options->callback)
    {
      stats_iteration.merit_value = theta;
      stats_iteration.alpha = tau;
      stats_iteration.status = 0;
      options->callback->collectStatsIteration(options->callback->env, n, z, F, err, &stats_iteration);
    }

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

newton_FBLSA_free:

  if (options->dWork == NULL)
  {
    free(H);
    free(JacThetaF_merit);
    free(F_merit);
    free(workV1);
    free(workV2);
    free(ipiv);
  }
  if (ls_data.previous_theta)
  {
    free(ls_data.previous_theta);
  }
}


