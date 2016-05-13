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

#include "Newton_Methods.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "SiconosLapack.h"
#include "ArmijoSearch.h"

#include "lcp_cst.h"
#include "NCP_cst.h"
#include "MCP_cst.h"
#include "VI_cst.h"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"


void newton_LSA(unsigned n, double *z, double *F, int *info, void* data, SolverOptions* options, functions_LSA* functions)
{
  /* size of the problem */
  assert(n>0);

  /* Checking the consistency of the functions_LSA struct */

  assert(functions->compute_F && "functions_LSA lacks compute_F");
  assert(functions->compute_F_merit && "functions_LSA lacks compute_F_merit");
  assert(functions->compute_error && "functions_LSA lacks compute_error");
  assert(((functions->compute_descent_direction) || \
        (functions->compute_RHS_desc && functions->compute_H_desc) || \
        functions->compute_H) && "functions_LSA lacks a way to compute a descente direction");
  assert(((!functions->compute_RHS_desc || !functions->compute_descent_direction) || \
        (functions->compute_JacTheta_merit || functions->compute_H)) && \
      "functions_LSA lacks a way to compute JacTheta_merit");


  unsigned int iter;


  int incx, incy;
  double theta, preRHS, tau, threshold;
  double theta_iter = 0.0;
  double err;

  double *workV1, *workV2;
  double *JacThetaF_merit, *F_merit;
  unsigned int itermax = options->iparam[0];
  double tol = options->dparam[0];

  incx = 1;
  incy = 1;

  /*output*/
  options->iparam[1] = 0;
  options->dparam[1] = 0.0;

  // Maybe there is a better way to initialize
//  for (unsigned int i = 0; i < n; i++) z[i] = 0.0;

  // if we keep the work space across calls
  if (options->iparam[SICONOS_IPARAM_PREALLOC] && (options->dWork == NULL))
  {
    options->dWork = (double *)calloc(4*n, sizeof(double));
    options->iWork = (int *)calloc(n, sizeof(int));
  }

  if (options->dWork)
  {
    F_merit = options->dWork;
    workV1 = F_merit + n;
    workV2 = workV1 + n;
    JacThetaF_merit = workV2 + n;
  }
  else
  {
    F_merit = (double *)calloc(n,  sizeof(double));
    JacThetaF_merit = (double *)calloc(n, sizeof(double));
    workV1 = (double *)calloc(n, sizeof(double));
    workV2 = (double *)calloc(n, sizeof(double));
  }

  newton_stats stats_iteration;
  if (options->callback)
  {
    stats_iteration.id = NEWTON_STATS_ITERATION;
  }

  assert(options->solverData);
  newton_LSA_param* params = (newton_LSA_param*)options->solverParameters;
  NumericsMatrix* H = ((newton_LSA_data*)options->solverData)->H;
  assert(params);
  assert(H);

  search_data ls_data;
  ls_data.compute_F = functions->compute_F;
  ls_data.compute_F_merit = functions->compute_F_merit;
  ls_data.z = z;
  ls_data.zc = workV2;
  ls_data.F = F;
  ls_data.F_merit = F_merit;
  ls_data.desc_dir = workV1;
  /** \todo this value should be settable by the user with a default value*/
  ls_data.alpha_min = 1e-12;
  ls_data.alpha0 = 2.0;
  ls_data.data = data;
  ls_data.set = NULL;
  ls_data.sigma = params->sigma;
  ls_data.searchtype = LINESEARCH;

  if (options->iparam[SICONOS_IPARAM_LSA_FORCE_ARCSEARCH])
  {
    assert(functions->get_set_from_problem_data && \
        "newton_LSA :: arc search selected but no et_set_from_problem_data provided!");
    ls_data.set = functions->get_set_from_problem_data(data);
    ls_data.searchtype = ARCSEARCH;
  }

  nm_ref_struct nm_ref_data;
  if (options->iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS] > 0)
  {
    fill_nm_data(&nm_ref_data, options->iparam);
    ls_data.nm_ref_data = &nm_ref_data;
  }
  else
  {
    ls_data.nm_ref_data = NULL;
  }

  // if error based on the norm of JacThetaF_merit, do something not too stupid
  // here
  JacThetaF_merit[0] = DBL_MAX;

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
    int info_dir_search = 0;

    functions->compute_F(data, z, F);

    DEBUG_PRINT("z ");
    DEBUG_EXPR_WE(for (unsigned int i = 0; i < n; ++i)
        { DEBUG_PRINTF("% 2.2e ", z[i]) }
        DEBUG_PRINT("\n"));

    DEBUG_PRINT("F ");
    DEBUG_EXPR_WE(for (unsigned int i = 0; i < n; ++i)
        { DEBUG_PRINTF("% 2.2e ", F[i]) }
        DEBUG_PRINT("\n"));

    /**************************************************************************
     * START COMPUTATION DESCENTE DIRECTION
     */
    if (functions->compute_descent_direction)
    {
      functions->compute_descent_direction(data, z, F, workV1, options);
    }
    else
    {
      // Construction of H and F_desc
      if (functions->compute_RHS_desc) // different merit function for the descent calc.(usually min)
      {
        functions->compute_H_desc(data, z, F, workV1, workV2, H);
        functions->compute_RHS_desc(data, z, F, F_merit);
      } /* No computation of JacThetaFF_merit, this will come later */
      else
      {
        functions->compute_H(data, z, F, workV1, workV2, H);
        functions->compute_F_merit(data, z, F, F_merit);
        NM_tgemv(1., H, F_merit, 0., JacThetaF_merit);
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
      info_dir_search = NM_gesv(H, workV1);
    }
    /**************************************************************************
     * END COMPUTATION DESCENTE DIRECTION
     */

    DEBUG_PRINT("d_k ");
    DEBUG_EXPR_WE(for (unsigned int i = 0; i < n; ++i)
        { DEBUG_PRINTF("% 2.10e ", workV1[i]) }
        DEBUG_PRINT("\n"));

    /**************************************************************************
     * START COMPUTATION JacTheta F_merit
     */
    // Computation of JacThetaF_merit
    // JacThetaF_merit = H^T * F_merit
    if (functions->compute_RHS_desc || functions->compute_descent_direction) // need to compute H and F_merit for the merit
    {
      // /!\ maide! workV1 cannot be used since it contains the descent
      // direction ...

      if (functions->compute_JacTheta_merit)
      {
        functions->compute_JacTheta_merit(data, z, F, F_merit, workV2, JacThetaF_merit, options);
      }
      else
      {
        functions->compute_H(data, z, F, F_merit, workV2, H);
        functions->compute_F_merit(data, z, F, F_merit);
        NM_tgemv(1., H, F_merit, 0., JacThetaF_merit);
      }
    }


    DEBUG_PRINT("JacThetaF_merit ");
    DEBUG_EXPR_WE(for (unsigned int i = 0; i < n; ++i)
        { DEBUG_PRINTF("% 2.2e ", JacThetaF_merit[i]) }
        DEBUG_PRINT("\n"));


    // xhub :: we should be able to continue even if DGESV fails!
    if (info_dir_search)
    {
      if (functions->compute_RHS_desc) // we are safe here
      {
        DEBUG_PRINT("functions->compute_RHS_desc : no min descent direction ! searching for FB descent direction\n");
        cblas_dcopy(n, F_merit, incx, workV1, incy);
        cblas_dscal(n, -1.0, workV1, incx);
        info_dir_search = NM_gesv(H, workV1);
      }
      else
      {
        if (verbose > 0)
        {
          printf("Problem in DGESV, info = %d\n", info_dir_search);
        }
        options->iparam[1] = iter;
        options->dparam[1] = theta;
        *info = 2;

        goto newton_LSA_free;
      }
    }

    if (info_dir_search == 0) /* direction search succeeded */
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
    else /* direction search failed, backup to gradient step*/
    {
      cblas_dcopy(n, JacThetaF_merit, incx, workV1, incy);
      cblas_dscal(n, -1.0, workV1, incx);
    }

    tau = 1.0;
    if ((theta_iter > params->sigma * theta) || (info_dir_search > 0 && functions->compute_RHS_desc)) // Newton direction not good enough or min part failed
    {
      if (verbose > 1)
        printf("newton_LSA :: pure Newton direction not acceptable theta_iter = %g > %g = theta\n", theta_iter, theta);

      // Computations for the line search
      // preRHS = <JacThetaF_merit, d>
      // TODO: find a better name for this variable
      preRHS = cblas_ddot(n, JacThetaF_merit, incx, workV1, incy);

      // TODO we should not compute this if min descent search has failed
      threshold = -params->rho*pow(cblas_dnrm2(n, workV1, incx), params->p);
      //threshold = -rho*cblas_dnrm2(n, workV1, incx)*cblas_dnrm2(n, JacThetaF_merit, incx);
      if (preRHS > threshold)
      {
        if (verbose > 1)
          printf("newton_LSA :: direction not acceptable %g > %g\n", preRHS, threshold);
        cblas_dcopy(n, JacThetaF_merit, incx, workV1, incy);
        cblas_dscal(n, -1.0, workV1, incx);
        preRHS = cblas_ddot(n, JacThetaF_merit, incx, workV1, incy);
        DEBUG_PRINT("steepest descent ! d_k ");
        DEBUG_EXPR_WE(for (unsigned int i = 0; i < n; ++i)
            { DEBUG_PRINTF("% 2.2e ", workV1[i]) }
            DEBUG_PRINT("\n"));

      }
      preRHS *= params->gamma;

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

newton_LSA_free:

  if (!options->dWork)
  {
    free(JacThetaF_merit);
    free(F_merit);
    free(workV1);
    free(workV2);
  }
  if (ls_data.nm_ref_data)
  {
    free_nm_data((nm_ref_struct*)ls_data.nm_ref_data);
  }
}

void newton_lsa_default_SolverOption(SolverOptions* options)
{
  options->iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS] = 0;
  options->iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS_M] = 0;
  options->dparam[SICONOS_DPARAM_LSA_ALPHA_MIN] = 1e-12;
}

void set_lsa_params_data(SolverOptions* options, NumericsMatrix* mat)
{
  assert(mat);
  if (!options->solverParameters)
  {
    options->solverParameters = malloc(sizeof(newton_LSA_param));
    newton_LSA_param* params = (newton_LSA_param*) options->solverParameters;

    // gamma in (0, 1) or (0, 0.5) ??
    // inconsistency between Facchinei--Pang and "A Theoretical and Numerical
    // Comparison of Some Semismooth Algorithm for Complementarity Problems"
    // The following values are taken from the latter.
    params->p = 2.1;
    params->sigma = .9;
    params->gamma = 1e-4;
    params->rho = 1e-8;
  }

  if (!options->solverData)
  {
    options->solverData = malloc(sizeof(newton_LSA_data));
    newton_LSA_data* sd = (newton_LSA_data*) options->solverData;
    sd->H = duplicateNumericsMatrix(mat);
  }
}

bool newton_LSA_check_solverId(int solverId)
{
  switch (solverId)
  {
    case SICONOS_NCP_NEWTON_FBLSA:
    case SICONOS_NCP_NEWTON_MINFBLSA:
    case SICONOS_MCP_NEWTON_FBLSA:
    case SICONOS_MCP_NEWTON_MINFBLSA:
    case SICONOS_LCP_NEWTON_FBLSA:
    case SICONOS_LCP_NEWTON_MINFBLSA:
    case SICONOS_VI_BOX_QI:
    case SICONOS_VI_BOX_AVI_LSA:
      return true;
    default:
      return false;
  }
}

void newton_LSA_free_solverOptions(SolverOptions* options)
{
  if(options->solverParameters)
  {
    free(options->solverParameters);
    options->solverParameters = NULL;
  }

  if (options->solverData)
  {
    newton_LSA_data* sd = (newton_LSA_data*) options->solverData;
    assert(sd->H);
    freeNumericsMatrix(sd->H);
    free(sd->H);
    free(sd);
    options->solverData = NULL;
  }

}
