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

#include "SiconosConfig.h"
#include "SiconosCompat.h"

#include "Newton_methods.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "numerics_verbose.h"
#include "SiconosLapack.h"
#include "ArmijoSearch.h"
#include "GoldsteinSearch.h"
#include "SolverOptions.h"
#include "lcp_cst.h"
#include "NCP_cst.h"
#include "MCP_cst.h"
#include "VI_cst.h"
#include "Friction_cst.h"

#include "hdf5_logger.h"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"

typedef double (*linesearch_fptr)(int n, double theta, double preRHS, search_data*);

#ifdef __cplusplus
using namespace std;
#endif

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

  char* solver_opt = getenv("SICONOS_SPARSE_SOLVER");
  if (solver_opt)
  {
    NM_setSparseSolver(H, atoi(solver_opt));
  }

  search_data ls_data;
  linesearch_fptr linesearch_algo;

  ls_data.compute_F = functions->compute_F;
  ls_data.compute_F_merit = functions->compute_F_merit;
  ls_data.z = z;
  ls_data.zc = workV2;
  ls_data.F = F;
  ls_data.F_merit = F_merit;
  ls_data.desc_dir = workV1;
  /** \todo this value should be settable by the user with a default value*/
  ls_data.alpha_min = options->dparam[SICONOS_DPARAM_LSA_ALPHA_MIN];
  ls_data.alpha0 = 2.0;
  ls_data.data = data;
  ls_data.set = NULL;
  ls_data.sigma = params->sigma;
  ls_data.searchtype = LINESEARCH;
  ls_data.extra_params = NULL;

  if (options->iparam[SICONOS_IPARAM_LSA_SEARCH_CRITERION] == SICONOS_LSA_GOLDSTEIN)
  {
    goldstein_extra_params* pG = (goldstein_extra_params*)malloc(sizeof(goldstein_extra_params));
    ls_data.extra_params = (void*) pG;
    search_Goldstein_params_init(pG);

/*    if (options->iparam[SICONOS_IPARAM_GOLDSTEIN_ITERMAX])
    {
      pG->iter_max = options->iparam[SICONOS_IPARAM_GOLDSTEIN_ITERMAX];
    }
    if (options->dparam[SICONOS_DPARAM_GOLDSTEIN_C])
    {
      pG->c = options->dparam[SICONOS_DPARAM_GOLDSTEIN_C];
    }
    if (options->dparam[SICONOS_DPARAM_GOLDSTEIN_ALPHAMAX])
    {
      pG->alpha_max = options->dparam[SICONOS_DPARAM_GOLDSTEIN_ALPHAMAX];
    }*/
    linesearch_algo = &linesearch_Goldstein2;
  }
  else if (options->iparam[SICONOS_IPARAM_LSA_SEARCH_CRITERION] == SICONOS_LSA_ARMIJO)
  {
    armijo_extra_params* pG = (armijo_extra_params*)malloc(sizeof(armijo_extra_params));
    ls_data.extra_params = (void*) pG;
    search_Armijo_params_init(pG);
    linesearch_algo = &linesearch_Armijo2;
  }
  else
  {
    fprintf(stderr, "Newton_LSA :: unknown linesearch specified");
    linesearch_algo = &linesearch_Armijo2;
  }

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

  unsigned log_hdf5 = SN_logh5_loglevel(SN_LOGLEVEL_NO);

  const char* hdf5_filename = getenv("SICONOS_HDF5_NAME");
  if (!hdf5_filename) hdf5_filename = "test.hdf5";
  SN_logh5* logger_s = NULL;
  if (log_hdf5)
  {
    logger_s = SN_logh5_init(hdf5_filename, itermax);
    SN_logh5_scalar_uinteger(0, "version", logger_s->file);
  }

  // Newton Iteration
  while ((iter < itermax) && (err > tol))
  {
    ++iter;
    int info_dir_search = 0;

    functions->compute_F(data, z, F);

    if (log_hdf5)
    {
      SN_logh5_new_iter(iter, logger_s);
      SN_LOG_LIGHT(log_hdf5,SN_logh5_vec_double(n, z, "z", logger_s->group));
      SN_LOG_VEC(log_hdf5,SN_logh5_vec_double(n, F, "F", logger_s->group));
    }

    /**************************************************************************
     * START COMPUTATION DESCENTE DIRECTION
     */
    if (functions->compute_descent_direction)
    {
      info_dir_search = functions->compute_descent_direction(data, z, F, workV1, options);
    }
    else
    {
      // Construction of H and F_desc
      if (functions->compute_RHS_desc) // different merit function for the descent calc.(usually min)
      {
        functions->compute_H_desc(data, z, F, workV1, workV2, H);
        functions->compute_RHS_desc(data, z, F, F_merit);
        if (log_hdf5)
        {
          SN_LOG_MAT(log_hdf5, SN_logh5_NM(H, "H_desc", logger_s));
          SN_LOG_VEC(log_hdf5, SN_logh5_vec_double(n, F_merit, "F_merit_desc", logger_s->group));
        }

      } /* No computation of JacThetaFF_merit, this will come later */
      else
      {
        functions->compute_H(data, z, F, workV1, workV2, H);
        functions->compute_F_merit(data, z, F, F_merit);
        NM_tgemv(1., H, F_merit, 0., JacThetaF_merit);
        if (log_hdf5)
        {
          SN_LOG_MAT(log_hdf5,SN_logh5_NM(H, "H", logger_s));
          SN_LOG_VEC(log_hdf5,SN_logh5_vec_double(n, F_merit, "F_merit", logger_s->group));
        }

      }

      // Find direction by solving H * d = -F_desc
      cblas_dcopy(n, F_merit, incx, workV1, incy);
      cblas_dscal(n, -1.0, workV1, incx);
      info_dir_search = NM_gesv(H, workV1, params->keep_H);
    }
    /**************************************************************************
     * END COMPUTATION DESCENTE DIRECTION
     */

    if (!info_dir_search && log_hdf5)  SN_LOG_VEC(log_hdf5,SN_logh5_vec_double(n, workV1, "desc_direction", logger_s->group));

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
        if (log_hdf5)
        {
          SN_LOG_MAT(log_hdf5,SN_logh5_NM(H, "H", logger_s));
          SN_LOG_VEC(log_hdf5,SN_logh5_vec_double(n, F_merit, "F_merit", logger_s->group));
        }
      }
    }

    if (log_hdf5)
    {
      SN_LOG_VEC(log_hdf5,SN_logh5_vec_double(n, JacThetaF_merit, "JacThetaF_merit", logger_s->group));
      SN_LOG_SCALAR(log_hdf5,SN_logh5_scalar_integer(info_dir_search, "info_dir_search_solve", logger_s->group));
    }

    // xhub :: we should be able to continue even if DGESV fails!
    if (info_dir_search)
    {
      if (functions->compute_RHS_desc) // we are safe here
      {
        DEBUG_PRINT("functions->compute_RHS_desc : no  descent direction found! searching for merit descent direction\n");
        cblas_dcopy(n, F_merit, incx, workV1, incy);
        cblas_dscal(n, -1.0, workV1, incx);
        info_dir_search = NM_gesv(H, workV1, params->keep_H);

        if (log_hdf5)
        {
          SN_LOG_SCALAR(log_hdf5,SN_logh5_scalar_integer(info_dir_search, "info_dir_search_solve_meritdesc", logger_s->group));
          if (!info_dir_search) SN_LOG_VEC(log_hdf5,SN_logh5_vec_double(n, workV1, "desc_merit_direction", logger_s->group));
        }
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
      theta_iter = INFINITY;
    }

    tau = 1.0;
    if ((theta_iter > params->sigma * theta) || (info_dir_search > 0 && functions->compute_RHS_desc)) // Newton direction not good enough or min part failed
    {
      if (log_hdf5)
      {
        SN_LOG_SCALAR(log_hdf5,SN_logh5_scalar_double(theta_iter, "theta_iter", logger_s->group));
        SN_LOG_SCALAR(log_hdf5,SN_logh5_scalar_double(params->sigma * theta, "theta_iter_threshold", logger_s->group));
      }

      if (verbose > 1)
        printf("newton_LSA :: pure Newton direction not acceptable theta_iter = %g > %g = theta\n", theta_iter, theta);

      // Computations for the line search
      // preRHS = <JacThetaF_merit, d>
      // TODO: find a better name for this variable
      preRHS = cblas_ddot(n, JacThetaF_merit, incx, workV1, incy);

      // TODO we should not compute this if min descent search has failed
      threshold = -params->rho*pow(cblas_dnrm2(n, workV1, incx), params->p);
      //threshold = -rho*cblas_dnrm2(n, workV1, incx)*cblas_dnrm2(n, JacThetaF_merit, incx);

      if (log_hdf5)
      {
        SN_LOG_SCALAR(log_hdf5,SN_logh5_scalar_double(preRHS, "preRHS_newton", logger_s->group));
        SN_LOG_SCALAR(log_hdf5,SN_logh5_scalar_double(threshold, "preRHS_threshold", logger_s->group));
      }

      if (params->check_dir_quality && preRHS > threshold)
      {
        if (verbose > 1)
          printf("newton_LSA :: direction not acceptable %g > %g\n", preRHS, threshold);

        cblas_dcopy(n, JacThetaF_merit, incx, workV1, incy);
        cblas_dscal(n, -1.0, workV1, incx);
        preRHS = cblas_ddot(n, JacThetaF_merit, incx, workV1, incy);
      }

      if (log_hdf5)
      {
        SN_LOG_SCALAR(log_hdf5,SN_logh5_scalar_double(preRHS, "preRHS", logger_s->group));
      }

      // Line search
      tau = (*linesearch_algo)(n, theta, preRHS, &ls_data);
    }

    if (isfinite(tau))
      cblas_daxpy(n , tau, workV1 , incx , z , incy);     //  z + tau*d --> z
    else
      cblas_daxpy(n, 1., workV1 , incx , z , incy);        // hack (restart)

    // Construction of the RHS for the next iterate
    functions->compute_F(data, z, F);
    functions->compute_F_merit(data, z, F, F_merit);

    // Merit Evaluation
    theta = cblas_dnrm2(n , F_merit , incx);
    theta = 0.5 * theta * theta;

    // Error Evaluation
    functions->compute_error(data, z, F, JacThetaF_merit, tol, &err);

    if (log_hdf5)
    {
      SN_logh5_scalar_double(err, "error", logger_s->group);
      SN_logh5_scalar_double(tau, "tau", logger_s->group);
      SN_logh5_scalar_double(theta, "theta", logger_s->group);
      SN_logh5_end_iter(logger_s);
    }


    if (options->callback)
    {
      stats_iteration.merit_value = theta;
      stats_iteration.alpha = tau;
      stats_iteration.status = 0;
      options->callback->collectStatsIteration(options->callback->env, n, z, F, err, &stats_iteration);
    }

  }

  options->iparam[1] = iter;
  options->dparam[1] = err;


  if (verbose > 0)
  {
    if (err > tol)
    {
      printf(" No convergence of the Newton algo after %d iterations\n" , iter);
      printf(" The residue is : %g \n", theta);
      *info = 1;
    }
    else
    {
      printf(" Convergence of the Newton algo after %d iterations\n" , iter);
      printf(" The residue is : %g \n", theta);
      *info = 0;
    }
  }
  else
  {
    if (err > tol) *info = 1;
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
  free_ls_data(&ls_data);

  if (log_hdf5)
  {
    SN_logh5_scalar_uinteger(iter, "nb_iter", logger_s->file);
    SN_logh5_scalar_double(err, "residual", logger_s->file);
    if (logger_s->group) SN_logh5_end_iter(logger_s);
    SN_logh5_end(logger_s);
  }
}

void newton_lsa_default_SolverOption(SolverOptions* options)
{
  options->iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS] = 0;
  options->iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS_M] = 0;
  options->dparam[SICONOS_DPARAM_LSA_ALPHA_MIN] = 0.;
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
    params->rho = 1e-8;

    params->keep_H = false;
    params->check_dir_quality = true;
  }

  if (!options->solverData)
  {
    options->solverData = malloc(sizeof(newton_LSA_data));
    newton_LSA_data* sd = (newton_LSA_data*) options->solverData;
    sd->H = NM_duplicate(mat);
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
    case SICONOS_FRICTION_3D_NSN_AC_TEST:
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
    NM_free(sd->H);
    free(sd->H);
    free(sd);
    options->solverData = NULL;
  }

}
