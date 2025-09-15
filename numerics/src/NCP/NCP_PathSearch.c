/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
#include "NCP_PathSearch.h"

#include <assert.h>  // for assert
#include <float.h>   // for DBL_EPSILON, DECIMAL_DIG
#include <math.h>    // for fmax, sqrt, fmin
#include <stdio.h>   // for printf, NULL
#include <stdlib.h>  // for malloc, free, exit, EXI...
#include <string.h>  // for memset

#include "ArmijoSearch.h"                     // for armijo_extra_params
#include "LCP_Solvers.h"                      // for lcp_compute_error, lcp_...
#include "LinearComplementarityProblem.h"     // for LinearComplementarityPr...
#include "NCP_Solvers.h"                      // for ncp_compute_error, ncp_...
#include "NMS.h"                              // for NMS_data, NMS, NMS_best...
#include "NSSTools.h"                         // for pos_part
#include "Newton_methods.h"                   // for functions_LSA, init_lsa...
#include "NonlinearComplementarityProblem.h"  // for NonlinearComplementarit...
#include "NumericsFwd.h"                      // for SolverOptions, Nonlinea...
#include "NumericsMatrix.h"                   // for NM_display, NM_clear
#include "PathSearch.h"                       // for pathsearch_data
#include "SiconosBlas.h"                      // for cblas_daxpy, cblas_ddot
#include "SiconosSets.h"                      // for set_set_id, SICONOS_SET...
#include "SolverOptions.h"                    // for SolverOptions, SICONOS_...
#include "lcp_cst.h"                          // for SICONOS_LCP_PIVOT_PATHS...
#include "lcp_pivot.h"                        // for LCP_PATHSEARCH_LEAVING_T
#include "line_search.h"                      // for search_data
#include "ncp_newton_FBLSA.h"                 // for FB_compute_F_ncp, FB_co...
#include "pivot-utils.h"                      // for lcp_pivot_diagnose_info
#include "sanitizer.h"                        // for cblas_dcopy_msan
#include "siconos_debug.h"                    // for DEBUG_PRINT, DEBUG_PRINTF

/** update the lcp subproblem: M, q and r
 * \param problem the NCP problem to solve
 * \param lcp_subproblem the lcp problem to fill
 * \param n size of the NCP problem
 * \param x_plus positive part of x
 * \param x current newton iterate
 * \param r value of the normal map
 */
static void ncp_pathsearch_update_lcp_data(NonlinearComplementarityProblem* problem,
                                           LinearComplementarityProblem* lcp_subproblem,
                                           unsigned n, double* restrict x_plus,
                                           double* restrict x, double* restrict r) {
  /* compute M = nabla F(x_plus) */
  problem->compute_nabla_F(problem->env, n, x_plus, problem->nabla_F);

  /* r = F_+(x) = F(x_+) + x - x_+ */
  /* the real q = q - r = x_+ - x - M x_plus */
  /* q = -M x_plus */
  NM_gemv(-1.0, problem->nabla_F, x_plus, 0.0, lcp_subproblem->q);

  /* first compute r = x - x_+ */
  cblas_dcopy(n, x, 1, r, 1);            /* r = x */
  cblas_daxpy(n, -1.0, x_plus, 1, r, 1); /* r -= x_plus */
  /* we factorized computations */
  cblas_daxpy(n, -1.0, r, 1, lcp_subproblem->q, 1); /* q -= x - x_plus */
  /* Finish r */
  problem->compute_F(problem->env, n, x_plus, x); /* compute F(x_plus) */
  cblas_daxpy(n, 1.0, x, 1, r, 1);                /* r += F(x_plus) */
}

/** Release memory used by SolverOptions members solverData and dWork.

    Internal stuff. Must be called at the end of ncp_pathsearch.
    This concerns only options parameters which management is
    specific to this solver. All others (iparam ...) are handled
    in solver_options_delete generic function.
 */
static void ncp_pathsearch_free(SolverOptions* opt) {
  if (opt->dWork) free(opt->dWork);
  opt->dWork = NULL;

  if (opt->solverData) {
    pathsearch_data* solverData_PathSearch = (pathsearch_data*)opt->solverData;
    free_NMS_data(solverData_PathSearch->data_NMS);
    free(solverData_PathSearch->lsa_functions);
    free(opt->solverData);
  }
  opt->solverData = NULL;
}

void ncp_pathsearch(NonlinearComplementarityProblem* problem, double* z, double* F, int* info,
                    SolverOptions* options) {
  /* Main step of the algorithm:
   * - compute jacobians
   * - call modified lemke
   */

  unsigned int n = problem->n;
  unsigned int preAlloc = options->iparam[SICONOS_IPARAM_PREALLOC];
  int itermax = options->iparam[SICONOS_IPARAM_MAX_ITER];

  DEBUG_EXPR(double merit_norm = 1.0;);
  double nn_tol = options->dparam[SICONOS_DPARAM_TOL];
  int nbiter = 0;

  /* declare a LinearComplementarityProblem on the stack*/
  LinearComplementarityProblem lcp_subproblem;
  lcp_subproblem.size = n;

  /* do some allocation if required
   * - nabla_F (used also as M for the LCP subproblem)
   * - q for the LCP subproblem
   *
   * Then fill the LCP subproblem
   */

  assert(problem->nabla_F);
  lcp_subproblem.M = problem->nabla_F;

  if (!preAlloc || (preAlloc && !options->dWork)) {
    options->dWork = (double*)malloc(4 * n * sizeof(double));
  }
  lcp_subproblem.q = options->dWork;
  double* x = &options->dWork[n];
  double* x_plus = &options->dWork[2 * n];
  double* r = &options->dWork[3 * n];

  NMS_data* data_NMS;
  functions_LSA* functions;

  if (!preAlloc || (preAlloc && !options->solverData)) {
    options->solverData = malloc(sizeof(pathsearch_data));
    pathsearch_data* solverData = (pathsearch_data*)options->solverData;

    /* do all the allocation */
    solverData->data_NMS = create_NMS_data(n, NM_DENSE, options->iparam, options->dparam);
    solverData->lsa_functions = (functions_LSA*)malloc(sizeof(functions_LSA));
    solverData->data_NMS->set = malloc(sizeof(positive_orthant));

    data_NMS = solverData->data_NMS;
    functions = solverData->lsa_functions;
    /* for use in NMS;  only those 3 functions are called */
    init_lsa_functions(functions, &FB_compute_F_ncp, &ncp_FB);
    functions->compute_H = &FB_compute_H_ncp;

    set_set_id(data_NMS->set, SICONOS_SET_POSITIVE_ORTHANT);

    /* fill ls_data */
    data_NMS->ls_data->compute_F = functions->compute_F;
    data_NMS->ls_data->compute_F_merit = functions->compute_F_merit;
    data_NMS->ls_data->z = NULL; /* XXX to check -- xhub */
    data_NMS->ls_data->zc = NMS_get_generic_workV(data_NMS->workspace, n);
    data_NMS->ls_data->F = NMS_get_F(data_NMS->workspace, n);
    data_NMS->ls_data->F_merit = NMS_get_F_merit(data_NMS->workspace, n);
    data_NMS->ls_data->desc_dir = NMS_get_dir(data_NMS->workspace, n);
    /** \todo this value should be settable by the user with a default value*/
    data_NMS->ls_data->alpha_min =
        fmin(data_NMS->alpha_min_watchdog, data_NMS->alpha_min_pgrad);
    data_NMS->ls_data->data = (void*)problem;
    data_NMS->ls_data->set = data_NMS->set;
    data_NMS->ls_data->sigma = options->dparam[SICONOS_DPARAM_NMS_SIGMA];
    armijo_extra_params* pG = (armijo_extra_params*)malloc(sizeof(armijo_extra_params));
    data_NMS->ls_data->extra_params = (void*)pG;
    search_Armijo_params_init(pG);
    pG->gamma = 1.;
    /* data_NMS->ls_data->searchtype is set in the NMS code */
  } else {
    pathsearch_data* solverData = (pathsearch_data*)options->solverData;
    data_NMS = solverData->data_NMS;
    functions = solverData->lsa_functions;
  }

  /* initial value for ref_merit */
  problem->compute_F(problem->env, n, z, F);
  functions->compute_F_merit(problem, z, F, data_NMS->ls_data->F_merit);

  data_NMS->ref_merit =
      .5 * cblas_ddot(n, data_NMS->ls_data->F_merit, 1, data_NMS->ls_data->F_merit, 1);
  data_NMS->merit_bestpoint = data_NMS->ref_merit;
  cblas_dcopy_msan(n, z, 1, NMS_checkpoint_0(data_NMS, n), 1);
  cblas_dcopy_msan(n, z, 1, NMS_checkpoint_T(data_NMS, n), 1);
  cblas_dcopy_msan(n, z, 1, NMS_bestpoint(data_NMS, n), 1);
  /* -------------------- end init ---------------------------*/

  int nms_failed = 0;
  double err = 10 * nn_tol;

  /* to check the solution */
  LinearComplementarityProblem lcp_subproblem_check;
  int check_lcp_solution = 1; /* XXX add config for that */

  double normal_norm2_newton_point;

  /* F is already computed here at z */

  while ((err > nn_tol) && (nbiter < itermax) && !nms_failed) {
    int force_watchdog_step = 0;
    int force_d_step_merit_check = 0;
    double check_ratio = 0.0;
    nbiter++;
    /* update M, q and r */

    /* First find x */
    ncp_pathsearch_compute_x_from_z(n, z, F, x);
    pos_part(n, x, x_plus); /* update x_plus */

    ncp_pathsearch_update_lcp_data(problem, &lcp_subproblem, n, x_plus, x, r);

    if (check_lcp_solution) {
      lcp_subproblem_check.size = n;
      lcp_subproblem_check.M = problem->nabla_F;
      lcp_subproblem_check.q = lcp_subproblem.q;
      // cblas_dcopy(n, x, 1, lcp_subproblem_check.q , 1);
      // NM_gemv(-1.0, problem->nabla_F, x_plus, 0.0, lcp_subproblem.q);
    }

    double norm_r2 = cblas_ddot(n, r, 1, r, 1);
    if (norm_r2 < DBL_EPSILON * DBL_EPSILON) /* ||r|| < 1e-15 */
    {
      DEBUG_PRINTF(
          "ncp_pathsearch :: ||r||  = %e < %e; path search procedure was successful!\n",
          norm_r2, DBL_EPSILON * DBL_EPSILON);
      (*info) = 0;
      ncp_compute_error(n, z, F, nn_tol,
                        &err); /* XXX F should be up-to-date, we should check only CC*/
      break;
    }

    /* end update M, q and r */

    lcp_pivot_covering_vector(&lcp_subproblem, x_plus, x, info, options->internalSolvers[0],
                              r);

    switch (*info) {
      case LCP_PIVOT_SUCCESS:
        DEBUG_PRINT("ncp_pathsearch :: path search procedure was successful!\n");
        if (check_lcp_solution) {
          double err_lcp = 0.0;
          cblas_daxpy(n, 1.0, r, 1, lcp_subproblem_check.q, 1);
          lcp_compute_error(&lcp_subproblem_check, x_plus, x, 1e-14, &err_lcp);
          double local_tol = fmax(1e-14, DBL_EPSILON * sqrt(norm_r2));
          printf("ncp_pathsearch :: lcp solved with error = %e; local_tol = %e\n", err_lcp,
                 local_tol);
          // assert(err_lcp < local_tol && "ncp_pathsearch :: lcp solved with very bad
          // precision");
          if (err_lcp > local_tol) {
            printf("ncp_pathsearch :: lcp solved with very bad precision\n");
            NM_display(lcp_subproblem.M);
            printf("z r q x_plus\n");
            for (unsigned i = 0; i < n; ++i)
              printf("%e %e %e %e\n", z[i], r[i], lcp_subproblem.q[i], x_plus[i]);
            options->internalSolvers[0]->iparam[SICONOS_LCP_IPARAM_PIVOTING_METHOD_TYPE] =
                0;  //
            // Note FP : no 0 in allowed values of the enum for PIVOTING METHOD
            lcp_pivot(&lcp_subproblem, x_plus, x, info, options->internalSolvers[0]);
            options->internalSolvers[0]->iparam[SICONOS_LCP_IPARAM_PIVOTING_METHOD_TYPE] =
                SICONOS_LCP_PIVOT_PATHSEARCH;
            lcp_compute_error(&lcp_subproblem_check, x_plus, x, 1e-14, &err_lcp);
            printf("ncp_pathsearch :: lcp resolved with error = %e; local_tol = %e\n", err_lcp,
                   local_tol);
          }

          /* XXX missing recompute x ?*/
          /* recompute the normal norm */
          problem->compute_F(problem->env, n, x_plus, r);
          cblas_daxpy(n, -1.0, x, 1, r, 1);
          normal_norm2_newton_point = cblas_ddot(n, r, 1, r, 1);
          if (normal_norm2_newton_point > norm_r2) {
            printf(
                "ncp_pathsearch :: lcp successfully solved, but the norm of the normal map "
                "increased! %e > %e\n",
                normal_norm2_newton_point, norm_r2);
            // assert(normal_norm2_newton_point <= norm_r2);
          } else {
            printf(
                "ncp_pathsearch :: lcp successfully solved, norm of the normal map decreased! "
                "%e < %e\n",
                normal_norm2_newton_point, norm_r2);
            // check_ratio = norm_r2/normal_norm2_newton_point;
          }
          if (50 * normal_norm2_newton_point < norm_r2) {
            force_d_step_merit_check = 1;
          } else if (10 * normal_norm2_newton_point < norm_r2) {
            //            check_ratio = sqrt(norm_r2/normal_norm2_newton_point);
          }
        }
        break;
      case LCP_PIVOT_RAY_TERMINATION:
        DEBUG_PRINT("ncp_pathsearch :: ray termination, let's fastened your seat belt!\n");
        break;
      case LCP_PATHSEARCH_LEAVING_T:
        DEBUG_PRINT("ncp_pathsearch :: leaving t, fastened your seat belt!\n");
        DEBUG_PRINTF("ncp_pathsearch :: max t value = %e\n",
                     options->internalSolvers[0]->dparam[2]); /* XXX fix 2 */
        /* try to retry solving the problem */
        /* XXX keep or not ? */
        /* recompute the normal norm */
        problem->compute_F(problem->env, n, x_plus, r);
        cblas_daxpy(n, -1.0, x, 1, r, 1);
        normal_norm2_newton_point = cblas_ddot(n, r, 1, r, 1);
        if (normal_norm2_newton_point > norm_r2) {
          printf(
              "ncp_pathsearch :: lcp successfully solved, but the norm of the normal map "
              "increased! %e > %e\n",
              normal_norm2_newton_point, norm_r2);
          // assert(normal_norm2_newton_point <= norm_r2);
        } else {
          printf(
              "ncp_pathsearch :: lcp successfully solved, norm of the normal map decreased! "
              "%e < %e\n",
              normal_norm2_newton_point, norm_r2);
          check_ratio = 5.0 * norm_r2 / normal_norm2_newton_point;
        }
        if (options->internalSolvers[0]->dparam[2] > 1e-5) break;
        memset(x_plus, 0, sizeof(double) * n);
        problem->compute_F(problem->env, n, x_plus, r);
        ncp_pathsearch_compute_x_from_z(n, x_plus, r, x);
        ncp_pathsearch_update_lcp_data(problem, &lcp_subproblem, n, x_plus, x, r);
        lcp_pivot_covering_vector(&lcp_subproblem, x_plus, x, info,
                                  options->internalSolvers[0], r);
        if (*info == LCP_PIVOT_SUCCESS) {
          DEBUG_PRINT("ncp_pathsearch :: Lemke start worked !\n");
          double err_lcp = 0.0;
          cblas_daxpy(n, 1.0, r, 1, lcp_subproblem_check.q, 1);
          lcp_compute_error(&lcp_subproblem_check, x_plus, x, 1e-14, &err_lcp);
          double local_tol = fmax(1e-14, DBL_EPSILON * sqrt(norm_r2));
          printf("ncp_pathsearch :: lcp solved with error = %e; local_tol = %e\n", err_lcp,
                 local_tol);
          assert(err_lcp < local_tol);
        } else {
          NM_display(lcp_subproblem.M);
          printf("z r q x_plus\n");
          for (unsigned i = 0; i < n; ++i)
            printf("%e %e %e %e\n", z[i], r[i], lcp_subproblem.q[i], x_plus[i]);
          DEBUG_PRINT("ncp_pathsearch :: Lemke start did not succeeded !\n");
          lcp_pivot_diagnose_info(*info);
          if (*info == LCP_PATHSEARCH_LEAVING_T) {
            DEBUG_PRINTF("ncp_pathsearch :: max t value after Lemke start = %e\n",
                         options->internalSolvers[0]->dparam[2]);
          }
          options->internalSolvers[0]->iparam[SICONOS_LCP_IPARAM_PIVOTING_METHOD_TYPE] = 0;
          lcp_pivot(&lcp_subproblem, x_plus, x, info, options->internalSolvers[0]);
          options->internalSolvers[0]->iparam[SICONOS_LCP_IPARAM_PIVOTING_METHOD_TYPE] =
              SICONOS_LCP_PIVOT_PATHSEARCH;
          double err_lcp = 0.0;
          lcp_compute_error(&lcp_subproblem, x_plus, x, 1e-14, &err_lcp);
          printf("ncp_pathsearch :: lemke start resolved with info = %d; error = %e\n", *info,
                 err_lcp);
          printf("x_plus x_minus\n");
          for (unsigned i = 0; i < n; ++i) printf("%e %e\n", x_plus[i], x[i]);
          /* recompute the normal norm */
          problem->compute_F(problem->env, n, x_plus, r);
          cblas_daxpy(n, -1.0, x, 1, r, 1);
          double normal_norm2_newton_point = cblas_ddot(n, r, 1, r, 1);
          if (normal_norm2_newton_point > norm_r2) {
            printf(
                "ncp_pathsearch :: lcp successfully solved, but the norm of the normal map "
                "increased! %e > %e\n",
                normal_norm2_newton_point, norm_r2);
            // assert(normal_norm2_newton_point <= norm_r2);
          } else {
            printf(
                "ncp_pathsearch :: lcp successfully solved, norm of the normal map decreased! "
                "%.*e < %.*e\n",
                DECIMAL_DIG, normal_norm2_newton_point, DECIMAL_DIG, norm_r2);
          }
          if (100 * normal_norm2_newton_point < norm_r2) {
            force_d_step_merit_check = 1;
          }
        }
        break;
      case LCP_PIVOT_NUL:
        printf("ncp_pathsearch :: kaboom, kaboom still more work needs to be done\n");
        lcp_pivot_diagnose_info(*info);
        //        exit(EXIT_FAILURE);
        force_watchdog_step = 1;
        break;
      case LCP_PATHSEARCH_NON_ENTERING_T:
        DEBUG_PRINT(
            "ncp_pathsearch :: non entering t, something is wrong here. Fix the f****** "
            "code!\n");
        assert(0 && "ncp_pathsearch :: non entering t, something is wrong here\n");
        force_watchdog_step = 1;
        break;
      default:
        printf("ncp_pathsearch :: unknown code returned by the path search\n");
        exit(EXIT_FAILURE);
    }

    nms_failed = NMS(data_NMS, problem, functions, z, x_plus, force_watchdog_step,
                     force_d_step_merit_check, check_ratio);
    /* at this point z has been updated */

    /* recompute the normal norm */
    problem->compute_F(problem->env, n, z, F);
    functions->compute_F_merit(problem, z, F, data_NMS->ls_data->F_merit);

    /* XXX is this correct ? */
    DEBUG_EXPR(merit_norm = .5 * cblas_ddot(n, data_NMS->ls_data->F_merit, 1,
                                            data_NMS->ls_data->F_merit, 1););
    ncp_compute_error(n, z, F, nn_tol,
                      &err); /* XXX F should be up-to-date, we should check only CC*/
    DEBUG_PRINTF("ncp_pathsearch :: iter = %d, ncp_error = %e; merit_norm^2 = %e\n", nbiter,
                 err, merit_norm);
  }

  options->iparam[SICONOS_IPARAM_ITER_DONE] = nbiter;
  options->dparam[SICONOS_DPARAM_RESIDU] = err;
  if (nbiter == itermax) {
    *info = 1;
  } else if (nms_failed) {
    *info = 2;
  } else {
    *info = 0;
  }

  DEBUG_PRINTF(
      "ncp_pathsearch procedure finished :: info = %d; iter = %d; ncp_error = %e; "
      "merit_norm^2 = %e\n",
      *info, nbiter, err, merit_norm);

  if (!preAlloc) {
    NM_clear(problem->nabla_F);
    free(problem->nabla_F);
    problem->nabla_F = NULL;

    ncp_pathsearch_free(options);
  }
}
