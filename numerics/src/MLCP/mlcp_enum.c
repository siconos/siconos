/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

/*
|A C| |u| |a| |0|
|   |*| |+| |=| |
|D B| |v| |b| |w|
0<z*v>0
dim(u)=mm
dim(v)=nn

**************************************************************************/

#include <assert.h>  // for assert
#include <float.h>   // for DBL_EPSILON
#include <math.h>    // for isinf, isnan
#ifndef __cplusplus
#include <stdbool.h>  // for false
#endif
#include <stdio.h>   // for printf
#include <stdlib.h>  // for malloc, free
#include <string.h>  // for memcpy

#include "MLCP_Solvers.h"    // for mlcp_compute_error
#include "NumericsMatrix.h"  // for NM_dense_display, Num...
#include "SiconosBlas.h"     // for cblas_dnrm2, cblas_ddot
#include "SiconosLapack.h"   // for DGELS, DGESV, DGESVD, lapack_int
#include "SolverOptions.h"   // for SolverOptions, SICONO...
#include "enum_tool.h"
#include "mlcp_cst.h"  // for SICONOS_IPARAM_MLCP_E...
#include "mlcp_enum.h"
#include "mlcp_enum_tool.h"  //
#include "numerics_errors.h"
#include "numerics_verbose.h"
#include "solver_registry.h"

/* #define DEBUG_MESSAGES */
#include "siconos_debug.h"
#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"  // for NV_display
#endif

/** Local, static functions **/

static void print_current_system(MixedLinearComplementarityProblem* problem,
                                 double* M_linear_system, double* q_linear_system) {
  int npm = problem->n + problem->m;
  int n_row = problem->M->size0;
  numerics_printf_verbose(2, "print_current_systemM:");
  NM_dense_display(M_linear_system, n_row, npm, 0);
  numerics_printf_verbose(
      2, "print_current_systemQ (ie -Q from mlcp because of linear system MZ=Q):");
  NM_dense_display(q_linear_system, n_row, 1, 0);
}

/* Minimum-norm least-squares solution of the square system A x = b (A: dim x dim, column
   major), by SVD, the singular values below dim * eps * s_max being discarded (the default
   rank threshold of LAPACK xGELSD and numpy lstsq). A is overwritten, b is overwritten by x.
   Fallback when DGESV finds the system of a case singular: redundant equality constraints
   (e.g. the 4 coplanar corners of a face of a rigid body resting on a plane, in a projection
   on the constraints) give a singular but consistent system, on which LU succeeds or fails
   depending on round-off (exact zero pivot or not). Whether the solution is accepted is
   decided afterwards by the usual checks (cone, mlcp_compute_error). Returns 0 on success. */
static int mlcp_enum_solve_svd(int dim, double* A, double* b) {
  double* S = (double*)malloc(dim * sizeof(double));
  double* U = (double*)malloc(dim * dim * sizeof(double));
  double* VT = (double*)malloc(dim * dim * sizeof(double));
  double* superb = (double*)malloc(dim * sizeof(double));
  double* tmp = (double*)malloc(dim * sizeof(double));
  lapack_int LAinfo = 0;
  DGESVD('A', 'A', dim, dim, A, dim, S, U, dim, VT, dim, superb, &LAinfo);
  if (!LAinfo) {
    const double threshold = dim * DBL_EPSILON * S[0];
    int rank = 0;
    for (int i = 0; i < dim; i++) {
      tmp[i] = 0.0;
      if (S[i] > threshold) {
        rank++;
        tmp[i] = cblas_ddot(dim, U + i * dim, 1, b, 1) / S[i];  // (U^T b)_i / s_i
      }
    }
    for (int j = 0; j < dim; j++)  // x = V tmp
      b[j] = cblas_ddot(dim, VT + j * dim, 1, tmp, 1);
    numerics_printf_verbose(1, "mlcp_enum: singular system solved by SVD, rank %d / %d", rank,
                            dim);
  }
  free(S);
  free(U);
  free(VT);
  free(superb);
  free(tmp);
  return (int)LAinfo;
}

/* An adaptation of the enum algorithm, to manage the case of MLCP-block formalization
 */
static void mlcp_enum_block(MixedLinearComplementarityProblem* problem, double* z, double* w,
                            int* info, SolverOptions* options) {
  DEBUG_BEGIN(" mlcp_enum_block(...)\n");

  int npm = (problem->n) + (problem->m);

  int NRHS = 1;
  int check;
  lapack_int LAinfo = 0;

  *info = 0;

  assert(problem->M);
  assert(problem->M->dense);
  assert(problem->q);

  int n_row = problem->M->size0;

  int n = problem->n;
  int m = problem->m;

  assert(problem->M->size1 == n + m);

  double* w_e = w;
  double* u = z;

  double tol = options->dparam[SICONOS_DPARAM_TOL];
  int itermax = options->iparam[SICONOS_IPARAM_MAX_ITER];
  int useDGELS = options->iparam[SICONOS_IPARAM_MLCP_ENUM_USE_DGELS];

  /*  LWORK = 2*npm; LWORK >= max( 1, MN + max( MN, NRHS ) ) where MN = min(M,N)*/
  // verbose=1;
  numerics_printf_verbose(1, "mlcp_enum_block BEGIN, n %d m %d tol %lf", n, m, tol);

  double* M_linear_system = options->dWork;
  /*  q_linear_system = M_linear_system + npm*npm;*/
  double* q_linear_system = M_linear_system + (n + m) * n_row;
  /*  q_linear_system_ref = q_linear_system  + m + n;*/
  double* q_linear_system_ref = q_linear_system + n_row;

  // double * work_DGELS = q_linear_system_ref + n_row;

  for (int row = 0; row < n_row; row++) q_linear_system_ref[row] = -problem->q[row];

  int* zw_indices = options->iWork;
  lapack_int* ipiv = zw_indices + m;
  int* indexInBlock = ipiv + m + n;
  if (m == 0) indexInBlock = 0;
  *info = 0;
  mlcp_enum_build_indexInBlock(problem, indexInBlock);

  unsigned long long int nbCase = enum_compute_nb_cases(problem->m);

  if (itermax < (int)nbCase) {
    printf("mlcp_enum_block itermax %i nbCase %llu ", itermax, nbCase);
    numerics_warning("mlcp_enum_block",
                     "all the cases will not be enumerated since itermax < nbCase)");
  }

  /* Pass 0: plain enumeration. Pass 1, only when pass 0 found no solution: the cases on
     which DGESV finds a singular system are solved again by SVD (mlcp_enum_solve_svd).
     Singular cases are frequent in the enumeration and the SVD is costly: it is not tried
     before all the regular cases have failed. */
  for (int svd_pass = 0; svd_pass < (useDGELS ? 1 : 2); svd_pass++) {
    EnumerationStruct* enum_struct = enum_init(problem->m);
    int iter = itermax;
    while (enum_next(zw_indices, problem->m, enum_struct) && iter-- > 0) {
      mlcp_enum_build_M_Block(zw_indices, M_linear_system, problem->M->dense, n, m, n_row,
                              indexInBlock);

      /* copy q_ref in q */
      memcpy(q_linear_system, q_linear_system_ref, n_row * sizeof(double));

      if (verbose > 1) print_current_system(problem, M_linear_system, q_linear_system);
      if (useDGELS) {
        DGELS(LA_NOTRANS, n_row, npm, NRHS, M_linear_system, n_row, q_linear_system, n_row,
              &LAinfo);
        numerics_printf_verbose(1, "Solution of dgels");
        {
          NM_dense_display(q_linear_system, n_row, 1, 0);
        }
      } else {
        DGESV(npm, NRHS, M_linear_system, npm, ipiv, q_linear_system, npm, &LAinfo);
        numerics_printf_verbose(1, "Solution of dgesv");
        if (LAinfo != 0) {
          numerics_printf_verbose(1, "DGESV FAILED");
          if (LAinfo > 0 && svd_pass) {  // singular system (LU overwrote it): retry by SVD
            mlcp_enum_build_M_Block(zw_indices, M_linear_system, problem->M->dense, n, m,
                                    n_row, indexInBlock);
            memcpy(q_linear_system, q_linear_system_ref, n_row * sizeof(double));
            LAinfo = mlcp_enum_solve_svd(npm, M_linear_system, q_linear_system);
          }
        } else
          numerics_printf_verbose(1, "DGESV SUCCEED");

        if (verbose > 1) {
          NM_dense_display(q_linear_system, n_row, 1, 0);
        }
      }
      if (!LAinfo) {
        if (useDGELS) {
          int cc = 0;
          int ii;

          for (ii = 0; ii < npm; ii++) {
            if (isnan(q_linear_system[ii]) || isinf(q_linear_system[ii])) {
              numerics_printf_verbose(1, "DGELS FAILED");
              cc = 1;
              break;
            }
          }
          if (cc) continue;

          if (n_row > npm) {
            double residual = cblas_dnrm2(n_row - npm, q_linear_system + npm, 1);

            if (residual > tol || isnan(residual) || isinf(residual)) {
              numerics_printf_verbose(
                  1, "DGELS, optimal point doesn't satisfy AX=b, residual = %e", residual);
              continue;
            }
            numerics_printf_verbose(1, "DGELS, optimal point residual = %e", residual);
          }
        }

        numerics_printf_verbose(1, "Solving linear system success, solution in cone?");
        if (verbose > 1) {
          NM_dense_display(q_linear_system, n_row, 1, 0);
        }

        check = 1;
        for (int row = 0; row < m; row++) {
          if (q_linear_system[indexInBlock[row]] < -tol) {
            check = 0;
            break; /*out of the cone!*/
          }
        }
        if (!check)
          continue;
        else {
          double err;
          mlcp_enum_fill_solution_Block(u, w_e, n, m, n_row, zw_indices, q_linear_system,
                                        indexInBlock);
          mlcp_compute_error(problem, z, w, tol, &err);
          /*because it happens the LU leads to an wrong solution witout raise any error.*/
          if (err > 10 * tol) {
            numerics_printf_verbose(1, "LU no-error, but mlcp_compute_error out of tol: %e!",
                                    err);
            continue;
          }
          numerics_printf_verbose(1, "mlcp_enum_block find a solution err = %e!", err);
          if (verbose) {
            mlcp_enum_display_solution_Block(u, w_e, n, m, n_row, indexInBlock);
          }
          // options->iparam[1]=scurrent-1;
          numerics_printf_verbose(1, "mlcp_enum_block END");
          options->dparam[SICONOS_DPARAM_RESIDU] = err;
          free(enum_struct);
          return;
        }
      } else {
        numerics_printf_verbose(1, "LU factorization failed:\n");
      }
    }
    free(enum_struct);
  }
  *info = 1;
  double err;
  mlcp_compute_error(problem, z, w, tol, &err);
  numerics_printf_verbose(1, "mlcp_enum_block failed! with error = %e\n", err);
  DEBUG_END(" mlcp_enum_block(...)\n");
}

/** End of static functions **/

int mlcp_enum_getNbIWork(MixedLinearComplementarityProblem* problem, SolverOptions* options) {
  CHECK_NULL(problem);
  return 2 * (problem->n + problem->m) + problem->m;
}

int mlcp_enum_getNbDWork(MixedLinearComplementarityProblem* problem, SolverOptions* options) {
  CHECK_NULL(problem);
  assert(problem->M);
  int LWORK = 0;
  if (options->iparam[SICONOS_IPARAM_MLCP_ENUM_USE_DGELS]) {
    LWORK = -1;
    // int info = 0;
    double dgelsSize = 0;
    // DGELS(problem->M->size0, problem->n + problem->m, 1, 0, problem->M->size0, 0,
    // problem->M->size0, &dgelsSize, LWORK, &info);
    LWORK = (int)dgelsSize;
  }
  return LWORK + 3 * (problem->M->size0) + (problem->n + problem->m) * (problem->M->size0);
}

void mlcp_enum(MixedLinearComplementarityProblem* problem, double* z, double* w, int* info,
               SolverOptions* options) {
  /* verbose=1; */
  int nbSol = 0;
  if (problem->blocksRows) {
    mlcp_enum_block(problem, z, w, info, options);
    return;
  }

  int* workingInt = options->iWork;
  int npm = (problem->n) + (problem->m);
  int NRHS = 1;
  lapack_int* ipiv;
  int check;
  lapack_int LAinfo = 0;
  *info = 0;

  /* sizes of the problem */
  int n_row = problem->M->size0;
  int n = problem->n;
  int m = problem->m;

  /* user parameters */
  double tol = options->dparam[SICONOS_DPARAM_TOL];
  int itermax = options->iparam[SICONOS_IPARAM_MAX_ITER];
  int useDGELS = options->iparam[SICONOS_IPARAM_MLCP_ENUM_USE_DGELS];

  /*OUTPUT param*/
  double* w_e = w;
  double* w_i = w + (n_row - problem->m); /*sW2 size :m */
  double* u = z;
  double* v = z + problem->n;

  /*  LWORK = 2*npm; LWORK >= max( 1, MN + max( MN, NRHS ) ) where MN = min(M,N)*/
  //  verbose=1;
  numerics_printf_verbose(1, "mlcp_enum BEGIN, n %d m %d tol %lf\n", n, m, tol);

  double* M_linear_system = options->dWork;
  /*  q_linear_system = M_linear_system + npm*npm;*/
  double* q_linear_system = M_linear_system + (n + m) * n_row;
  /*  q_linear_system_ref = q_linear_system + m + n;*/
  double* q_linear_system_ref = q_linear_system + n_row;

  // double * work_DGELS = q_linear_system_ref + n_row;

  for (int row = 0; row < n_row; row++) q_linear_system_ref[row] = -problem->q[row];

  int* zw_indices = workingInt;
  ipiv = zw_indices + m;
  /* Pass 0: plain enumeration. Pass 1, only when pass 0 found no solution: the cases on
     which DGESV finds a singular system are solved again by SVD (mlcp_enum_solve_svd).
     Singular cases are frequent in the enumeration and the SVD is costly: it is not tried
     before all the regular cases have failed. */
  for (int svd_pass = 0; svd_pass < (useDGELS ? 1 : 2); svd_pass++) {
    EnumerationStruct* enum_struct = enum_init(problem->m);
    int iter = itermax;
    while (enum_next(zw_indices, problem->m, enum_struct) && iter-- > 0) {
      mlcp_enum_build_M(zw_indices, M_linear_system, problem->M->dense, n, m, n_row);

      /* copy q_ref in q */
      memcpy(q_linear_system, q_linear_system_ref, n_row * sizeof(double));

      if (verbose > 1) print_current_system(problem, M_linear_system, q_linear_system);

      if (useDGELS) {
        DGELS(LA_NOTRANS, n_row, npm, NRHS, M_linear_system, n_row, q_linear_system, n_row,
              &LAinfo);
        numerics_printf_verbose(1, "Solution of dgels\n");
        if (verbose > 1) {
          NM_dense_display(q_linear_system, n_row, 1, 0);
        }
      } else {
        DGESV(npm, NRHS, M_linear_system, npm, ipiv, q_linear_system, npm, &LAinfo);
        numerics_printf_verbose(1, "Solution of dgesv\n");
        if (LAinfo > 0 && svd_pass) {  // singular system (LU overwrote it): retry by SVD
          mlcp_enum_build_M(zw_indices, M_linear_system, problem->M->dense, n, m, n_row);
          memcpy(q_linear_system, q_linear_system_ref, n_row * sizeof(double));
          LAinfo = mlcp_enum_solve_svd(npm, M_linear_system, q_linear_system);
        }
        if (verbose > 1) {
          NM_dense_display(q_linear_system, n_row, 1, 0);
        }
      }
      if (!LAinfo) {
        if (useDGELS) {
          int cc = 0;
          int ii;
          for (ii = 0; ii < npm; ii++) {
            if (isnan(q_linear_system[ii]) || isinf(q_linear_system[ii])) {
              numerics_printf_verbose(1, "DGELS FAILED\n");
              cc = 1;
              break;
            }
          }
          if (cc) continue;

          if (n_row > npm) {
            double residual = cblas_dnrm2(n_row - npm, q_linear_system + npm, 1);

            if (residual > tol || isnan(residual) || isinf(residual)) {
              numerics_printf_verbose(
                  1, "DGELS, optimal point doesn't satisfy AX=b, residual = %e\n", residual);
              continue;
            }
            numerics_printf_verbose(1, "DGELS, optimal point residual = %e\n", residual);
          }
        }

        numerics_printf_verbose(1, "Solving linear system success, solution in cone?\n");
        if (verbose > 1) {
          NM_dense_display(q_linear_system, n_row, 1, 0);
        }

        check = 1;
        for (int row = 0; row < m; row++) {
          if (q_linear_system[n + row] < -tol) {
            check = 0;
            break; /*out of the cone!*/
          }
        }
        if (!check)
          continue;
        else {
          double err;
          mlcp_enum_fill_solution(u, v, w_e, w_i, n, m, n_row, zw_indices, q_linear_system);
          mlcp_compute_error(problem, z, w, tol, &err);
          /*because it happens the LU leads to an wrong solution witout raise any error.*/
          if (err > 10 * tol) {
            numerics_printf_verbose(1, "LU no-error, but mlcp_compute_error out of tol: %e!\n",
                                    err);
            continue;
          }
          nbSol++;
          numerics_printf_verbose(1, "mlcp_enum find a solution, err=%e !\n", err);
          if (verbose > 1) {
            mlcp_enum_display_solution(u, v, w_e, w_i, n, m, n_row);
          }
          numerics_printf_verbose(1, "mlcp_enum END");
          free(enum_struct);
          return;
        }
      } else {
        numerics_printf_verbose(1, "LU factorization failed:\n");
      }
    }
    free(enum_struct);
  }
  *info = 1;
  numerics_printf_verbose(1, "mlcp_enum failed nbSol=%i!\n", nbSol);
}

void mlcp_enum_set_default(SolverOptions* options) {
  options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000000;
  options->dparam[SICONOS_IPARAM_MLCP_ENUM_USE_DGELS] = 0;
  options->filterOn = false;
}

/* ===========================================================================
 * Solver Registration
 * ===========================================================================
 * This registers SICONOS_MLCP_ENUM in the global solver registry.
 */

static int mlcp_enum_init_wrap(void* problem, SolverOptions* options) {
  /* set_default already called by solver_options_create */
  (void)problem;
  (void)options;
  return NUMERICS_OK;
}

static int mlcp_enum_solve_wrap(void* problem, double* z, double* w, SolverOptions* options) {
  int info = NUMERICS_OK;
  mlcp_enum((MixedLinearComplementarityProblem*)problem, z, w, &info, options);
  return info;
}

REGISTER_SOLVER(SICONOS_MLCP_ENUM, "MLCP_ENUM",
                "Enumerative solver for Mixed Linear Complementarity Problems",
                mlcp_enum_init_wrap, mlcp_enum_solve_wrap, NULL, /* free function */
                NULL,                                            /* error function */
                mlcp_enum_set_default,                           /* set_default */
                100000000,                                       /* default_max_iter */
                1e-6,                                            /* default_tol */
                0 /* is_local_solver */);
