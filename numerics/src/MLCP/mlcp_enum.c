/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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

/*
|A C| |u| |a| |0|
|   |*| |+| |=| |
|D B| |v| |b| |w|
0<z*v>0
dim(u)=mm
dim(v)=nn

**************************************************************************/


#include <assert.h>                             // for assert
#include <math.h>                               // for isinf, isnan
#ifndef __cplusplus
#include <stdbool.h>                       // for false
#endif
#include <stdio.h>                              // for printf
#include <string.h>                             // for memcpy

#include "MLCP_Solvers.h"                       // for mlcp_compute_error
#include "MixedLinearComplementarityProblem.h"  // for MixedLinearComplement...
#include "NumericsMatrix.h"                     // for NM_dense_display, Num...
#include "SiconosBlas.h"                        // for cblas_dnrm2
#include "SiconosLapack.h"                      // for DGELS, DGESV, lapack_int
#include "SolverOptions.h"                      // for SolverOptions, SICONO...
#include "mlcp_cst.h"                           // for SICONOS_IPARAM_MLCP_E...
#include "enum_tool.h"
#include "mlcp_enum_tool.h"                     //
#include "mlcp_enum.h"
#include "numerics_verbose.h"                   // for verbose
#include "SiconosConfig.h"                      // for MLCP_DEBUG // IWYU pragma: keep

/* #define DEBUG_MESSAGES */
#include "debug.h"
#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"                     // for NV_display
#endif


/** Local, static functions **/

static void print_current_system(MixedLinearComplementarityProblem* problem,
                               double * M_linear_system,
                               double * q_linear_system)
{
  int npm = problem->n + problem->m;
  int n_row = problem->M->size0;
  numerics_printf_verbose(2,"print_current_systemM:");
  NM_dense_display(M_linear_system, n_row, npm, 0);
  numerics_printf_verbose(2,"print_current_systemQ (ie -Q from mlcp because of linear system MZ=Q):");
  NM_dense_display(q_linear_system, n_row, 1, 0);
}

/* An adaptation of the enum algorithm, to manage the case of MLCP-block formalization
 */
static void mlcp_enum_block(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options)
{
  DEBUG_BEGIN(" mlcp_enum_block(...)\n");

  int npm = (problem->n) + (problem->m);


  int NRHS = 1;
  int check;
  lapack_int LAinfo = 0;

  *info = 0;

  assert(problem->M);
  assert(problem->M->matrix0);
  assert(problem->q);

  int n_row = problem->M->size0;

  int n = problem->n;
  int m = problem->m;

  assert(problem->M->size1 == n+m);

  double * w_e  = w;
  double * u = z;


  double tol = options->dparam[SICONOS_DPARAM_TOL];
  int itermax = options->iparam[SICONOS_IPARAM_MAX_ITER];
  int useDGELS = options->iparam[SICONOS_IPARAM_MLCP_ENUM_USE_DGELS];

  /*  LWORK = 2*npm; LWORK >= max( 1, MN + max( MN, NRHS ) ) where MN = min(M,N)*/
  //verbose=1;
  numerics_printf_verbose(1,"mlcp_enum_block BEGIN, n %d m %d tol %lf", n, m, tol);

  double * M_linear_system = options->dWork;
  /*  q_linear_system = M_linear_system + npm*npm;*/
  double * q_linear_system = M_linear_system + (n + m) * n_row;
  /*  q_linear_system_ref = q_linear_system  + m + n;*/
  double * q_linear_system_ref = q_linear_system  + n_row;

  // double * work_DGELS = q_linear_system_ref + n_row;

  for(int row = 0; row < n_row; row++)
    q_linear_system_ref[row] =  - problem->q[row];

  int * zw_indices = options->iWork;
  lapack_int * ipiv = zw_indices + m;
  int * indexInBlock = ipiv + m + n;
  if(m == 0)
    indexInBlock = 0;
  *info = 0;
  mlcp_enum_build_indexInBlock(problem, indexInBlock);

  EnumerationStruct * enum_struct = enum_init(problem->m);

  unsigned long long int nbCase =  enum_compute_nb_cases(problem->m);

  if(itermax < (int)nbCase)
  {
    numerics_warning("mlcp_enum_block", "all the cases will not be enumerated since itermax < nbCase)");
  }

  while(enum_next(zw_indices, problem->m, enum_struct) && itermax-- > 0)
  {
    mlcp_enum_build_M_Block(zw_indices, M_linear_system, problem->M->matrix0, n, m, n_row, indexInBlock);

    /* copy q_ref in q */
    memcpy(q_linear_system, q_linear_system_ref, n_row * sizeof(double));

    if(verbose >1)
      print_current_system(problem, M_linear_system, q_linear_system);
    if(useDGELS)
    {
      DGELS(LA_NOTRANS,n_row, npm, NRHS, M_linear_system, n_row, q_linear_system, n_row,&LAinfo);
      numerics_printf_verbose(1,"Solution of dgels");
      {
        NM_dense_display(q_linear_system, n_row, 1, 0);
      }
    }
    else
    {
      DGESV(npm, NRHS, M_linear_system, npm, ipiv, q_linear_system, npm, &LAinfo);
      numerics_printf_verbose(1,"Solution of dgesv");
      if(LAinfo != 0)
        numerics_printf_verbose(1,"DGESV FAILED");
      else
        numerics_printf_verbose(1,"DGESV SUCCEED");

      if(verbose > 1)
      {
        NM_dense_display(q_linear_system, n_row, 1, 0);
      }
    }
    if(!LAinfo)
    {
      if(useDGELS)
      {
        int cc = 0;
        int ii;

        for(ii = 0; ii < npm; ii++)
        {
          if(isnan(q_linear_system[ii]) || isinf(q_linear_system[ii]))
          {
            numerics_printf_verbose(1,"DGELS FAILED");
            cc = 1;
            break;
          }
        }
        if(cc)
          continue;

        if(n_row > npm)
        {
          double residual = cblas_dnrm2(n_row - npm, q_linear_system + npm, 1);

          if(residual > tol || isnan(residual) || isinf(residual))
          {
            numerics_printf_verbose(1,"DGELS, optimal point doesn't satisfy AX=b, residual = %e", residual);
            continue;
          }
          numerics_printf_verbose(1,"DGELS, optimal point residual = %e", residual);
        }
      }


      numerics_printf_verbose(1,"Solving linear system success, solution in cone?");
      if(verbose > 1)
      {
        NM_dense_display(q_linear_system, n_row, 1, 0);
      }

      check = 1;
      for(int row = 0 ; row < m; row++)
      {
        if(q_linear_system[indexInBlock[row]] < - tol)
        {
          check = 0;
          break;/*out of the cone!*/
        }
      }
      if(!check)
        continue;
      else
      {
        double err;
        mlcp_enum_fill_solution_Block(u, w_e, n, m, n_row, zw_indices, q_linear_system, indexInBlock);
        mlcp_compute_error(problem, z, w, tol, &err);
        /*because it happens the LU leads to an wrong solution witout raise any error.*/
        if(err > 10 * tol)
        {
          numerics_printf_verbose(1,"LU no-error, but mlcp_compute_error out of tol: %e!", err);
          continue;
        }
        numerics_printf_verbose(1,"mlcp_enum_block find a solution err = %e!", err);
        if(verbose)
        {
          mlcp_enum_display_solution_Block(u, w_e, n, m, n_row, indexInBlock);
        }
        // options->iparam[1]=scurrent-1;
        numerics_printf_verbose(1,"mlcp_enum_block END");
        options->dparam[SICONOS_DPARAM_RESIDU] = err;
        return;
      }
    }
    else
    {

      numerics_printf_verbose(1,"LU factorization failed:\n");

    }
  }
  *info = 1;
  numerics_printf_verbose(1,"mlcp_enum_block failed!\n");
  DEBUG_END(" mlcp_enum_block(...)\n");
}


/** End of static functions **/

int mlcp_enum_getNbIWork(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{
  if(!problem)
    return 0;
  return 2 * (problem->n + problem->m) + problem->m;
}

int mlcp_enum_getNbDWork(MixedLinearComplementarityProblem* problem, SolverOptions* options)
{
  if(!problem)
    return 0;
  assert(problem->M);
  int LWORK = 0;
  if(options->iparam[SICONOS_IPARAM_MLCP_ENUM_USE_DGELS])
  {
    LWORK = -1;
    //int info = 0;
    double dgelsSize = 0;
    //DGELS(problem->M->size0, problem->n + problem->m, 1, 0, problem->M->size0, 0, problem->M->size0, &dgelsSize, LWORK, &info);
    LWORK = (int) dgelsSize;
  }
  return LWORK + 3 * (problem->M->size0) + (problem->n + problem->m) * (problem->M->size0);
}

void mlcp_enum(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options)
{
  /* verbose=1; */
  int nbSol = 0;
  if(problem->blocksRows)
  {
    mlcp_enum_block(problem, z, w, info, options);
    return;
  }

  int * workingInt = options->iWork;
  int npm = (problem->n) + (problem->m);
  int NRHS = 1;
  lapack_int * ipiv;
  int check;
  lapack_int LAinfo = 0;
  *info = 0;

  /* sizes of the problem */
  int n_row = problem->M->size0;
  int n  = problem->n;
  int m = problem->m;

  /* user parameters */
  double tol = options->dparam[SICONOS_DPARAM_TOL];
  int itermax = options->iparam[SICONOS_IPARAM_MAX_ITER];
  int useDGELS = options->iparam[SICONOS_IPARAM_MLCP_ENUM_USE_DGELS];


  /*OUTPUT param*/
  double * w_e = w;
  double * w_i = w + (n_row - problem->m); /*sW2 size :m */
  double * u = z;
  double * v = z + problem->n;

  /*  LWORK = 2*npm; LWORK >= max( 1, MN + max( MN, NRHS ) ) where MN = min(M,N)*/
  //  verbose=1;
  numerics_printf_verbose(1,"mlcp_enum BEGIN, n %d m %d tol %lf\n", n, m, tol);

  double * M_linear_system = options->dWork;
  /*  q_linear_system = M_linear_system + npm*npm;*/
  double * q_linear_system = M_linear_system + (n + m) * n_row;
  /*  q_linear_system_ref = q_linear_system + m + n;*/
  double * q_linear_system_ref = q_linear_system + n_row;

  // double * work_DGELS = q_linear_system_ref + n_row;

  for(int row = 0; row < n_row; row++)
    q_linear_system_ref[row] =  - problem->q[row];

  int * zw_indices = workingInt;
  ipiv = zw_indices + m;
  EnumerationStruct * enum_struct = enum_init(problem->m);

  while(enum_next(zw_indices, problem->m, enum_struct) && itermax-- > 0)
  {
    mlcp_enum_build_M(zw_indices, M_linear_system, problem->M->matrix0, n, m, n_row);

    /* copy q_ref in q */
    memcpy(q_linear_system, q_linear_system_ref, n_row * sizeof(double));

    if(verbose > 1)
      print_current_system(problem, M_linear_system, q_linear_system);

    if(useDGELS)
    {
      DGELS(LA_NOTRANS,n_row, npm, NRHS, M_linear_system, n_row, q_linear_system, n_row, &LAinfo);
      numerics_printf_verbose(1,"Solution of dgels\n");
      if(verbose > 1)
      {
        NM_dense_display(q_linear_system, n_row, 1, 0);
      }
    }
    else
    {
      DGESV(npm, NRHS, M_linear_system, npm, ipiv, q_linear_system, npm, &LAinfo);
      numerics_printf_verbose(1,"Solution of dgesv\n");
      if(verbose > 1)
      {
        NM_dense_display(q_linear_system, n_row, 1, 0);
      }
    }
    if(!LAinfo)
    {
      if(useDGELS)
      {
        int cc = 0;
        int ii;
        for(ii = 0; ii < npm; ii++)
        {
          if(isnan(q_linear_system[ii]) || isinf(q_linear_system[ii]))
          {
            numerics_printf_verbose(1,"DGELS FAILED\n");
            cc = 1;
            break;
          }
        }
        if(cc)
          continue;

        if(n_row > npm)
        {
          double residual = cblas_dnrm2(n_row - npm, q_linear_system + npm, 1);

          if(residual > tol || isnan(residual) || isinf(residual))
          {
            numerics_printf_verbose(1,"DGELS, optimal point doesn't satisfy AX=b, residual = %e\n", residual);
            continue;
          }
          numerics_printf_verbose(1,"DGELS, optimal point residual = %e\n", residual);
        }
      }

      numerics_printf_verbose(1,"Solving linear system success, solution in cone?\n");
      if(verbose > 1)
      {
        NM_dense_display(q_linear_system, n_row, 1, 0);
      }

      check = 1;
      for(int row = 0 ; row < m; row++)
      {
        if(q_linear_system[n + row] < - tol)
        {
          check = 0;
          break;/*out of the cone!*/
        }
      }
      if(!check)
        continue;
      else
      {
        double err;
        mlcp_enum_fill_solution(u, v, w_e, w_i, n, m, n_row, zw_indices, q_linear_system);
        mlcp_compute_error(problem, z, w, tol, &err);
        /*because it happens the LU leads to an wrong solution witout raise any error.*/
        if(err > 10 * tol)
        {
          numerics_printf_verbose(1,"LU no-error, but mlcp_compute_error out of tol: %e!\n", err);
          continue;
        }
        nbSol++;
        numerics_printf_verbose(1,"mlcp_enum find a solution, err=%e !\n", err);
        if(verbose >1)
        {
          mlcp_enum_display_solution(u, v, w_e, w_i, n, m, n_row);
        }
        numerics_printf_verbose(1,"mlcp_enum END");
        return;
      }
    }
    else
    {
      numerics_printf_verbose(1,"LU factorization failed:\n");
    }
  }
  *info = 1;
  numerics_printf_verbose(1,"mlcp_enum failed nbSol=%i!\n", nbSol);
}

void mlcp_enum_set_default(SolverOptions* options)
{
  options->iparam[SICONOS_IPARAM_MAX_ITER] = 10000000;
  options->dparam[SICONOS_IPARAM_MLCP_ENUM_USE_DGELS] = 0;
  options->filterOn = false;

}
