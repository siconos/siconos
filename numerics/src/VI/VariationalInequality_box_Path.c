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

#include <stdio.h>                          // for printf
#include "NumericsFwd.h"                    // for SolverOptions, Variationa...
#include "VariationalInequality_Solvers.h"  // for vi_box_path
#include "SiconosConfig.h" // for HAVE_PATHFERRIS // IWYU pragma: keep



#ifdef HAVE_PATHFERRIS
#include "VariationalInequality.h"
#include "SolverOptions.h"

#include <limits.h>
#include <assert.h>

#if defined(__cplusplus)
extern "C"
{
#endif
#include "PATH_SDK/include/MCP_Interface.h"

#include "PATH_SDK/include/Path.h"
#include "PATH_SDK/include/PathOptions.h"

#include "PATH_SDK/include/Macros.h"
#include "PATH_SDK/include/Output_Interface.h"
#include "PATH_SDK/include/Options.h"
#if defined(__cplusplus)
}
#endif

#include "SiconosSets.h"
#include "PathAlgebra.h"

#include "Path_interface.h"
#include "NumericsMatrix.h"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"

#if defined(__cplusplus)
#undef restrict
#define restrict __restrict
#endif

static CB_FUNC(void) PATH_problem_size(void* restrict id, int* restrict n, int* restrict nnz)
{
  SN_generic_path_env* env = (SN_generic_path_env*) id;
  *n = env->n;
  *nnz = env->nnz;
  return;
}

static CB_FUNC(void) vi_box_PATH_bounds(void* restrict id, int n, double* restrict z, double* restrict lb, double* restrict ub)
{
  SN_generic_path_env* env = (SN_generic_path_env*) id;

  for(unsigned i = 0; i < (unsigned)n; ++i)
  {
    z[i] = env->z[i];
    lb[i] = env->lb[i];
    ub[i] = env->ub[i];
  }
  return;
}

static CB_FUNC(int) vi_box_PATH_function_eval(void* id, int n, double*z, double *f)
{
  VariationalInequality* vi_box = (VariationalInequality*)((SN_generic_path_env*) id)->problem;

  vi_box->F(vi_box, n, z, f);
  return 0;
}

static CB_FUNC(int) vi_box_PATH_jacobian_eval(void *id, int n, double *z, int wantf,
    double *f, int *nnz,
    int *col_start, int *col_len,
    int *row, double *data)
{
  int err = 0;
  SN_generic_path_env* env = (SN_generic_path_env*) id;
  VariationalInequality* vi_box = (VariationalInequality*)env->problem;

  if(wantf)
  {
    //err +=
    vi_box->F(vi_box, n, z, f);
  }

  // err += to be added
  vi_box->compute_nabla_F(vi_box, n, z, vi_box->nabla_F);

  /* Write the jacobian in a Path-Sparse format */
  assert(vi_box->nabla_F->matrix0);
  assert(vi_box->nabla_F->storageType == NM_DENSE);
  convertToPathSparse(n, n, vi_box->nabla_F->matrix0, col_start, col_len, row, data);
  (*nnz) = env->nnz;
  return err;
}



void vi_box_path(VariationalInequality* problem, double *z, double* F, int *info, SolverOptions* options)
{
  assert(problem);
  unsigned n = problem->size;
  assert(n > 0);
  int nnz = n*n;

  assert(((generic_set*)problem->set)->id == SICONOS_SET_BOX);

  box_constraints* box_set = (box_constraints*) problem->set;

  SN_generic_path_env PATH_env =
  {
    n,
    nnz,
    z,
    F,
    box_set->lb,
    box_set->ub,
    problem
  };


  MCP_Interface mcp_interface =
  {
    &PATH_env,
    &PATH_problem_size, &vi_box_PATH_bounds,
    &vi_box_PATH_function_eval, &vi_box_PATH_jacobian_eval,
    NULL, /* hessian evaluation */
    NULL, NULL,
    NULL, NULL,
    NULL
  };

  SN_path_interface(&mcp_interface, z, F, info, options->dparam, options->iparam);

  return;
}

#else

void vi_box_path(VariationalInequality* problem, double *z, double* F, int *info, SolverOptions* options)
{
  printf("vi_box_path :: Path was not configured at compile time!\n");
}

#endif
