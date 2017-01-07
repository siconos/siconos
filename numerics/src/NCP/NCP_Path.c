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

#include <stdio.h>

#include "NonlinearComplementarityProblem.h"
#include "SolverOptions.h"
#include "NCP_Solvers.h"
#include "sn_error_handling.h"

#ifdef HAVE_PATHFERRIS

#include <limits.h>
#include <assert.h>

#include "NumericsMatrix.h"
#include "numerics_verbose.h"

#include "PATH_SDK/include/MCP_Interface.h"
#include "PATH_SDK/include/Path.h"
#include "PATH_SDK/include/PathOptions.h"
#include "PATH_SDK/include/Macros.h"
#include "PATH_SDK/include/Output_Interface.h"
#include "PATH_SDK/include/Options.h"

#include "PathAlgebra.h"

#include "Path_interface.h"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"

#ifdef __cplusplus
#undef restrict
#define restrict __restrict
#endif

static CB_FUNC(void) ncp_PATH_problem_size(void* restrict id, int* restrict n, int* restrict nnz)
{
  SN_generic_path_env* env = (SN_generic_path_env*) id;
  *n = env->n;
  *nnz = env->nnz;
  return;
}

static CB_FUNC(void) ncp_PATH_bounds(void* restrict id, int n, double* restrict z, double* restrict lb, double* restrict ub)
{
  SN_generic_path_env* env = (SN_generic_path_env*) id;

  for (unsigned i = 0; i < (unsigned )n; ++i) {
    z[i] = env->z[i];
    lb[i] = 0.;
    ub[i] = 1e20; /* Magic value, see PATH documentation */
  }
  return;
}

static CB_FUNC(int) ncp_PATH_function_eval(void* id, int n, double*z, double *f)
{
  NonlinearComplementarityProblem* ncp = (NonlinearComplementarityProblem*)((SN_generic_path_env*) id)->problem;

  ncp->compute_F(ncp->env, n, z, f);
  return 0;
}

static CB_FUNC(int) ncp_PATH_jacobian_eval(void *id, int n, double *z, int wantf, 
                                        double *f, int *nnz,
                                        int *col_start, int *col_len, 
                                        int *row, double *data)
{
  int err = 0;
  SN_generic_path_env* env = (SN_generic_path_env*) id;
  NonlinearComplementarityProblem* ncp = (NonlinearComplementarityProblem*)env->problem;

  if (wantf) {
    //err += 
    ncp->compute_F(ncp->env, n, z, f);
  }

  // err += to be added
  ncp->compute_nabla_F(ncp->env, n, z, ncp->nabla_F);

  /* Write jacobianFGlocker in a Path-Sparse format */
  convertToPathSparse(n, n, ncp->nabla_F->matrix0, col_start, col_len, row, data);
  (*nnz) = env->nnz;
  return err;
}



void ncp_path(NonlinearComplementarityProblem* problem, double *z, double* F, int *info , SolverOptions* options)
{
  assert(problem);
  unsigned n = problem->n;
  assert(n > 0);
  int nnz = n*n;


  SN_generic_path_env PATH_env = {
    n,
    nnz,
    z,
    F,
    NULL,
    NULL,
    problem
  };


  MCP_Interface mcp_interface = {
    &PATH_env,
    &ncp_PATH_problem_size, &ncp_PATH_bounds,
    &ncp_PATH_function_eval, &ncp_PATH_jacobian_eval,
    NULL, /* hessian evaluation */
    NULL, NULL,
    NULL, NULL,
    NULL
  };

  SN_path_interface(&mcp_interface, z, F, info, options->dparam, options->iparam);

  return;
}

#else

void ncp_path(NonlinearComplementarityProblem* problem, double *z, double* F, int *info , SolverOptions* options)
{
  sn_fatal_error(SN_NOT_COMPILED_ERROR, "ncp_path :: Path was not configured at compile time!\n");
}

#endif
