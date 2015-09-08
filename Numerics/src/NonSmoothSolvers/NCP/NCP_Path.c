/* Siconos-Numerics, Copyright INRIA 2005-2015
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
#include "NonlinearComplementarityProblem.h"
#include "SolverOptions.h"

#include "NumericsConfig.h"

#include "NCP_Solvers.h"

#ifdef HAVE_PATHFERRIS

#include <limits.h>
#include <assert.h>

#include "MCP_Interface.h"

#include "Path.h"
#include "PathOptions.h"

#include "Macros.h"
#include "Output_Interface.h"
#include "Options.h"

#include "PathAlgebra.h"

#include "Path_interface.h"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
#include "debug.h"

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

  for (unsigned i = 0; i < n; ++i) {
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
  printf("ncp_path :: Path was not configured at compile time!\n");
}

#endif
