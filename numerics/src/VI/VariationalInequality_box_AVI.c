/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

#include <assert.h>                         // for assert
#include <stdlib.h>                         // for malloc, free, NULL
#include "Newton_methods.h"                 // for functions_LSA, init_lsa_f...
#include "NumericsFwd.h"                    // for RelayProblem, SolverOptions
#include "NumericsMatrix.h"                 // for NM_assert, NM_create_from...
#include "RelayProblem.h"                   // for RelayProblem, freeRelay_p...
#include "Relay_Solvers.h"                  // for relay_avi_caoferris
#include "relay_cst.h"                      // for SICONOS_RELAY_AVI_CAOFERRIS
#include "SiconosBlas.h"                    // for cblas_daxpy
#include "SiconosSets.h"                    // for box_constraints
#include "SolverOptions.h"                  // for SolverOptions, solver_opt...
#include "VI_Newton.h"                      // for VI_compute_F, VI_compute_...
#include "VariationalInequality.h"          // for VariationalInequality
#include "VariationalInequality_Solvers.h"  // for variationalInequality_BOX...
#include "sanitizer.h"                      // for cblas_dcopy_msan

typedef struct {
  NumericsMatrix* mat;
  RelayProblem* relay_pb;
} vi_box_AVI_LSA_data;

static int vi_compute_decent_dir_by_avi(void* problem, double* z, double* F, double* descent_dir, SolverOptions* options)
{
  VariationalInequality* vi_pb = (VariationalInequality*) problem;
  int n = vi_pb->size;
  vi_pb->F(vi_pb, n, z, F);
  RelayProblem* relay_pb = ((vi_box_AVI_LSA_data*)options->solverData)->relay_pb;

  NM_assert(NM_DENSE, relay_pb->M);

  vi_pb->compute_nabla_F(vi_pb, n, z, relay_pb->M);

  cblas_dcopy_msan(n, F, 1, relay_pb->q, 1);
  NM_gemv(-1.0, relay_pb->M, z, 1.0, relay_pb->q);

  int local_info = 0;
  relay_avi_caoferris(relay_pb, descent_dir, F, &local_info, options->internalSolvers[0]);

  cblas_daxpy(n, -1.0, z, 1, descent_dir, 1);

  return local_info;

//  avi_caoferris_stage3(avi_options, x, s_vec, problem->size, A,
//  options->internalSolver);
//  for (unsigned int i = 0; i<n; ++i) x[i] = s_vec[A[i]-1] + problem->lb[i];
}

void * vi_get_set(void* problem)
{
  return ((VariationalInequality*) problem)->set;
}


/** Release memory used by SolverOptions member solverData.

    Internal stuff. Must be called at the end of ncp_pathsearch.
    This concerns only options parameters which management is
    specific to this solver. All others (iparam ...) are handled
    in solver_options_delete generic function.
 */ 
static void vi_box_AVI_free(SolverOptions* options)
{
  if(options->solverData)
    {
      vi_box_AVI_LSA_data* sData = (vi_box_AVI_LSA_data*)options->solverData;
      NM_clear(sData->mat);
      free(sData->mat);
      sData->mat = NULL;
      sData->relay_pb->lb = NULL;
      sData->relay_pb->ub = NULL;
      freeRelay_problem(sData->relay_pb);
      free(sData);
    }
  options->solverData = NULL;
}

void vi_box_AVI_LSA(VariationalInequality* problem, double* z, double* F, int* info, SolverOptions* options)
{

  int n = problem->size;

  if (!options->solverData)
  {
    RelayProblem* relay_pb = (RelayProblem*)malloc(sizeof(RelayProblem));
    relay_pb->size = n;
    relay_pb->M = NM_create_from_data(NM_DENSE, n, n, malloc(n * n * sizeof(double)));;
    relay_pb->q = (double*) malloc(n * sizeof(double));

    box_constraints* box = (box_constraints*) problem->set;
    relay_pb->lb = box->lb;
    relay_pb->ub = box->ub;
    vi_box_AVI_LSA_data* sData = (vi_box_AVI_LSA_data*)malloc(sizeof(vi_box_AVI_LSA_data));
    sData->mat = (NumericsMatrix*)NM_duplicate(problem->nabla_F);
    sData->relay_pb = relay_pb;
    options->solverData = sData;
  }

  functions_LSA functions_AVI_LSA;
  init_lsa_functions(&functions_AVI_LSA, &VI_compute_F, &VI_compute_F_box_Qi);
  functions_AVI_LSA.compute_H = &VI_compute_H_box_Qi;
  functions_AVI_LSA.compute_error = &VI_compute_error_box;
  functions_AVI_LSA.compute_descent_direction = &vi_compute_decent_dir_by_avi;
  functions_AVI_LSA.get_set_from_problem_data = &vi_get_set;  
  set_lsa_params_data(options, problem->nabla_F);
  newton_LSA(problem->size, z, F, info, (void *)problem, options, &functions_AVI_LSA);

  vi_box_AVI_free(options);

}

void variationalInequality_BOX_AVI_set_default(SolverOptions* options)
{
  options->iparam[SICONOS_IPARAM_LSA_FORCE_ARCSEARCH] = 1;

  assert(options->numberOfInternalSolvers == 1);
  options->internalSolvers[0] = solver_options_create(SICONOS_RELAY_AVI_CAOFERRIS);
}

