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

#include "Newton_methods.h"                      // for functions_LSA, init_...
#include "NumericsFwd.h"                         // for VariationalInequality
#include "Qi_merit.h"                            // for Jac_F_Qi, phi_Qi
#include "SiconosSets.h"                         // for box_constraints
#include "SolverOptions.h"                       // for SolverOptions
#include "VI_Newton.h"                           // for VI_compute_F, VI_com...
#include "VariationalInequality.h"               // for VariationalInequality
#include "VariationalInequality_Solvers.h"       // for variationalInequalit...
#include "VariationalInequality_computeError.h"  // for variationalInequalit...

void VI_compute_F(void* data_opaque, double* x, double* F) {
  VariationalInequality* problem = (VariationalInequality*)data_opaque;
  problem->F(problem, problem->size, x, F);
}

void VI_compute_error_box(void* data_opaque, double* x, double* F, double* Jac_F_merit,
                          double tol, double* err) {
  VariationalInequality* problem = (VariationalInequality*)data_opaque;
  variationalInequality_compute_error_box(problem, x, F, tol, err);
}

void VI_compute_F_box_Qi(void* data_opaque, double* x, double* F, double* Fbox) {
  VariationalInequality* problem = (VariationalInequality*)data_opaque;
  phi_Qi(problem->size, x, F, Fbox, ((box_constraints*)problem->set)->lb,
         ((box_constraints*)problem->set)->ub);
}

void VI_compute_H_box_Qi(void* data_opaque, double* x, double* F, double* workV1,
                         double* workV2, NumericsMatrix* H) {
  VariationalInequality* problem = (VariationalInequality*)data_opaque;
  problem->compute_nabla_F(problem, problem->size, x, problem->nabla_F);

  Jac_F_Qi(problem->size, x, F, workV1, workV2, problem->nabla_F,
           ((box_constraints*)problem->set)->lb, ((box_constraints*)problem->set)->ub, H);
}

void* vi_get_set(void* problem); /*XXX */

void variationalInequality_box_newton_QiLSA(VariationalInequality* problem, double* x,
                                            double* F, int* info, SolverOptions* options) {
  functions_LSA functions_QiLSA;
  init_lsa_functions(&functions_QiLSA, &VI_compute_F, &VI_compute_F_box_Qi);
  functions_QiLSA.compute_H = &VI_compute_H_box_Qi;
  functions_QiLSA.compute_error = &VI_compute_error_box;
  functions_QiLSA.get_set_from_problem_data = &vi_get_set;

  set_lsa_params_data(options, problem->nabla_F);
  newton_LSA(problem->size, x, F, info, (void*)problem, options, &functions_QiLSA);
}

void variationalInequality_BOX_QI_set_default(SolverOptions* options) {
  options->iparam[SICONOS_IPARAM_STOPPING_CRITERION] = SICONOS_STOPPING_CRITERION_USER_ROUTINE;
  options->iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS] = 0;
  options->iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS_M] = 0;
  options->iparam[SICONOS_IPARAM_LSA_FORCE_ARCSEARCH] = 1;
  options->dparam[SICONOS_DPARAM_LSA_ALPHA_MIN] = 1e-16;
}
