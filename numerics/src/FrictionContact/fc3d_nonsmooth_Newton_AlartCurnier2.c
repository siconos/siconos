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
#include <assert.h>                              // for assert
#include <math.h>                                // for sqrt
#ifndef __cplusplus
#include <stdbool.h>                             // for false, bool, true
#endif
#include <stdio.h>                               // for size_t, NULL
#include <stdlib.h>                              // for free, calloc, malloc
#include "SiconosBlas.h"                         // for cblas_dcopy, cblas_d...
#include "AlartCurnierGenerated.h"               // for fc3d_AlartCurnierFun...
#include "FrictionContactProblem.h"              // for FrictionContactProblem
#include "Friction_cst.h"                        // for SICONOS_FRICTION_3D_...
#include "Newton_methods.h"                      // for newton_LSA, function...
#include "NumericsFwd.h"                         // for FrictionContactProblem
#include "NumericsMatrix.h"                      // for NM_gemv, NumericsMatrix
#include "SolverOptions.h"                       // for SolverOptions, solve...
#include "VI_cst.h"                              // for SICONOS_VI_ERROR_EVA...
#include "fc3d_AlartCurnier_functions.h"         // for computeAlartCurnierJ...
#include "fc3d_Solvers.h"                        // for fc3d_VI_ExtraGradient
#include "fc3d_compute_error.h"                  // for fc3d_compute_error
#include "fc3d_nonsmooth_Newton_AlartCurnier.h"  // for AlartCurnierParams
#include "fc3d_nonsmooth_Newton_solvers.h"       // for fc3d_nonsmooth_Newto...
#include "line_search.h"                         // for SICONOS_LSA_GOLDSTEIN
#include "numerics_verbose.h"                    // for numerics_error

typedef struct
{
  FrictionContactProblem* problem;
  fc3d_nonsmooth_Newton_solvers* equation;
  double* rho;
  double* Ax;
  double* Bx;
  double normq;
  bool AwpB_data_computed;
} FC3D_Newton_data;

static void FC3D_compute_F(void* data_opaque, double* reaction, double* velocity)
{
  FrictionContactProblem* problem = ((FC3D_Newton_data*) data_opaque)->problem;
  // velocity <- M*reaction + qfree
  cblas_dcopy(problem->M->size1, problem->q, 1, velocity, 1);
  NM_gemv(1., problem->M, reaction, 1., velocity);
}

static void FC3D_compute_error(void* data_opaque, double* reaction, double* velocity, double* Jac_F_merit, double tol, double* err)
{
  FC3D_Newton_data* dat = (FC3D_Newton_data*) data_opaque;
  FrictionContactProblem* problem = dat->problem;

  fc3d_compute_error(problem, reaction, velocity, tol, NULL, dat->normq, err);
}

static void FC3D_compute_F_merit(void* data_opaque, double* reaction, double* velocity, double* F)
{
  FC3D_Newton_data* dat = (FC3D_Newton_data*) data_opaque;
  FrictionContactProblem* problem = dat->problem;
  fc3d_nonsmooth_Newton_solvers* equation = dat->equation;

  equation->function(equation->data, problem->M->size0, reaction, velocity, equation->problem->mu, dat->rho, F, dat->Ax, dat->Bx);
  dat->AwpB_data_computed = true;
}

static void FC3D_compute_AWpB(void* data_opaque, double* reaction, double* velocity, double* workV1, double* workV2, NumericsMatrix* AWpB)
{
  FC3D_Newton_data* dat = (FC3D_Newton_data*) data_opaque;
  assert(dat->AwpB_data_computed);
  FrictionContactProblem* problem = dat->problem;
  // AW + B
  computeAWpB(dat->Ax, problem->M, dat->Bx, AWpB);
  dat->AwpB_data_computed = false;

}





void fc3d_nonsmooth_Newton_AlartCurnier2(
  FrictionContactProblem* problem,
  double *reaction,
  double *velocity,
  int *info,
  SolverOptions *options)
{
  assert(problem);
  assert(reaction);
  assert(velocity);
  assert(info);
  assert(options);

  assert(problem->dimension == 3);

  assert(options->iparam);
  assert(options->dparam);

  assert(problem->q);
  assert(problem->mu);
  assert(problem->M);

  assert(!options->iparam[4]); // only host

  AlartCurnierParams acparams;

  switch(options->iparam[SICONOS_FRICTION_3D_NSN_FORMULATION])
  {
  case SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD:
  {
    acparams.computeACFun3x3 = &computeAlartCurnierSTD;
    break;
  }
  case SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_STD:
  {
    acparams.computeACFun3x3 = &computeAlartCurnierJeanMoreau;
    break;
  };
  case SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_GENERATED:
  {
    acparams.computeACFun3x3 = &fc3d_AlartCurnierFunctionGenerated;
    break;
  }
  case SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_GENERATED:
  {
    acparams.computeACFun3x3 = &fc3d_AlartCurnierJeanMoreauFunctionGenerated;
    break;
  }
  }

  fc3d_nonsmooth_Newton_solvers equation;

  equation.problem = problem;
  equation.data = (void *) &acparams;
  equation.function = &nonsmoothEqnAlartCurnierFun;

  /*************************************************************************
   * START NEW STUFF
   */
  size_t problemSize = problem->M->size0;
  size_t _3problemSize = problemSize + problemSize + problemSize;
  FC3D_Newton_data opaque_data;
  opaque_data.problem = problem;
  opaque_data.equation = &equation;
  opaque_data.rho = (double*)calloc(problemSize, sizeof(double));
  for(size_t i = 0; i < problemSize; ++i) opaque_data.rho[i] = 1.;
  opaque_data.Ax = (double*)calloc(_3problemSize, sizeof(double));
  opaque_data.Bx = (double*)calloc(_3problemSize, sizeof(double));
  opaque_data.normq = cblas_dnrm2(problemSize, problem->q, 1);
  opaque_data.AwpB_data_computed = false;

  functions_LSA functions_AC;
  init_lsa_functions(&functions_AC, &FC3D_compute_F, &FC3D_compute_F_merit);
  functions_AC.compute_H = &FC3D_compute_AWpB;
  functions_AC.compute_error = &FC3D_compute_error;
  functions_AC.get_set_from_problem_data = NULL;

  set_lsa_params_data(options, problem->M);
  newton_LSA_param* params = (newton_LSA_param*) options->solverParameters;
  params->check_dir_quality = false;

  options->iparam[SICONOS_IPARAM_LSA_SEARCH_CRITERION] = SICONOS_LSA_GOLDSTEIN;
//  options->iparam[SICONOS_IPARAM_LSA_SEARCH_CRITERION] = SICONOS_LSA_ARMIJO;
  /*************************************************************************
   * END NEW STUFF
   */

  if(options->iparam[SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY] ==  SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_VI_EG_NSN)
  {
    SolverOptions * options_vi_eg = solver_options_create(SICONOS_FRICTION_3D_VI_EG);
    options_vi_eg->iparam[SICONOS_IPARAM_MAX_ITER] = 50;
    options_vi_eg->dparam[SICONOS_DPARAM_TOL] = sqrt(options->dparam[SICONOS_DPARAM_TOL]);
    options_vi_eg->iparam[SICONOS_VI_IPARAM_ERROR_EVALUATION] = SICONOS_VI_ERROR_EVALUATION_LIGHT;
    fc3d_VI_ExtraGradient(problem, reaction, velocity, info, options_vi_eg);
    solver_options_delete(options_vi_eg);
    options_vi_eg = NULL;

    newton_LSA(problemSize, reaction, velocity, info, (void *)&opaque_data, options, &functions_AC);
  }
  else if(options->iparam[SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY] ==  SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NO)
  {
    newton_LSA(problemSize, reaction, velocity, info, (void *)&opaque_data, options, &functions_AC);
  }
  else
  {
    numerics_error("fc3d_nonsmooth_Newton_AlartCurnier","Unknown nsn hybrid solver");
  }

  free(opaque_data.rho);
  free(opaque_data.Ax);
  free(opaque_data.Bx);
}
