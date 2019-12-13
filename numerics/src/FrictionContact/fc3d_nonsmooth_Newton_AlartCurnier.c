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

#include "fc3d_nonsmooth_Newton_AlartCurnier.h"
#include <assert.h>                         // for assert
#include <math.h>                           // for sqrt
#include <stddef.h>                         // for NULL
#include "AlartCurnierGenerated.h"          // for fc3d_AlartCurnierFunction...
#include "FrictionContactProblem.h"         // for FrictionContactProblem
#include "Friction_cst.h"                   // for SICONOS_FRICTION_3D_NSN_H...
#include "NumericsFwd.h"                    // for SolverOptions, FrictionCo...
#include "SolverOptions.h"                  // for SolverOptions, solver_opt...
#include "VI_cst.h"                         // for SICONOS_VI_ERROR_EVALUATI...
#include "fc3d_AlartCurnier_functions.h"    // for computeAlartCurnierJeanMo...
#include "fc3d_Solvers.h"                   // for fc3d_VI_ExtraGradient
#include "fc3d_nonsmooth_Newton_solvers.h"  // for fc3d_nonsmooth_Newton_sol...
#include "numerics_verbose.h"               // for numerics_error

void fc3d_AlartCurnierFunction(
  unsigned int problemSize,
  AlartCurnierFun3x3Ptr computeACFun3x3,
  double *reaction,
  double *velocity,
  double *mu,
  double *rho,
  double *result,
  double *A,
  double *B)
{
  assert(reaction);
  assert(velocity);
  assert(rho);
  assert(mu);

  //assert(problemSize / 3 > 0);
  assert(problemSize % 3 == 0);

  unsigned int i;
  for (i = 0; i < problemSize; i += 3)
  {

    computeACFun3x3(reaction, velocity, *mu, rho, result, A, B);

    reaction += 3;
    velocity += 3;
    mu++;
    rho += 3;

    if (result)
      result += 3;

    if (A)
      A += 9;

    if (B)
      B += 9;

  }
}


void nonsmoothEqnAlartCurnierFun(void* arg,
                                 unsigned int problemSize,
                                 double* reaction,
                                 double* velocity,
                                 double* mu,
                                 double* rho,
                                 double* result,
                                 double* A,
                                 double* B)
{
  AlartCurnierParams* acparams_p = (AlartCurnierParams *) arg;

  fc3d_AlartCurnierFunction(problemSize,
                                         acparams_p->computeACFun3x3,
                                         reaction,
                                         velocity,
                                         mu,
                                         rho,
                                         result,
                                         A,
                                         B);
}




void fc3d_nonsmooth_Newton_AlartCurnier(
  FrictionContactProblem* problem,
  double *reaction,
  double *velocity,
  int *info,
  SolverOptions *options)
{
  /* verbose=1; */
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

  switch (options->iparam[SICONOS_FRICTION_3D_NSN_FORMULATION])
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

  if(options->iparam[SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY] ==  SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_VI_EG_NSN)
  {
    SolverOptions * options_vi_eg = solver_options_create(SICONOS_FRICTION_3D_VI_EG);
    options_vi_eg->iparam[SICONOS_IPARAM_MAX_ITER] = 50;
    options_vi_eg->dparam[SICONOS_DPARAM_TOL] = sqrt(options->dparam[SICONOS_DPARAM_TOL]);
    options_vi_eg->iparam[SICONOS_VI_IPARAM_ERROR_EVALUATION] = SICONOS_VI_ERROR_EVALUATION_LIGHT;
    fc3d_VI_ExtraGradient(problem, reaction , velocity , info , options_vi_eg);

    fc3d_nonsmooth_Newton_solvers_solve(&equation, reaction, velocity, info, options);
    solver_options_delete(options_vi_eg);
    options_vi_eg = NULL;

  }
  else if (options->iparam[SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY] ==  SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NO)
  {
    fc3d_nonsmooth_Newton_solvers_solve(&equation, reaction, velocity, info, options);
  }
  else
  {
    numerics_error("fc3d_nonsmooth_Newton_AlartCurnier","Unknown nsn hybrid solver");
  }
}

void fc3d_nsn_ac_set_default(SolverOptions* options)
{
  options->iparam[3] = 100000; /* nzmax*/
  options->iparam[5] = 1;
  options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] = 1;
  options->iparam[SICONOS_FRICTION_3D_NSN_RHO_STRATEGY] = SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPECTRAL_NORM;
  options->dparam[SICONOS_FRICTION_3D_NSN_RHO] = 1;      /* default rho */

  options->iparam[SICONOS_FRICTION_3D_NSN_MPI_COM] = -1;     /* mpi com fortran */
  options->iparam[SICONOS_FRICTION_3D_NSN_FORMULATION] = SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_GENERATED;
  /* 0 STD AlartCurnier, 1 JeanMoreau, 2 STD generated, 3 JeanMoreau generated */
  options->iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH] = SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE;     /* 0 GoldsteinPrice line search, 1 FBLSA */
  options->iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH_MAX_ITER] = 100;   /* max iter line search */
  options->iparam[SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY] = SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NO;
#ifdef WITH_MUMPS
  options->iparam[SICONOS_FRICTION_3D_NSN_LINEAR_SOLVER] = SICONOS_FRICTION_3D_NSN_USE_MUMPS;
#else
  options->iparam[SICONOS_FRICTION_3D_NSN_LINEAR_SOLVER] = SICONOS_FRICTION_3D_NSN_USE_CSLUSOL;
#endif

}
