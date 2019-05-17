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

#include "debug.h"
#include "numerics_verbose.h"
#include "op3x3.h"
#include "SparseBlockMatrix.h"
#include "fc3d_Solvers.h"
#include "FrictionContactProblem.h"
#include "fc3d_compute_error.h"
#include "AlartCurnierGenerated.h"
#include "fc3d_nonsmooth_Newton_solvers.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "Friction_cst.h"
#include "VI_cst.h"
#include "SiconosLapack.h"

#define DEBUG_MESSAGES
#include "debug.h"

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
  case 0:
  {
    acparams.computeACFun3x3 = &computeAlartCurnierSTD;
    break;
  }
  case 1:
  {
    acparams.computeACFun3x3 = &computeAlartCurnierJeanMoreau;
    break;
  };
  case 2:
  {
    acparams.computeACFun3x3 = &fc3d_AlartCurnierFunctionGenerated;
    break;
  }
  case 3:
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
    SolverOptions * options_vi_eg =(SolverOptions *)malloc(sizeof(SolverOptions));
    fc3d_VI_ExtraGradient_setDefaultSolverOptions(options_vi_eg);
    options_vi_eg->iparam[SICONOS_IPARAM_MAX_ITER] = 50;
    options_vi_eg->dparam[SICONOS_DPARAM_TOL] = sqrt(options->dparam[0]);
    options_vi_eg->iparam[SICONOS_VI_IPARAM_ERROR_EVALUATION] = SICONOS_VI_ERROR_EVALUATION_LIGHT;
    fc3d_VI_ExtraGradient(problem, reaction , velocity , info , options_vi_eg);

    fc3d_nonsmooth_Newton_solvers_solve(&equation, reaction, velocity, info, options);
    solver_options_delete(options_vi_eg);
    free(options_vi_eg);
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

int fc3d_nonsmooth_Newton_AlartCurnier_setDefaultSolverOptions(
  SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the default solver options for the NSN_AC Solver\n");
  }

  options->solverId = SICONOS_FRICTION_3D_NSN_AC;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 20;
  options->dSize = 20;
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  options->dWork = NULL;
  solver_options_nullify(options);
  options->iparam[SICONOS_IPARAM_MAX_ITER] = 200;
  options->iparam[3] = 100000; /* nzmax*/
  options->iparam[5] = 1;
  options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] = 1;      /* erritermax */

  options->dparam[SICONOS_DPARAM_TOL] = 1e-3;
  options->iparam[SICONOS_FRICTION_3D_NSN_RHO_STRATEGY] = SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPECTRAL_NORM;
  options->dparam[SICONOS_FRICTION_3D_NSN_RHO] = 1;      /* default rho */

  options->iparam[8] = -1;     /* mpi com fortran */
  options->iparam[SICONOS_FRICTION_3D_NSN_FORMULATION] = SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_GENERATED;
  /* 0 STD AlartCurnier, 1 JeanMoreau, 2 STD generated, 3 JeanMoreau generated */
  options->iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH] = SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE;     /* 0 GoldsteinPrice line search, 1 FBLSA */
  options->iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH_MAXITER] = 100;   /* max iter line search */

#ifdef WITH_MUMPS
  options->iparam[13] = 1;
#else
  options->iparam[13] = 0;     /* Linear solver used at each Newton iteration. 0: cs_lusol, 1 mumps */
#endif

  options->internalSolvers = NULL;

  return 0;
}
