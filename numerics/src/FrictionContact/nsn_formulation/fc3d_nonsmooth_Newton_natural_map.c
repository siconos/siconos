// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "fc3d_nonsmooth_Newton_natural_map.h"

#include <assert.h>  // for assert
#include <math.h>    // for sqrt
#include <stdio.h>   // for printf, NULL
#include <stdlib.h>  // for calloc, free, malloc

#include "FrictionContactProblem.h"   // IWYU pragma: keep
#include "FrictionContact_options.h"  // for FC3D_NSN_NM
#include "NaturalMapGenerated.h"      // for fc3d_NaturalMapFunctionGe...
#include "fc3d_Solvers.h"
#include "fc3d_nonsmooth_Newton_solvers.h"  // for fc3d_nonsmooth_Newton_sol...
#include "naming_conventions.h"
#include "numerics_verbose.h"

void fc3d_NaturalMapFunction(unsigned int problemSize, NaturalMapFun3x3Ptr computeACFun3x3,
                             double* reaction, double* velocity, double* mu, double* rho,
                             double* result, double* A, double* B) {
  assert(reaction);
  assert(velocity);
  assert(rho);
  assert(mu);

  assert(problemSize / 3 > 0);
  assert(problemSize % 3 == 0);

  unsigned int i;
  for (i = 0; i < problemSize; i += 3) {
    computeACFun3x3(reaction, velocity, *mu, rho, result, A, B);

    reaction += 3;
    velocity += 3;
    mu++;
    rho += 3;

    if (result) result += 3;

    if (A) A += 9;

    if (B) B += 9;
  }
}

int fc3d_nonsmooth_Newton_NaturalMap_compute_error(FrictionContactProblem* problem,
                                                   double* reaction, double* velocity,
                                                   double tolerance, SolverOptions* options,
                                                   double* error) {
  double* A = NULL;
  double* B = NULL;

  unsigned int problemSize = 3 * problem->numberOfContacts;

  double* rho = (double*)malloc(problemSize * sizeof(double));
  double* F = (double*)malloc(problemSize * sizeof(double));

  NaturalMapFun3x3Ptr computeACFun3x3;

  switch (options->iparam[SICONOS_FRICTION_3D_NSN_FORMULATION]) {
    case SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD: {
      computeACFun3x3 = &fc3d_NaturalMapFunctionGenerated;
      break;
    }
  }

  fc3d_NaturalMapFunction(problemSize, computeACFun3x3, reaction, velocity, problem->mu, rho,
                          F, A, B);

  *error = 0.;
  for (unsigned int i = 0; i < problemSize; i += 3) {
    *error += sqrt(F[i] * F[i] + F[i + 1] * F[i + 1] + F[i + 2] * F[i + 2]);
  }

  *error /= (problem->numberOfContacts + 1);

  free(F);
  free(rho);

  if (*error > tolerance) {
    if (verbose > 1)
      printf(" Numerics - fc3d_compute_error: error = %g > tolerance = %g.\n", *error,
             tolerance);
    return 1;
  } else {
    return 0;
  }
}

void fc3d_nonsmooth_newton_jean_moreau_set_default(SolverOptions* options) {
  options->iparam[1] = 1;
  options->iparam[3] = 100000; /* nzmax*/
  options->iparam[5] = 1;
  SOLVER_MAX_ITER(options) = 1; /* erritermax */
  options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] = 1;

  options->dparam[SICONOS_FRICTION_3D_NSN_RHO] = 1.;

  options->iparam[SICONOS_FRICTION_3D_NSN_MPI_COM] = -1;
  options->iparam[SICONOS_FRICTION_3D_NSN_FORMULATION] =
      SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD;
  options->iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH] =
      SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE;
  /* Levenberg--Marquardt fallback when the Newton iteration stalls: off by default (see
   * nonsmooth_newton_levenberg_marquardt.h) */
  options->iparam[SICONOS_FRICTION_3D_NSN_LM_FALLBACK] = 0;
  options->dparam[SICONOS_FRICTION_3D_NSN_LM_COEF] = 1.0;
  options->iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH_MAX_ITER] = 100;

#ifdef WITH_MUMPS
  options->iparam[SICONOS_FRICTION_3D_NSN_LINEAR_SOLVER] = SICONOS_FRICTION_3D_NSN_USE_MUMPS;
#else
  options->iparam[SICONOS_FRICTION_3D_NSN_LINEAR_SOLVER] = SICONOS_FRICTION_3D_NSN_USE_CSLUSOL;
#endif
}

typedef struct {
  NaturalMapFun3x3Ptr computeACFun3x3;
} NaturalMapParams;

int nonsmoothEqnNaturalMapFun(void* arg, unsigned int problemSize, double* reaction,
                              double* velocity, double* mu, double* rho, double* result,
                              double* A, double* B);
int nonsmoothEqnNaturalMapFun(void* arg, unsigned int problemSize, double* reaction,
                              double* velocity, double* mu, double* rho, double* result,
                              double* A, double* B) {
  NaturalMapParams* acparams_p = (NaturalMapParams*)arg;

  fc3d_NaturalMapFunction(problemSize, acparams_p->computeACFun3x3, reaction, velocity, mu,
                          rho, result, A, B);
  return 0;
}

void fc3d_nonsmooth_Newton_NaturalMap(FrictionContactProblem* problem, double* reaction,
                                      double* velocity, int* info, SolverOptions* options) {
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

  assert(!options->iparam[4]);  // only host

  NaturalMapParams acparams;

  switch (options->iparam[SICONOS_FRICTION_3D_NSN_FORMULATION]) {
    case SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD: {
      acparams.computeACFun3x3 = &fc3d_NaturalMapFunctionGenerated;
      break;
    }
  }

  fc3d_nonsmooth_Newton_solvers equation;

  equation.problem = problem;
  equation.data = (void*)&acparams;
  equation.function = &nonsmoothEqnNaturalMapFun;

  fc3d_nonsmooth_Newton_solvers_solve(&equation, reaction, velocity, info, options);
}

/* Registration is in fc3d_nonsmooth_Newton_solvers.c to avoid duplicate registrations */
