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

#include "fc3d_nonsmooth_Newton_natural_map.h"
#include <assert.h>                         // for assert
#include <math.h>                           // for sqrt
#include <stdio.h>                          // for printf, NULL
#include <stdlib.h>                         // for calloc, free, malloc
#include "FrictionContactProblem.h"         // for FrictionContactProblem
#include "Friction_cst.h"                   // for SICONOS_FRICTION_3D_NSN_NM
#include "NaturalMapGenerated.h"            // for fc3d_NaturalMapFunctionGe...
#include "SolverOptions.h"                  // for SolverOptions, solver_opt...
#include "fc3d_nonsmooth_Newton_solvers.h"  // for fc3d_nonsmooth_Newton_sol...
#include "numerics_verbose.h"               // for verbose
#include "fc3d_Solvers.h"

void fc3d_NaturalMapFunction(
  unsigned int problemSize,
  NaturalMapFun3x3Ptr computeACFun3x3,
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

  assert(problemSize / 3 > 0);
  assert(problemSize % 3 == 0);

  unsigned int i;
  for(i = 0; i < problemSize; i += 3)
  {

    computeACFun3x3(reaction, velocity, *mu, rho, result, A, B);

    reaction += 3;
    velocity += 3;
    mu++;
    rho += 3;

    if(result)
      result += 3;

    if(A)
      A += 9;

    if(B)
      B += 9;

  }

}


int fc3d_nonsmooth_Newton_NaturalMap_compute_error(
  FrictionContactProblem* problem,
  double *z, double *w, double tolerance,
  SolverOptions * options, double * error)
{

  double *A = NULL;
  double *B = NULL;

  unsigned int problemSize = 3 * problem->numberOfContacts;

  double *rho = (double*) malloc(problemSize*sizeof(double));
  double *F = (double *) malloc(problemSize*sizeof(double));

  NaturalMapFun3x3Ptr computeACFun3x3;

  switch(options->iparam[SICONOS_FRICTION_3D_NSN_FORMULATION])
  {
  case SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD:
  {

    computeACFun3x3 = &fc3d_NaturalMapFunctionGenerated;
    break;
  }
  }

  fc3d_NaturalMapFunction(
    problemSize,
    computeACFun3x3,
    z, w,
    problem->mu, rho,
    F, A, B);

  *error=0.;
  for(unsigned int i=0; i<problemSize;
      i+=3)
  {
    *error += sqrt(F[i]*F[i] + F[i+1]*F[i+1] + F[i+2]*F[i+2]);
  }

  *error /= (problem->numberOfContacts + 1);

  free(F);
  free(rho);

  if(*error > tolerance)
  {
    if(verbose > 1)
      printf(" Numerics - fc3d_compute_error: error = %g > tolerance = %g.\n",
             *error, tolerance);
    return 1;
  }
  else
  {
    return 0;
  }
}

void fc3d_nsn_nm_set_default(SolverOptions* options)
{
  options->iparam[1] = 1;
  options->iparam[3] = 100000; /* nzmax*/
  options->iparam[5] = 1;
  options->iparam[SICONOS_IPARAM_MAX_ITER] = 1;      /* erritermax */
  options->iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] = 1;

  options->dparam[SICONOS_FRICTION_3D_NSN_RHO] = 1.;

  options->iparam[SICONOS_FRICTION_3D_NSN_MPI_COM] = -1;
  options->iparam[SICONOS_FRICTION_3D_NSN_FORMULATION] = SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD;
  options->iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH] = SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE;
  options->iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH_MAX_ITER] = 100;

#ifdef WITH_MUMPS
  options->iparam[SICONOS_FRICTION_3D_NSN_LINEAR_SOLVER ] = SICONOS_FRICTION_3D_NSN_USE_MUMPS;
#else
  options->iparam[SICONOS_FRICTION_3D_NSN_LINEAR_SOLVER ] = SICONOS_FRICTION_3D_NSN_USE_CSLUSOL;
#endif
}


typedef struct
{
  NaturalMapFun3x3Ptr computeACFun3x3;
} NaturalMapParams;

void nonsmoothEqnNaturalMapFun(void* arg,
                               unsigned int problemSize,
                               double* reaction,
                               double* velocity,
                               double* mu,
                               double* rho,
                               double* result,
                               double* A,
                               double* B);
void nonsmoothEqnNaturalMapFun(void* arg,
                               unsigned int problemSize,
                               double* reaction,
                               double* velocity,
                               double* mu,
                               double* rho,
                               double* result,
                               double* A,
                               double* B)
{
  NaturalMapParams* acparams_p = (NaturalMapParams *) arg;

  fc3d_NaturalMapFunction(problemSize,
                          acparams_p->computeACFun3x3,
                          reaction,
                          velocity,
                          mu,
                          rho,
                          result,
                          A,
                          B);
}




void fc3d_nonsmooth_Newton_NaturalMap(
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

  NaturalMapParams acparams;

  switch(options->iparam[SICONOS_FRICTION_3D_NSN_FORMULATION])
  {
  case SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD:
  {
    acparams.computeACFun3x3 = &fc3d_NaturalMapFunctionGenerated;
    break;
  }
  }

  fc3d_nonsmooth_Newton_solvers equation;

  equation.problem = problem;
  equation.data = (void *) &acparams;
  equation.function = &nonsmoothEqnNaturalMapFun;

  fc3d_nonsmooth_Newton_solvers_solve(&equation, reaction, velocity, info,
                                      options);

}
