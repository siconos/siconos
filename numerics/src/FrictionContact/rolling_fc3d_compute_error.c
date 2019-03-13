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


#include "FrictionContactProblem.h"
#include "SolverOptions.h"
#include "rolling_fc3d_compute_error.h"
#include "rolling_fc3d_projection.h"
#include "projectionOnRollingCone.h"
#include "SiconosLapack.h"
#include "numerics_verbose.h"

#include <math.h>
#include <assert.h>
#include <float.h>

/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"
#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"
#endif


void rolling_fc3d_unitary_compute_and_add_error(
  double* restrict r,
  double* restrict u,
  double mu,
  double mur,
  double* restrict error,
  double * worktmp)
{
  DEBUG_BEGIN("rolling_fc3d_unitary_compute_and_add_error(...)\n");
  DEBUG_EXPR(NV_display(r,5););
  DEBUG_EXPR(NV_display(u,5););
  DEBUG_PRINTF(" tilde u[0] = %f\n", u[0] + mu *  hypot(u[1], u[2]) + mur * hypot(u[3], u[4] ));
  /* Compute the modified local velocity */
  //normUT = hypot(u[1], u[2]); // i.e sqrt(u[ic3p1]*u[ic3p1]+u[ic3p2]*u[ic3p2]);
  worktmp[0] = r[0] -  u[0] - mu *  hypot(u[1], u[2]) - mur * hypot(u[3], u[4]);
  worktmp[1] = r[1] -  u[1] ;
  worktmp[2] = r[2] -  u[2] ;
  worktmp[3] = r[3] -  u[3] ;
  worktmp[4] = r[4] -  u[4] ;
  DEBUG_PRINT("r-rho tilde v  before projection");DEBUG_EXPR(NV_display(worktmp,5););
  projectionOnRollingCone(worktmp, mu, mur);
  DEBUG_PRINT("after projection");DEBUG_EXPR(NV_display(worktmp,5););
  worktmp[0] = r[0] -  worktmp[0];
  worktmp[1] = r[1] -  worktmp[1];
  worktmp[2] = r[2] -  worktmp[2];
  worktmp[3] = r[3] -  worktmp[3];
  worktmp[4] = r[4] -  worktmp[4];
  DEBUG_EXPR(NV_display(worktmp,5););
  *error +=
    worktmp[0] * worktmp[0] +
    worktmp[1] * worktmp[1] +
    worktmp[2] * worktmp[2] +
    worktmp[3] * worktmp[3] +
    worktmp[4] * worktmp[4];
  DEBUG_END("rolling_fc3d_unitary_compute_and_add_error(...)\n");
}

int rolling_fc3d_compute_error(
  RollingFrictionContactProblem* problem,
  double *reaction , double *velocity, double tolerance,
  SolverOptions * options, double norm, double * error)
{
  DEBUG_BEGIN("rolling_fc3d_compute_error(...)\n");
  assert(problem);
  assert(reaction);
  assert(velocity);
  assert(error);

  /* Computes velocity = Mreaction + q */
  int incx = 1, incy = 1;
  int nc = problem->numberOfContacts;
  int n = nc * 5;
  double *mu = problem->mu;
  double *mur = problem->mu_r;

  cblas_dcopy(n , problem->q , incx , velocity , incy); // velocity <-q
  // Compute the current velocity
  //NM_prod_mv_3x3(n, n, problem->M, reaction, velocity);
  NM_gemv(1.0, problem->M, reaction,
               1.0,
               velocity);
  DEBUG_EXPR(NV_display(problem->q,n););
  DEBUG_EXPR(NV_display(velocity,n););
  DEBUG_EXPR(NV_display(reaction,n););
  
  *error = 0.;
  int ic, ic5;
  double worktmp[5];
  for (ic = 0, ic5 = 0 ; ic < nc ; ic++, ic5 += 5)
  {
    rolling_fc3d_unitary_compute_and_add_error(reaction + ic5, velocity + ic5, mu[ic], mur[ic], error, worktmp);
    DEBUG_PRINTF("squared absolute error = %12.8e contact =%i nc= %i\n", *error, ic, nc);
  }
  *error = sqrt(*error);
  DEBUG_PRINTF("absolute error = %12.8e\n", *error);
  /* Compute relative error with respect to norm */
  DEBUG_PRINTF("norm = %12.8e\n", norm);
  if (fabs(norm) > DBL_EPSILON)
    *error /= norm;
  /* *error = *error / (norm + 1.0); old version */
  DEBUG_PRINTF("relative error = %12.8e\n", *error);
  DEBUG_END("rolling_fc3d_compute_error(...)\n");
  if (*error > tolerance)
    return 1;

  return 0;
}


