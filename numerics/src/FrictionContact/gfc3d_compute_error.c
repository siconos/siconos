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

#include "gfc3d_compute_error.h"

#include "GlobalFrictionContactProblem.h"
#include "gfc3d_Solvers.h"
#include "projectionOnCone.h"
#include "SiconosLapack.h"
#include <math.h>
#include <assert.h>
#include <float.h>
#include "sanitizer.h"
#include "numerics_verbose.h"
#include "NumericsMatrix.h"
int gfc3d_compute_error(GlobalFrictionContactProblem* problem, double* restrict reaction , double* restrict velocity, double* restrict globalVelocity, double tolerance, double* restrict error)
{

  /* Checks inputs */
  if (problem == NULL || reaction == NULL || velocity == NULL || globalVelocity == NULL)
    numerics_error("gfc3d_compute_error", "null input");

  gfc3d_init_workspace(problem);
  NumericsMatrix* factorized_M = problem->workspace->factorized_M;
  double* globalVelocitytmp = problem->workspace->globalVelocity;

  /* Computes error = dnorm2( GlobalVelocity -M^-1( q + H reaction)*/
  int nc = problem->numberOfContacts;
  int m = nc * 3;
  int n = problem->M->size0;
  double *mu = problem->mu;
  double *q = problem->q;
  NumericsMatrix *H = problem->H;

  cblas_dcopy_msan(n, q, 1, globalVelocitytmp, 1);

  NM_gemv(1.0, H, reaction, 1.0, globalVelocitytmp);

  CHECK_RETURN(!NM_gesv_expert(factorized_M, globalVelocitytmp, NM_KEEP_FACTORS));

  cblas_daxpy(n , -1.0 , globalVelocity , 1 , globalVelocitytmp, 1);

  /* We first accumulate the square terms and at the end we take the square
   * root */
  *error = cblas_ddot(n, globalVelocitytmp, 1, globalVelocitytmp, 1);

  cblas_dcopy(m, problem->b, 1, velocity, 1);
  NM_tgemv(1, H, globalVelocity, 1, velocity);

  double worktmp[3];
  double normUT;
  double rho = 1.0;
  for (int ic = 0 ; ic < nc ; ic++)
  {
    /* Compute the modified local velocity */
    normUT = sqrt(velocity[ic * 3 + 1] * velocity[ic * 3 + 1] + velocity[ic * 3 + 2] * velocity[ic * 3 + 2]);
    worktmp[0] = reaction[ic * 3] - rho * (velocity[ic * 3] + mu[ic] * normUT);
    worktmp[1] = reaction[ic * 3 + 1] - rho * velocity[ic * 3 + 1] ;
    worktmp[2] = reaction[ic * 3 + 2] - rho * velocity[ic * 3 + 2] ;
    projectionOnCone(worktmp, mu[ic]);
    worktmp[0] = reaction[ic * 3] -  worktmp[0];
    worktmp[1] = reaction[ic * 3 + 1] -  worktmp[1];
    worktmp[2] = reaction[ic * 3 + 2] -  worktmp[2];
    *error +=  worktmp[0] * worktmp[0] + worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2];
  }

  /* Done, taking the square root */
  *error = sqrt(*error);

  /* Computes error */
  double norm_q = cblas_dnrm2(n , problem->q , 1);
  *error = *error / (norm_q + 1.0);

  if (*error > tolerance)
  {
    /*       if (verbose > 0) printf(" Numerics - gfc3d_compute_error failed: error = %g > tolerance = %g.\n",*error, tolerance); */
    return 1;
  }
  else
    return 0;
}
