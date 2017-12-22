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
#include "fc3d_compute_error.h"

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

/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"
int gfc3d_compute_error(GlobalFrictionContactProblem* problem,
                        double*  reaction , double*  velocity,
                        double*  globalVelocity,
                        double tolerance,  double norm, double* restrict error)

{

  /* Checks inputs */
  if (problem == NULL || globalVelocity == NULL)
    numerics_error("gfc3d_compute_error", "null input");

  gfc3d_init_workspace(problem);

  double* tmp = problem->workspace->globalVelocity;

  
  /* Computes error = dnorm2( GlobalVelocity -M^-1( q + H reaction)*/
  int nc = problem->numberOfContacts;
  int m = nc * 3;
  int n = problem->M->size0;
  double *mu = problem->mu;
  double *q = problem->q;
  
  NumericsMatrix *H = problem->H;
  NumericsMatrix *M = problem->M;


  cblas_dcopy_msan(n, q, 1, tmp , 1);
  if (nc >0)
  {
    NM_gemv(1.0, H, reaction, 1.0, tmp);
  }

  NM_gemv(-1.0, M, globalVelocity, 1.0, tmp);
  *error = cblas_dnrm2(n,tmp,1);
  *error = *error* *error;
  DEBUG_PRINTF("square norm of -M v + H R + q = %e\n", *error);
  
  /* CHECK_RETURN(!NM_gesv_expert(problem->M, globalVelocity, NM_KEEP_FACTORS)); */
 

  if (nc >0)
  {
    /* Checks inputs */
    if (reaction == NULL || velocity == NULL)
      numerics_error("gfc3d_compute_error", "null input");

    cblas_dcopy(m, problem->b, 1, velocity, 1);
    NM_tgemv(1, H, globalVelocity, 1, velocity);

    double worktmp[3];
    for (int ic = 0 ; ic < nc ; ic++)
    {
      fc3d_unitary_compute_and_add_error(&reaction[ic * 3], &velocity[ic * 3], mu[ic],
                                         error,  worktmp);
    }
  }

  DEBUG_PRINTF("square of the error = %e\n", *error);
  
  /* Done, taking the square root */
  *error = sqrt(*error);

  DEBUG_PRINTF("error before normalization = %e\n", *error);
  
  DEBUG_PRINTF("norm = %12.8e\n", norm);
  if (fabs(norm) > DBL_EPSILON)
    *error /= norm;

  if (verbose)
  {
    if (tolerance * norm <  DBL_EPSILON)
      numerics_warning("gfc3d_compute_error", "The required relative precision (tolerance * norm = %e) is below  the machine accuracy", tolerance * norm);
  }

  
  if (*error > tolerance)
  {
    /*       if (verbose > 0) printf("Numerics - gfc3d_compute_error failed: error = %g > tolerance = %g.\n",*error, tolerance); */
    return 1;
  }
  else
    return 0;
}
