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
#include "fc2d_compute_error.h"
#include <float.h>                   // for DBL_EPSILON
#include <math.h>                    // for fabs, sqrt
#include <stdio.h>                   // for printf
#include "FrictionContactProblem.h"  // for FrictionContactProblem
#include "NumericsMatrix.h"          // for NM_gemv
#include "numerics_verbose.h"        // for numerics_error, verbose
#include "SiconosBlas.h"                   // for cblas_dcopy

#define SGN(x) ((x) < 0 ? -1 : (x) > 0 ? 1 : 0)

int fc2d_compute_error(
  FrictionContactProblem* problem,
  double *z,
  double *w,
  double tolerance,
  double norm,
  double * error)
{

  /* Checks inputs */
  if(! problem || ! z || ! w)
    numerics_error("fc2d_compute_error", "null input for problem and/or z and/or w");

  int nc = problem->numberOfContacts;

  int n = nc * 2;

  int ic, iN, iT;

  double *mu = problem->mu;

  double tmp[2];

  double normT;

  cblas_dcopy(n, problem->q, 1, w, 1); // w <-q
  NM_gemv(1.0, problem->M, z, 1.0, w);

  *error = 0.;

  /* Num. Methods For Nonsmooth Dynamics, A.13 P 528 */
  /* DesaxceFeng98 */
  /* K* -) x _|_ y  (- K  <=>  x = projK(x-rho.y) for all rho>0 */

  for(ic = 0, iN = 0, iT = 1 ; ic < nc ; ++ic, ++iN, ++iN, ++iT, ++iT)
  {
    /* Compute the modified local velocity */
    tmp[0] = z[iN] - (w[iN] + mu[ic] * fabs(w[iT]));  /* rho=1 */
    tmp[1] = z[iT] - w[iT];                     /* rho=1 */

    /* projection */
    normT = fabs(tmp[1]);
    if(mu[ic]*normT <= -tmp[0])
    {
      tmp[0] = 0.;
      tmp[1] = 0.;
    }
    else if(normT > mu[ic]*tmp[0])
    {
      /* solve([sqrt((r1-mu*ra)^2+(r0-ra)^2)=abs(mu*r0-r1)/sqrt(mu*mu+1)],[ra]) */
      tmp[0] = (mu[ic] * normT + tmp[0]) / (mu[ic] * mu[ic] + 1);
      tmp[1] = mu[ic] * tmp[0] * SGN(tmp[1]);
    }

    tmp[0] = z[iN] -  tmp[0];
    tmp[1] = z[iT] -  tmp[1];
    *error += tmp[0] * tmp[0] + tmp[1] * tmp[1];

  }

  *error = sqrt(*error);
  if(fabs(norm) > DBL_EPSILON)
    *error /= norm;

  if(*error > tolerance)
  {
    if(verbose > 1) printf(" Numerics - fc2d_compute_error failed: error = %g > tolerance = %g.\n", *error, tolerance);
    return 1;
  }
  else
    return 0;
}
