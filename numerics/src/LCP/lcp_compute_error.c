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
#include "SiconosBlas.h"
#include "LinearComplementarityProblem.h"
#include "LCP_Solvers.h"
#include "NumericsMatrix.h"
#include "numerics_verbose.h"
#include <math.h>
#include <assert.h>
#include <float.h>
/* void lcp_compute_error_only(unsigned int n, double* restrict z , double* restrict w, double* restrict error) */
/* { */
/*   /\* Checks complementarity *\/ */

/*   *error = 0.; */
/*   double zi, wi; */
/*   for (unsigned int i = 0 ; i < n ; i++) */
/*   { */
/*     zi = z[i]; */
/*     wi = w[i]; */
/*     if (zi < 0.0) */
/*     { */
/*       *error += -zi; */
/*       if (wi < 0.0) *error += zi * wi; */
/*     } */
/*     if (wi < 0.0) *error += -wi; */
/*     if ((zi > 0.0) && (wi > 0.0)) *error += zi * wi; */
/*   } */
/* } */

void lcp_compute_error_only(unsigned int n, double* restrict z , double* restrict w, double* restrict error)
{
  *error = 0.;
  double zi, wi;
  for (unsigned int i = 0 ; i < n ; i++)
  {
    *error += pow(z[i] - fmax(0,(z[i] - w[i])),2);
  }
  *error =sqrt(*error);
}

int lcp_compute_error(LinearComplementarityProblem* problem, double *z , double *w, double tolerance, double * error)
{
  /* Checks inputs */
  if (problem == NULL || z == NULL || w == NULL)
    numerics_error("lcp_compute_error", "null input for problem and/or z and/or w");

  /* Computes w = Mz + q */
  int incx = 1, incy = 1;
  unsigned int n = problem->size;
  cblas_dcopy(n , problem->q , incx , w , incy);  // w <-q
  NM_gemv(1.0, problem->M, z, 1.0, w);
  double norm_q = cblas_dnrm2(n , problem->q , incx);
  lcp_compute_error_only(n, z, w, error);
  if (fabs(norm_q) > DBL_EPSILON)
    *error /= norm_q;

  
  //*error = *error / (norm_q + 1.0); /* Need some comments on why this is needed */
  if (*error > tolerance)
  {
    if (verbose > 0) printf(" Numerics - lcp_compute_error : error = %g > tolerance = %g.\n", *error, tolerance);
    return 1;
  }
  else
    return 0;

}
