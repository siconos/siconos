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
#include "RelayProblem.h"

#include "Relay_Solvers.h"

#include "sanitizer.h"
#include "numerics_verbose.h"
#include "NumericsMatrix.h"

void project_on_box(int n, double* restrict z, double* restrict lb, double* restrict ub)
{

  for (int i = 0 ; i < n ; ++i)
  {
    if (z[i] < lb[i]) z[i] = lb[i];
    else if (z[i] > ub[i]) z[i] = ub[i];
  }
}

int relay_compute_error(RelayProblem* problem, double* restrict z , double* restrict w, double tolerance, double* restrict error)
{
  /* Checks inputs */
  assert(problem);
  assert(z);
  assert(w);

  /* Computes w = Mz + q */
  int n = problem->size;
  cblas_dcopy(n , problem->q , 1 , w , 1);  // w <-q
  prodNumericsMatrix(n, n, 1.0, problem->M, z, 1.0, w);
  double * ztmp = (double*)malloc(n * sizeof(double));
  cblas_dcopy_msan(n , z , 1 , ztmp, 1);  // ztmp <-z

  double rho = -1.0;
  cblas_daxpy(n, rho, w, 1, ztmp, 1);    //ztmp <- ztmp - rho w

  project_on_box(n, ztmp , problem->lb, problem->ub);
  cblas_daxpy(n, -1.0, z, 1, ztmp, 1);    //ztmp <- ztmp -z


  *error = cblas_dnrm2(n , ztmp , 1);


  /* Computes error */
  double norm_q = cblas_dnrm2(n , problem->q , 1);
  *error = *error / (norm_q + 1.0);
  free(ztmp);
  if (*error > tolerance)
  {
    if (verbose > 0) printf(" Numerics - relay_compute_error: error = %g > tolerance = %g.\n", *error, tolerance);
    return 1;
  }
  else
    return 0;
}
