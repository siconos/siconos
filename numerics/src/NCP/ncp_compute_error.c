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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>
#include "SiconosBlas.h"
#include "NCP_Solvers.h"

#include "LCP_Solvers.h"

/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"


int ncp_compute_error(int n, double* z, double * F, double tol, double* err)
{
  DEBUG_BEGIN("ncp_compute_error(int n, double* z, double * F, double tol, double* err)\n")
  lcp_compute_error_only(n, z, F, err);

  DEBUG_PRINTF("ncp_compute_error err= %g, tol =%g \n", *err, tol);

  if (*err >= tol)
  {
    DEBUG_END("ncp_compute_error(int n, double* z, double * F, double tol, double* err)\n");
    return 1;
  }
  else
  {
    DEBUG_END("ncp_compute_error(int n, double* z, double * F, double tol, double* err)\n");
    return 0;
  }
}
