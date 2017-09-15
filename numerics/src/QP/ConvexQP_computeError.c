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


#include "ConvexQP.h"
#include "SolverOptions.h"
#include "ConvexQP_computeError.h"
#include "NumericsMatrix.h"
#include "SiconosLapack.h"
#include "SiconosSets.h"

#include <math.h>
#include <assert.h>
#include <float.h>

/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"
#include "numerics_verbose.h"

int convexQP_computeError(
  ConvexQP* problem,
  double *z , double *w, double tolerance,
  SolverOptions * options, double * error)
{

  assert(problem);
  assert(z);
  assert(w);
  assert(error);
  
  int incx = 1;
  int n = problem->size;

  *error = 0.;
  if (!options->dWork)
  {
    options->dWork = (double*)calloc(2*n,sizeof(double));
  }
  double *ztmp =  options->dWork;
  double *wtmp =  &(options->dWork[n]);

  
  if (!problem->istheNormConvexQPset)
  {
    problem->normConvexQP= cblas_dnrm2(n , problem->q , 1);
    DEBUG_PRINTF("problem->norm ConvexQP= %12.8e\n", problem->normConvexQP);
    problem->istheNormConvexQPset=1;
  }

  double norm_q =problem->normConvexQP;
  DEBUG_PRINTF("norm_q = %12.8e\n", norm_q);

  /* q --> w */
  cblas_dcopy(n , problem->q , 1 , w, 1);

  /* M z + q --> w */
  NM_gemv(1.0, problem->M, z, 1.0, w);
  
  cblas_dcopy(n , z , 1 , ztmp, 1);
  cblas_daxpy(n, -1.0, w , 1, ztmp , 1) ;

  problem->ProjectionOnC(problem,ztmp,wtmp);

  cblas_daxpy(n, -1.0, z , 1, wtmp , 1) ;
  *error = cblas_dnrm2(n , wtmp , incx);

  /* Computes error */
  if (fabs(norm_q) > DBL_EPSILON)
    *error /= norm_q;

  DEBUG_PRINTF("error = %e\n",*error);
  if (*error > tolerance)
  {
    if (verbose > 1)
      printf(" Numerics - convexQP_compute_error: error = %g > tolerance = %g.\n",
             *error, tolerance);
    return 1;
  }
  else
    return 0;
}


