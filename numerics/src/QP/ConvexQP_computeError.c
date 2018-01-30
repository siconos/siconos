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
#include "NumericsVector.h"
#include <math.h>
#include <assert.h>
#include <float.h>

/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "debug.h"
#include "numerics_verbose.h"

int convexQP_compute_error_reduced(
  ConvexQP* problem,
  double *z , double *w, double tolerance,
  SolverOptions * options,
  double norm, double * error)
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

  /* q --> w */
  cblas_dcopy(n , problem->q , 1 , w, 1);

  /* M z + q --> w */
  NM_gemv(1.0, problem->M, z, 1.0, w);
  DEBUG_EXPR(NV_display(w,n));
  DEBUG_EXPR(NM_display(problem->M));
  cblas_dcopy(n , z , 1 , ztmp, 1);
  cblas_daxpy(n, -1.0, w , 1, ztmp , 1) ;
  DEBUG_EXPR(NV_display(w,n));
  problem->ProjectionOnC(problem,ztmp,wtmp);

  cblas_daxpy(n, -1.0, z , 1, wtmp , 1) ;
  *error = cblas_dnrm2(n , wtmp , incx);

  /* Computes error */
  if (fabs(norm) > DBL_EPSILON)
    *error /= norm;

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



int convexQP_compute_error(
  ConvexQP* problem,
  double *z , double * xi,
  double *w , double * u,
  double tolerance,
  double scaling,
  SolverOptions * options,
  double norm,
  double * error)
{

  assert(problem);
  assert(z);
  assert(u);
  assert(xi);
  assert(error);

  int incx = 1;
  int n = problem->size;
  int m = problem->m;


  double tmp = 0.0;

  *error = 0.;

  double *utmp =  options->dWork;
  double *utmp1 = &(options->dWork[m]) ;
  double *wtmp =  &(options->dWork[m+m]);




  DEBUG_EXPR(NV_display(z,n));
  DEBUG_EXPR(NV_display(xi,m));
  DEBUG_EXPR(NV_display(u,m));

  /* q --> w */
  cblas_dcopy(n , problem->q , 1 , w, 1);

  /* M z + q --> w */
  NM_gemv(1.0, problem->M, z, 1.0, w);

  DEBUG_EXPR(NV_display(w,n));

  /* Check that w= A^T xi */
  DEBUG_EXPR(NV_display(xi,m));
  
  cblas_dcopy(n , w , 1 , wtmp, 1);
  NM_tgemv(-scaling, problem->A, xi, 1.0, wtmp);
  DEBUG_EXPR(NV_display(wtmp,n));
  tmp = cblas_dnrm2(n , wtmp , incx);
  *error = tmp * tmp;
  DEBUG_PRINTF("square norm of Mz + q - A^T xi  = %e\n", *error);
  DEBUG_PRINTF("error = %e\n",*error);

  if (!problem->A)
  {
    cblas_dcopy(n,z,1,utmp,1);
  }
  else
  {
    /* b --> u */
    cblas_dcopy(m , problem->b , 1 , utmp, 1);
    /* A z + b --> u */
    NM_gemv(1.0, problem->A, z, 1.0, utmp);
  }

  cblas_daxpy(m, -1.0, u , 1, utmp , 1) ;
  tmp = cblas_dnrm2(m , utmp , incx);
  *error += tmp * tmp;

  /* Check that - xi \in \partial \Psi_C(u) */
  cblas_dcopy(m , u , 1 , utmp, 1);
  cblas_daxpy(m, -1.0, xi , 1, utmp , 1) ;

  problem->ProjectionOnC(problem,utmp,utmp1);

  DEBUG_EXPR(NV_display(utmp,m));
  DEBUG_EXPR(NV_display(utmp1,m));

  cblas_daxpy(m, -1.0, u , 1, utmp1 , 1) ;
  DEBUG_EXPR(NV_display(utmp1,m));

  tmp= cblas_dnrm2(m , utmp1 , incx);
  *error += tmp * tmp;
  DEBUG_PRINTF("error = %e\n",*error);
  /* Done, taking the square root */
  *error = sqrt(*error);

  DEBUG_PRINTF("norm = %12.8e\n", norm);
  if (fabs(norm) > DBL_EPSILON)
    *error /= norm;

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
