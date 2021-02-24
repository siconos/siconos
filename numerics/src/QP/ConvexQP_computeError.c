/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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

#include "ConvexQP_computeError.h"
#include <assert.h>            // for assert
#include <float.h>             // for DBL_EPSILON
#include <math.h>              // for fabs, sqrt
#include <stdio.h>             // for printf
#include <stdlib.h>            // for calloc
#include "ConvexQP.h"          // for ConvexQP
#include "NumericsMatrix.h"    // for NM_gemv, NM_tgemv
#include "SolverOptions.h"     // for SolverOptions
/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "siconos_debug.h"             // for DEBUG_EXPR, DEBUG_PRINTF
#include "numerics_verbose.h"  // for verbose
#include "SiconosBlas.h"             // for cblas_dcopy, cblas_daxpy, cblas_dnrm2

int convexQP_compute_error_reduced(
  ConvexQP* problem,
  double *z, double *w, double tolerance,
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
  if(!options->dWork)
  {
    options->dWork = (double*)calloc(2*n,sizeof(double));
  }
  double *ztmp =  options->dWork;
  double *wtmp =  &(options->dWork[n]);

  /* q --> w */
  cblas_dcopy(n, problem->q, 1, w, 1);

  /* M z + q --> w */
  NM_gemv(1.0, problem->M, z, 1.0, w);
  DEBUG_EXPR(NV_display(w,n));
  DEBUG_EXPR(NM_display(problem->M));
  cblas_dcopy(n, z, 1, ztmp, 1);
  cblas_daxpy(n, -1.0, w, 1, ztmp, 1) ;
  DEBUG_EXPR(NV_display(w,n));
  problem->ProjectionOnC(problem,ztmp,wtmp);

  cblas_daxpy(n, -1.0, z, 1, wtmp, 1) ;
  *error = cblas_dnrm2(n, wtmp, incx);

  /* Computes error */
  if(fabs(norm) > DBL_EPSILON)
    *error /= norm;

  DEBUG_PRINTF("error = %e\n",*error);
  if(*error > tolerance)
  {
    numerics_printf_verbose(2," Numerics - convexQP_compute_error: error = %g > tolerance = %g.\n",
                            *error, tolerance);
    return 1;
  }
  else
    return 0;
}



int convexQP_compute_error(
  ConvexQP* problem,
  double *z, double * xi,
  double *w, double * u,
  double tolerance,
  double scaling,
  SolverOptions * options,
  double norm_q,
  double norm_b,
  double * error)
{
  DEBUG_BEGIN("convexQP_compute_error(...)\n");
  assert(problem);
  assert(z);
  assert(u);
  assert(xi);
  assert(error);

  int incx = 1;
  size_t n = problem->size;
  size_t m = problem->m;

  DEBUG_PRINTF("scaling= %12.8e\n", scaling);
  DEBUG_PRINTF("norm_b = %12.8e\n", norm_b);
  DEBUG_PRINTF("norm_q = %12.8e\n", norm_q);

  double norm_xi = cblas_dnrm2(m,xi,1);
  DEBUG_PRINTF("norm of xi %e\n", norm_xi);
  double norm_u = cblas_dnrm2(m,u,1);
  DEBUG_PRINTF("norm of u %e\n", norm_u);
  DEBUG_PRINTF("norm of z  %e\n", cblas_dnrm2(n,z,1));

  double tmp = 0.0;

  *error = 0.;

  if(!options->dWork || options->dWorkSize < 2*m+n)
  {
    options->dWork = (double *)calloc(2*n,sizeof(double));
    options->dWorkSize = 2*m+n;
  }

  double *tmp_m =  options->dWork;
  double *tmp_m1 = &(options->dWork[m]) ;

  double * tmp_n =  &(options->dWork[m+m]);

  /****************************************/
  /* error in Mz + q - rho A^T xi =0      */
  /****************************************/

  /* q --> w */
  cblas_dcopy(n, problem->q, 1, w, 1);

  /* M z + q --> w */
  NM_gemv(1.0, problem->M, z, 0.0, w);
  double norm_Mz =  cblas_dnrm2(n,w,1);
  DEBUG_PRINTF("norm of Mz %e\n", norm_Mz);
  cblas_daxpy(n, 1.0, problem->q, 1, w, 1);

  /* Check that w= A^T xi */
  if(!problem->A)
  {
    cblas_daxpy(n, scaling, xi, 1, tmp_n, 1);
  }
  else
  {
    NM_tgemv(scaling, problem->A, xi, 0.0, tmp_n);
  }
  double norm_rhoATxi =  cblas_dnrm2(n,tmp_n,1);
  DEBUG_PRINTF("norm of rhoATxi %e, ATxi = %e \n", norm_rhoATxi, norm_rhoATxi/scaling);

  cblas_daxpy(n, -1.0, w, 1, tmp_n, 1);
  *error = cblas_dnrm2(n, tmp_n, incx);
  DEBUG_PRINTF("absolute error of Mz + q - rho A^T xi  = %e\n", *error);

  double relative_scaling = fmax(norm_q, fmax(norm_Mz,norm_rhoATxi));
  if(fabs(relative_scaling) > DBL_EPSILON)
    *error = *error/relative_scaling;
  DEBUG_PRINTF("relatice error of Mz + q - rho A^T xi  = %e\n", *error);


  /****************************************/
  /* error in u = A z + b                 */
  /****************************************/


  if(!problem->A)
  {
    cblas_dcopy(n,z,1,tmp_m,1);
  }
  else
  {
    /* A z + b --> u */
    NM_gemv(1.0, problem->A, z, 0.0, tmp_m);
  }
  double norm_Az= cblas_dnrm2(m, tmp_m, 1);
  /* b --> u */
  cblas_daxpy(m, 1.0, problem->b, 1, tmp_m, 1);


  cblas_daxpy(m, -1.0, u, 1, tmp_m, 1) ;
  tmp = cblas_dnrm2(m, tmp_m, incx);
  DEBUG_PRINTF("absolute error of u = A z + b  %e\n", tmp);

  relative_scaling = fmax(norm_b, norm_Az);
  if(fabs(relative_scaling) > DBL_EPSILON)
    tmp = tmp/relative_scaling;

  DEBUG_PRINTF("relative error of u = A z + b  %e\n", tmp);
  *error += tmp ;

  /****************************************/
  /* error in  - xi \in \partial \Psi_C(u)*/
  /****************************************/

  /* Check that - xi \in \partial \Psi_C(u) */
  cblas_dcopy(m, u, 1, tmp_m, 1);
  cblas_daxpy(m, -1.0, xi, 1, tmp_m, 1) ;

  problem->ProjectionOnC(problem,tmp_m,tmp_m1);

  cblas_daxpy(m, -1.0, u, 1, tmp_m1, 1) ;

  tmp= cblas_dnrm2(m, tmp_m1, incx);
  DEBUG_PRINTF("absolute error in complementarity= %e\n", tmp);

  relative_scaling = fmax(norm_u, norm_xi);
  if(fabs(relative_scaling) > DBL_EPSILON)
    tmp = tmp/relative_scaling;
  DEBUG_PRINTF("relative error in complementarity= %e\n", tmp);

  *error += tmp ;
  DEBUG_PRINTF("error = %e\n",*error);

  if(*error > tolerance)
  {
    if(verbose > 1)
      printf(" Numerics - convexQP_compute_error: error = %g > tolerance = %g.\n",
             *error, tolerance);
    DEBUG_END("convexQP_compute_error(...)\n");
    return 1;
  }
  else
  {
    DEBUG_END("convexQP_compute_error(...)\n");
    return 0;
  }
}
