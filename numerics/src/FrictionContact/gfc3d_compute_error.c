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

#include "gfc3d_compute_error.h"
#include <float.h>                         // for DBL_EPSILON
#include <math.h>                          // for fabs, sqrt
#include <stdlib.h>                        // for NULL, calloc
#include "GlobalFrictionContactProblem.h"  // for GlobalFrictionContactProblem
#include "NumericsMatrix.h"                // for NM_gemv, NM_tgemv, Numeric...
#include "SolverOptions.h"                 // for SolverOptions
#include "fc3d_compute_error.h"            // for fc3d_unitary_compute_and_a...
#include "numerics_verbose.h"              // for numerics_error, numerics_w...
#include "sanitizer.h"                     // for cblas_dcopy_msan
#include "SiconosBlas.h"                         // for cblas_dcopy, cblas_dnrm2
#include "projectionOnCone.h"

/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "siconos_debug.h"

#define MIN_RELATIVE_SCALING sqrt(DBL_EPSILON)

int gfc3d_compute_error(GlobalFrictionContactProblem* problem,
                        double*  reaction, double*  velocity,
                        double*  globalVelocity,
                        double tolerance,
                        SolverOptions * options,
                        double norm_q, double norm_b,
                        double* restrict error)

{
  DEBUG_BEGIN("gfc3d_compute_error(...)\n");
  /* Checks inputs */
  if(problem == NULL || globalVelocity == NULL)
    numerics_error("gfc3d_compute_error", "null input");

  /* Computes error = dnorm2( GlobalVelocity -M^-1( q + H reaction)*/
  int nc = problem->numberOfContacts;
  int m = nc * 3;
  size_t n = problem->M->size0;
  double *mu = problem->mu;
  double *q = problem->q;

  DEBUG_PRINTF("norm_b = %12.8e\n", norm_b);
  DEBUG_PRINTF("norm_q = %12.8e\n", norm_q);
  double norm_r= cblas_dnrm2(m,reaction,1);
  DEBUG_PRINTF("norm of reaction %e\n", cblas_dnrm2(m,reaction,1));
  DEBUG_PRINTF("norm of global velocity  %e\n", cblas_dnrm2(n,globalVelocity,1));

  /* DEBUG_EXPR(NV_display(globalVelocity,n)); */
  /* DEBUG_EXPR(NV_display(reaction,m)); */
  /* DEBUG_EXPR(NV_display(velocity,m)); */

  NumericsMatrix *H = problem->H;
  NumericsMatrix *M = problem->M;

  if(!options->dWork || options->dWorkSize < 2*n)
  {
    options->dWork = (double *)calloc(2*n,sizeof(double));
    options->dWorkSize = 2*n;
  }
  double* tmp = options->dWork;
  double* tmp_1 = &options->dWork[n];

  cblas_dcopy_msan(n, q, 1, tmp_1, 1);
  if(nc >0)
  {
    NM_gemv(1.0, H, reaction, 0.0, tmp);
  }
  double norm_Hr = cblas_dnrm2(n,tmp,1);
  DEBUG_PRINTF("norm of H r %e\n", norm_Hr);


  cblas_daxpy(n, 1.0, tmp, 1, tmp_1, 1);


  NM_gemv(-1.0, M, globalVelocity, 0.0, tmp);
  double norm_Mv = cblas_dnrm2(n,tmp,1);
  DEBUG_PRINTF("norm of M v %e\n", norm_Mv);

  cblas_daxpy(n, 1.0, tmp, 1, tmp_1, 1);

  double relative_scaling = fmax(norm_q, fmax(norm_Mv,norm_Hr));
  *error = cblas_dnrm2(n,tmp_1,1);
  DEBUG_PRINTF("absolute error  of -M v + H R + q = %e\n", *error);
  if(fabs(relative_scaling) > MIN_RELATIVE_SCALING)
    *error = *error/relative_scaling;

  DEBUG_PRINTF("relative error  of -M v + H R + q = %e\n", *error);

  /* CHECK_RETURN(!NM_gesv_expert(problem->M, globalVelocity, NM_KEEP_FACTORS)); */

  double error_complementarity =0.0;

  /* Checks inputs */
  if(reaction == NULL || velocity == NULL)
    numerics_error("gfc3d_compute_error", "null input");


  /* we re-compute local velocity */
  /* the error in the equation u = H^T v +b is then accuaret to machine precision */

  cblas_dcopy(m, problem->b, 1, velocity, 1);
  NM_tgemv(1, H, globalVelocity, 1.0, velocity);
  double norm_u = cblas_dnrm2(m,velocity,1);
  DEBUG_PRINTF("norm of velocity %e\n", norm_u);

  double worktmp[3];
  for(int ic = 0 ; ic < nc ; ic++)
  {
    fc3d_unitary_compute_and_add_error(&reaction[ic * 3], &velocity[ic * 3], mu[ic],
                                       &error_complementarity,  worktmp);
  }

  error_complementarity = sqrt(error_complementarity);

  DEBUG_PRINTF("absolute error in complementarity= %e\n", error_complementarity);

  relative_scaling = fmax(norm_u, norm_r);
  if(fabs(relative_scaling) > MIN_RELATIVE_SCALING)
    error_complementarity = error_complementarity/relative_scaling;
  DEBUG_PRINTF("relative error in complementarity= %e\n", error_complementarity);


  *error += error_complementarity;
  DEBUG_PRINTF("relative error = %e\n", *error);
  if(*error > tolerance)
  {
    DEBUG_END("gfc3d_compute_error(...)\n");
    return 1;
  }
  else
  {
    DEBUG_END("gfc3d_compute_error(...)\n");
    return 0;
  }

}
static inline  void fc3d_unitary_compute_and_add_error_convex(double* restrict r, double* restrict u, double mu, double* restrict error, double * worktmp)
{

  //double normUT;
  //double worktmp[3];
  /* Compute the modified local velocity */
  /* worktmp[0] = r[0] - u[0] - mu *  hypot(u[1], u[2]); */
  worktmp[0] = r[0] -  u[0] ;
  worktmp[1] = r[1] -  u[1] ;
  worktmp[2] = r[2] -  u[2] ;
  projectionOnCone(worktmp, mu);
  worktmp[0] = r[0] -  worktmp[0];
  worktmp[1] = r[1] -  worktmp[1];
  worktmp[2] = r[2] -  worktmp[2];
  *error +=  worktmp[0] * worktmp[0] + worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2];

}
int gfc3d_compute_error_convex(GlobalFrictionContactProblem* problem,
                               double*  reaction, double*  velocity,
                               double*  globalVelocity,
                               double tolerance,
                               SolverOptions * options,
                               double norm_q, double norm_b,
                               double* restrict error)

{
  DEBUG_BEGIN("gfc3d_compute_error_convex(...)\n");
  /* Checks inputs */
  if(problem == NULL || globalVelocity == NULL)
    numerics_error("gfc3d_compute_error", "null input");

  /* Computes error = dnorm2( GlobalVelocity -M^-1( q + H reaction)*/
  int nc = problem->numberOfContacts;
  int m = nc * 3;
  size_t n = problem->M->size0;
  double *mu = problem->mu;
  double *q = problem->q;

  DEBUG_PRINTF("norm_b = %12.8e\n", norm_b);
  DEBUG_PRINTF("norm_q = %12.8e\n", norm_q);
  double norm_r= cblas_dnrm2(m,reaction,1);
  DEBUG_PRINTF("norm of reaction %e\n", cblas_dnrm2(m,reaction,1));
  DEBUG_PRINTF("norm of global velocity  %e\n", cblas_dnrm2(n,globalVelocity,1));

  /* DEBUG_EXPR(NV_display(globalVelocity,n)); */
  /* DEBUG_EXPR(NV_display(reaction,m)); */
  /* DEBUG_EXPR(NV_display(velocity,m)); */

  NumericsMatrix *H = problem->H;
  NumericsMatrix *M = problem->M;

  int d_work_size= 2*n > m ? 2*n : m;
  if(!options->dWork || options->dWorkSize < d_work_size)
  {
    if (options->dWork)
      free(options->dWork);
    options->dWork = (double *)calloc(d_work_size,sizeof(double));
    options->dWorkSize = d_work_size;
  }

  double* tmp = options->dWork;
  double* tmp_1 = &options->dWork[n];

  cblas_dcopy_msan(n, q, 1, tmp_1, 1);
  if(nc >0)
  {
    NM_gemv(1.0, H, reaction, 0.0, tmp);
  }
  double norm_Hr = cblas_dnrm2(n,tmp,1);
  DEBUG_PRINTF("norm of H r %e\n", norm_Hr);

  cblas_daxpy(n, 1.0, tmp, 1, tmp_1, 1);

  NM_gemv(-1.0, M, globalVelocity, 0.0, tmp);
  double norm_Mv = cblas_dnrm2(n,tmp,1);
  DEBUG_PRINTF("norm of M v %e\n", norm_Mv);

  cblas_daxpy(n, 1.0, tmp, 1, tmp_1, 1);

  double relative_scaling = fmax(norm_q, fmax(norm_Mv,norm_Hr));
  *error = cblas_dnrm2(n,tmp_1,1);
  DEBUG_PRINTF("absolute error  of -M v + H R + q = %e\n", *error);
  if(fabs(relative_scaling) > MIN_RELATIVE_SCALING)
    *error = *error/relative_scaling;

  DEBUG_PRINTF("relative error  of -M v + H R + q = %e\n", *error);
  //printf("relative error  of -M v + H R + q = %e\n", *error);
  /* CHECK_RETURN(!NM_gesv_expert(problem->M, globalVelocity, NM_KEEP_FACTORS)); */

  double error_complementarity =0.0;

  /* Checks inputs */
  if(reaction == NULL || velocity == NULL)
    numerics_error("gfc3d_compute_error", "null input");


  /* we re-compute local velocity */
  /* the error in the equation u = H^T v +b is then accuaret to machine precision */

  /* cblas_dcopy(m, problem->b, 1, velocity, 1); */
  /* NM_tgemv(1, H, globalVelocity, 1.0, velocity); */
  /* double norm_u = cblas_dnrm2(m,velocity,1); */
  /* DEBUG_PRINTF("norm of velocity %e\n", norm_u); */

  /* computation of the relative feasibility error = relative norm of the primal residual*/
  /* |H*v+w-u|/max{|H*v|, |w|, |u|}*/
  double *primal_residual = tmp;
  NM_tgemv(1, H, globalVelocity, 0.0, primal_residual);
  double norm_Hv = cblas_dnrm2(m, primal_residual, 1);
  DEBUG_PRINTF(" norm_Hv  =  %e\n", norm_Hv);

  cblas_daxpy(m, 1.0, problem->b, 1, primal_residual, 1);
  cblas_daxpy(m, -1.0, velocity, 1, primal_residual, 1);
  double norm_w = cblas_dnrm2(m, problem->b, 1);
  double norm_u = cblas_dnrm2(m, velocity, 1);
  DEBUG_PRINTF(" norm_u  = %e, norm_b  = %e\n", norm_u, norm_b);


  relative_scaling = fmax(norm_Hv, fmax(norm_w, norm_u));
  double norm_primal_residual = cblas_dnrm2(m, primal_residual, 1);

  DEBUG_PRINTF(" norm primal_residual  = %e\n", norm_primal_residual);
  if(fabs(relative_scaling) > MIN_RELATIVE_SCALING)
    *error += norm_primal_residual/relative_scaling;
  else
    *error += norm_primal_residual;

  DEBUG_PRINTF(" error with primal_residual  = %e\n", *error);

  double worktmp[3];
  for(int ic = 0 ; ic < nc ; ic++)
  {
    fc3d_unitary_compute_and_add_error_convex(&reaction[ic * 3], &velocity[ic * 3], mu[ic],
        &error_complementarity,  worktmp);
  }

  error_complementarity = sqrt(error_complementarity);

  //  error_complementarity = cblas_ddot(m, velocity, 1, reaction, 1)/m;

  DEBUG_PRINTF("absolute error in complementarity= %e\n", error_complementarity);

  relative_scaling = fmax(norm_u, norm_r);
  //relative_scaling = fmax(1.0, relative_scaling);
  if(fabs(relative_scaling) > MIN_RELATIVE_SCALING)
    error_complementarity = error_complementarity/relative_scaling;
  DEBUG_PRINTF("relative error in complementarity= %e\n", error_complementarity);

  *error += error_complementarity;

  DEBUG_PRINTF("relative error = %e\n", *error);
  numerics_printf_verbose(1,"---- GFC3D - Compute Error Convex case");
  if(*error > tolerance)

    {
    DEBUG_END("gfc3d_compute_error_convex(...)\n");
    return 1;
  }
  else
  {
    DEBUG_END("gfc3d_compute_error_convex(...)\n");
    return 0;
  }

}
