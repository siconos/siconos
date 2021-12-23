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

#include "grfc3d_compute_error.h"
#include <assert.h>                               // for assert
#include <float.h>                                // for DBL_EPSILON
#include <math.h>                                 // for sqrt, fabs
#include <stdlib.h>                               // for NULL, calloc
#include "SiconosBlas.h"                          // for cblas_dcopy
#include "sanitizer.h"                            // for cblas_dcopy_msan
#include "NumericsMatrix.h"                       // for NM_gemv, NM_tgemv, Numeric...
#include "GlobalRollingFrictionContactProblem.h"  // for GlobalRollingFrictionContactPro...
#include "siconos_debug.h"                        // for DEBUG_EXPR, DEBUG_PRINTF
#include "projectionOnRollingCone.h"              // for projectionOnRollingCone
#include "numerics_verbose.h"                     // for numerics_error, numerics_w...
#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"
#endif

#define MIN_RELATIVE_SCALING sqrt(DBL_EPSILON)

void grfc3d_unitary_compute_and_add_error(
  double* r,
  double* u,
  double mu,
  double mur,
  double*  error,
  double * worktmp,
  int problemIsNotConvex)
{
  DEBUG_BEGIN("grfc3d_unitary_compute_and_add_error(...)\n");

  if (problemIsNotConvex == 0)
  {
    worktmp[0] = r[0] - u[0];
  }
  else
    worktmp[0] = r[0] - u[0]
                      - mu  * sqrt(u[1] * u[1] + u[2] * u[2])
                      - mur * sqrt(u[3] * u[3] + u[4] * u[4]);

  worktmp[1] = r[1] -  u[1] ;
  worktmp[2] = r[2] -  u[2] ;
  worktmp[3] = r[3] -  u[3] ;
  worktmp[4] = r[4] -  u[4] ;

  projectionOnRollingCone(worktmp, mu, mur);

  worktmp[0] = r[0] -  worktmp[0];
  worktmp[1] = r[1] -  worktmp[1];
  worktmp[2] = r[2] -  worktmp[2];
  worktmp[3] = r[3] -  worktmp[3];
  worktmp[4] = r[4] -  worktmp[4];

  *error +=
    worktmp[0] * worktmp[0] +
    worktmp[1] * worktmp[1] +
    worktmp[2] * worktmp[2] +
    worktmp[3] * worktmp[3] +
    worktmp[4] * worktmp[4];

  DEBUG_END("grfc3d_unitary_compute_and_add_error(...)\n");
}

int grfc3d_compute_error(GlobalRollingFrictionContactProblem* problem,
                               double*  reaction, double*  velocity,
                               double*  globalVelocity,
                               double tolerance,
                               double* error, int problemIsNotConvex)
{
  DEBUG_BEGIN("grfc3d_compute_error(...)\n");
  /* Checks inputs */
  if(problem == NULL || globalVelocity == NULL || velocity == NULL || reaction == NULL || error == NULL)
    numerics_error("grfc3d_compute_error", "null input");

  int incx = 1, incy = 1;
  int nc = problem->numberOfContacts;
  assert(nc > 0);
  int nd = nc * 5;
  size_t m = problem->M->size0;
  double *mu = problem->mu;
  double *mur = problem->mu_r;

  NumericsMatrix *H = problem->H;
  NumericsMatrix *M = problem->M;

  double norm_r = cblas_dnrm2(nd,reaction,1);
  double norm_u = cblas_dnrm2(m,velocity,1);
  double norm_q = cblas_dnrm2(m, problem->q, 1);        // = f
  double norm_b = cblas_dnrm2(nd, problem->b, 1);       // = w

  /* --- Relative primal residual = |-Mv + Hr + q|/max{|Mv|, |Hr|, |q|} --- */
  double* tmp_m_1 = (double *)calloc(m,sizeof(double));
  double* tmp_m_2 = (double *)calloc(m,sizeof(double));


  cblas_dcopy_msan(m, problem->q, 1, tmp_m_2, 1);        // tmp_m_2 = q
  // if(tmp_m_2)
  // {
  //   printf("\nDEBUG q vector:\n");
  //   for(unsigned int i = 0; i < m; i++) printf("q[ %i ] = %12.8e\n", i, tmp_m_2[i]);
  // }

  NM_gemv(1.0, H, reaction, 0.0, tmp_m_1);        // tmp_m_1 = Hr
  double norm_Hr = cblas_dnrm2(m,tmp_m_1,1);
  cblas_daxpy(m, 1.0, tmp_m_1, 1, tmp_m_2, 1);      // tmp_m_2 = Hr + q
  NM_gemv(-1.0, M, globalVelocity, 0.0, tmp_m_1); // tmp_m_1 = -Mv
  double norm_Mv = cblas_dnrm2(m,tmp_m_1,1);
  cblas_daxpy(m, 1.0, tmp_m_1, 1, tmp_m_2, 1);      // tmp_m_2 = -Mv + Hr + q

  // if(tmp_m_2)
  // {
  //   printf("\nDEBUG error_primal:\n");
  //   for(unsigned int i = 0; i < m; i++) printf("error_primal[ %i ] = %12.16e\n", i, tmp_m_2[i]);
  // }

  double error_primal = cblas_dnrm2(m,tmp_m_2,1);   // error_primal = |-Mv + Hr - q|
// printf("\n\n#################### 001 error = %12.16e ####################\n", *error);
  double relative_scaling = fmax(norm_q, fmax(norm_Mv, norm_Hr));
  if(relative_scaling > MIN_RELATIVE_SCALING)
    *error = error_primal/relative_scaling;         // error = |-Mv + Hr - q|/max{|Mv|, |Hr|, |q|}
  else
    *error = error_primal;
// printf("#################### 002 error = %12.16e ####################\n", *error);
  free(tmp_m_1);
  free(tmp_m_2);


  /* --- Relative dual residual = |H'v + b - u|/max{|H'v|, |b|, |u|} --- */
  double* tmp_nd_1 = (double *)calloc(nd,sizeof(double));
  double* tmp_nd_2 = (double *)calloc(nd,sizeof(double));

  cblas_dcopy_msan(nd, problem->b, 1, tmp_nd_2, 1);        // tmp_nd_2 = b
  NM_tgemv(1.0, H, globalVelocity, 0.0, tmp_nd_1);        // tmp_nd_1 = H'v
  double norm_HTv = cblas_dnrm2(nd,tmp_nd_1,1);
  cblas_daxpy(nd, 1.0, tmp_nd_1, 1, tmp_nd_2, 1);      // tmp_nd_2 = H'v + b
  cblas_daxpy(nd, -1.0, velocity, 1, tmp_nd_2, 1);      // tmp_nd_2 = H'v + b - u
  double error_dual = cblas_dnrm2(nd,tmp_nd_2,1);       // error_dual = |H'v + b - u|

  relative_scaling = fmax(norm_u, fmax(norm_b, norm_HTv));
  if(relative_scaling > MIN_RELATIVE_SCALING)
    *error += error_dual/relative_scaling;         // error = |H'v + b - u|/max{|H'v|, |b|, |u|}
  else
    *error += error_dual;
// printf("#################### 003 error = %12.16e ####################\n", *error);
  free(tmp_nd_1);
  free(tmp_nd_2);


  /* --- Projection error = |r - projectionOnRollingCone(r-u)|/max{|r|, |u|} for convex case --- */
  /* --- Projection error = |r - projectionOnRollingCone(r-u-mu*|uT|-mur*||wR)|/max{|r|, |u|} for non-convex case --- */
  double error_complementarity = 0.0;
  double worktmp[5];
  for(int ic = 0 ; ic < nc ; ic++)
  {
    grfc3d_unitary_compute_and_add_error(&reaction[ic * 5], &velocity[ic * 5],
      mu[ic], mur[ic] , &error_complementarity,  worktmp, problemIsNotConvex);
  }

  error_complementarity = sqrt(error_complementarity);

  relative_scaling = fmax(norm_u, norm_r);
  if(relative_scaling > MIN_RELATIVE_SCALING)
    *error += error_complementarity/relative_scaling;
  else
    *error += error_complementarity;

// printf("#################### 004 error = %12.16e ####################\n\n", *error);
  //numerics_printf_verbose(1,"---- GRFC3D - Compute Error ");
  DEBUG_END("grfc3d_compute_error(...)\n");

  if(*error > tolerance)
    return 1;

  return 0;
}
