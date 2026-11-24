// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "plasticity_2d_compute_error.h"

#include <assert.h>  // for assert
#include <float.h>   // for DBL_EPSILON
#include <math.h>    // for sqrt, fabs
#include <stddef.h>  // for NULL

#include "NumericsMatrix.h"     // for NM_prod_mv_3x3, NM_gemv
#include "PlasticityProblem.h"  // for Mohrcoulomb2dproblem
#include "SiconosBlas.h"        // for cblas_dcopy, cblas_dnrm2
#include "numerics_errors.h"
#include "numerics_verbose.h"
#include "plasticity_2d_stress_coordinates.h"
#include "plasticity_2d_vonmises.h"
#include "projectionOnCone.h"  // for projectionOnCone
/* #define DEBUG_NOCOLOR */
/* #define DEBUG_STDOUT */
/* #define DEBUG_MESSAGES */
#include "siconos_debug.h"  // for DEBUG_PRINTF, DEBUG_EXPR, DEBUG_...
#ifdef DEBUG_MESSAGES
#include "NumericsVector.h"
#endif

void plasticity_2d_unitary_compute_and_add_error(double* restrict r, double* restrict u,
                                                 double eta, double theta,
                                                 double* restrict error, double* worktmp) {
  // double normUT;
  // double worktmp[3];
  /* Compute the modified local velocity */
  /* worktmp[0] = r[0] - u[0] - eta *  hypot(u[1], u[2]); */
  worktmp[0] = r[0] - u[0] - theta * sqrt(u[1] * u[1] + u[2] * u[2]);
  worktmp[1] = r[1] - u[1];
  worktmp[2] = r[2] - u[2];
  projectionOnCone(worktmp, eta);
  worktmp[0] = r[0] - worktmp[0];
  worktmp[1] = r[1] - worktmp[1];
  worktmp[2] = r[2] - worktmp[2];
  *error += worktmp[0] * worktmp[0] + worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2];
}

int plasticity_2d_DruckerPrager_zero_is_solution(const double q[3], double eta, double theta) {
  double nT = sqrt(q[1] * q[1] + q[2] * q[2]);
  double scale = fmax(fabs(q[0]), nT);
  return (eta * nT <= q[0] + theta * nT + DBL_EPSILON * scale) ? 1 : 0;
}

int plasticity_2d_compute_error(PlasticityProblem* problem, double* z, double* w,
                                double tolerance, SolverOptions* options, double norm,
                                double* error) {
  DEBUG_BEGIN("plasticity_2d_compute_error(...)\n");
  CHECK_NULL(problem);
  CHECK_NULL(z);
  CHECK_NULL(w);
  CHECK_NULL(error);

  if (problem->model_type != PLASTICITY_MODEL_DRUCKER_PRAGER &&
      problem->model_type != PLASTICITY_MODEL_VON_MISES) {
    *error = 0.0;
    DEBUG_END("plasticity_2d_compute_error(...)\n");
    return numerics_error("plasticity_2d_compute_error", "unknown plasticity model %d",
                          problem->model_type);
  }

  /* Computes w = Mz + q */
  int incx = 1, incy = 1;
  int nc = problem->numberOfCones;
  int n = nc * 3;
  int drucker_prager = (problem->model_type == PLASTICITY_MODEL_DRUCKER_PRAGER);

  /* Compute the current velocity */
  cblas_dcopy(n, problem->q, incx, w, incy);  // w <-q
  NM_prod_mv_3x3(n, n, problem->M, z, w);     // w = Mz +q

  DEBUG_PRINTF("norm of the reaction %e\n", cblas_dnrm2(n, z, 1));
  DEBUG_PRINTF("norm of the velocity %e\n", cblas_dnrm2(n, w, 1));
  DEBUG_PRINTF("norm of q = %12.8e\n", norm);
  /* DEBUG_EXPR(NV_display(problem->q,n);); */
  /* DEBUG_EXPR(NV_display(w,n);); */
  /* DEBUG_EXPR(NV_display(z,n);); */

  *error = 0.;
  int ic, ic3;
  double worktmp[3];
  PlasticityStressCoordinates coordinates = problem->stress_coordinates;
  for (ic = 0, ic3 = 0; ic < nc; ic++, ic3 += 3) {
    /* the error is computed in cone coordinates (sigma_c = P s, w_c = P^{-T} w_s) */
    double r[3], v[3];
    plasticity_2d_stress_to_cone(coordinates, 1, z + ic3, r);
    plasticity_2d_strainrate_to_cone(coordinates, 1, w + ic3, v);
    if (drucker_prager) {
      /* with a cohesion c, the error is computed for sigma + c/eta e_N (cone without
       * cohesion, see Plasticity_DruckerPrager_model) */
      const double* c = problem->model.drucker_prager->cohesion;
      double eta = problem->model.drucker_prager->eta[ic];
      if (c && c[ic] != 0.0) r[0] += c[ic] / eta;
      plasticity_2d_unitary_compute_and_add_error(
          r, v, eta, problem->model.drucker_prager->theta[ic], error, worktmp);
    } else
      plasticity_2d_vonmises_unitary_compute_and_add_error(
          r, v, problem->model.von_mises->sigma_y[ic], error, worktmp);
    /*DEBUG_PRINTF("absolute error = %12.8e contact =%i nc= %i\n", *error, ic, nc);*/
  }
  *error = sqrt(*error);
  DEBUG_PRINTF("absolute error in complementarity = %12.8e\n", *error);

  /* Compute relative error */
  // double norm_r = cblas_dnrm2(n, z, 1);
  // double norm_u = cblas_dnrm2(n, w, 1);
  //  double relative_scaling = fmax(norm, fmax(norm_r, norm_u));
  /* double relative_scaling = fmax(norm_r,norm_w); */
  /* double relative_scaling = norm; */

  // if (fabs(relative_scaling) > DBL_EPSILON) *error /= relative_scaling;

  DEBUG_PRINTF("relative error in complementarity = %12.8e\n", *error);
  DEBUG_END("plasticity_2d_compute_error(...)\n");
  if (*error > tolerance) return 1;

  return 0;
}
