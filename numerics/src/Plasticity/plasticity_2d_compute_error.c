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
#include <stdlib.h>  // for malloc

#include "NumericsMatrix.h"     // for NM_prod_mv_3x3, NM_gemv
#include "PlasticityProblem.h"  // for Mohrcoulomb2dproblem
#include "Plasticity_options.h"
#include "SOC_scaling.h"
#include "SiconosBlas.h"  // for cblas_dcopy, cblas_dnrm2
#include "SolverOptions.h"
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

double plasticity_2d_scaling(int norm_id, int reduction, const double* r, const double* v,
                             int m, const double* eta, double alpha) {
  int mean = (reduction == PLASTICITY_SCALING_REDUCTION_MEAN);
  int mixed = (reduction == PLASTICITY_SCALING_REDUCTION_MIXED);
  double s = 1.0;
  switch (norm_id) {
    case PLASTICITY_SCALING_EUCLIDEAN:
      s = mean    ? SOC_geom_scaling_euclidean(r, v, m, 3)
          : mixed ? SOC_mixed_scaling_euclidean(r, v, m, 3)
                  : SOC_scaling_euclidean(r, v, m, 3);
      break;
    case PLASTICITY_SCALING_CONE:
      s = mean    ? SOC_geom_scaling_cone(r, v, m, 3)
          : mixed ? SOC_mixed_scaling_cone(r, v, m, 3)
                  : SOC_scaling_cone(r, v, m, 3);
      break;
    case PLASTICITY_SCALING_CONE_ETA:
      s = mean    ? SOC_geom_scaling_cone_mu(r, v, m, eta, 3)
          : mixed ? SOC_mixed_scaling_cone_mu(r, v, m, eta, 3)
                  : SOC_scaling_cone_mu(r, v, m, eta, 3);
      break;
    case PLASTICITY_SCALING_CONE_WEIGHTED:
      s = mean    ? SOC_geom_scaling_cone_weighted(r, v, m, 3, alpha)
          : mixed ? SOC_mixed_scaling_cone_weighted(r, v, m, 3, alpha)
                  : SOC_scaling_cone_weighted(r, v, m, 3, alpha);
      break;
    default:
      break;
  }
  /* safeguard, as in fc3d (the mixed scales are not safeguarded by SOC_scaling) */
  return (s > 1.0) ? s : 1.0;
}

double plasticity_2d_DruckerPrager_local_scaled_error(const double M[9], const double q[3],
                                                      double eta, double theta,
                                                      const double r[3]) {
  double w[3], v[3], P[3];
  for (int i = 0; i < 3; i++) w[i] = q[i] + M[i] * r[0] + M[i + 3] * r[1] + M[i + 6] * r[2];
  double nM = 0.0; /* Frobenius norm: cheap, within a factor sqrt(3) of ||M||_2 */
  for (int i = 0; i < 9; i++) nM += M[i] * M[i];
  nM = sqrt(nM);
  double rho = (nM > 0.0) ? 1.0 / nM : 1.0;
  v[0] = rho * (w[0] + theta * sqrt(w[1] * w[1] + w[2] * w[2]));
  v[1] = rho * w[1];
  v[2] = rho * w[2];
  for (int i = 0; i < 3; i++) P[i] = r[i] - v[i];
  projectionOnCone(P, eta);
  double F = sqrt((r[0] - P[0]) * (r[0] - P[0]) + (r[1] - P[1]) * (r[1] - P[1]) +
                  (r[2] - P[2]) * (r[2] - P[2]));
  double nr = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
  double nv = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  return F / fmax(1.0, nr + nv);
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

  /* Scaling, as in fc3d_compute_error (options->iparam[PLASTICITY_SCALING_*]): the residual
   * of each element (PER_CONE) or of the whole vector (PER_VECTOR) is divided by a scale
   * combining the norms of r and v (each safeguarded by max(1, .)). */
  int norm_id = PLASTICITY_SCALING_EUCLIDEAN;
  int scope = PLASTICITY_SCALING_SCOPE_PER_CONE;
  /* defaults (options == NULL): as in plasticity_2d_nsgs_set_default */
  int reduction = PLASTICITY_SCALING_REDUCTION_MIXED;
  double alpha = 1.0;
  int rho_mode = PLASTICITY_SCALING_RHO_ELEMENT;
  if (options) {
    rho_mode = options->iparam[PLASTICITY_SCALING_RHO];
    norm_id = options->iparam[PLASTICITY_SCALING_NORM];
    scope = options->iparam[PLASTICITY_SCALING_SCOPE];
    reduction = options->iparam[PLASTICITY_SCALING_REDUCTION];
    alpha = options->dparam[PLASTICITY_SCALING_ALPHA];
  }
  /* rho of the natural map (balance of the stress and of the strain rate): from the norms of
   * the diagonal blocks T^T M_ee T in the variables of the solvers */
  double* rho = NULL;
  if (rho_mode != PLASTICITY_SCALING_RHO_ONE) {
    rho = (double*)malloc((size_t)nc * sizeof(double));
    if (!rho) return numerics_error("plasticity_2d_compute_error", "bad alloc");
    double T[9], Tinv[9];
    int change = plasticity_2d_change_of_variables(problem, T, Tinv);
    double rho_min = INFINITY;
    for (int ic = 0; ic < nc; ic++) {
      double buffer[9], *B = buffer;
      NM_extract_diag_block3(problem->M, ic, &B);
      double By[9]; /* column-major T^T B T, T row-major */
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
          double sum = 0.0;
          if (change) {
            for (int k = 0; k < 3; k++)
              for (int l = 0; l < 3; l++) sum += T[3 * k + i] * B[k + 3 * l] * T[3 * l + j];
          } else
            sum = B[i + 3 * j];
          By[i + 3 * j] = sum;
        }
      double nB = 0.0; /* Frobenius norm (as in the local error of the local solvers) */
      for (int i = 0; i < 9; i++) nB += By[i] * By[i];
      nB = sqrt(nB);
      rho[ic] = (nB > 0.0) ? 1.0 / nB : 1.0;
      rho_min = fmin(rho_min, rho[ic]);
    }
    if (rho_mode == PLASTICITY_SCALING_RHO_GLOBAL)
      for (int ic = 0; ic < nc; ic++) rho[ic] = rho_min;
  }

  /* Per element, in the variables of the solvers (cone coordinates, cohesion shift, plane
   * stress scaling): r the stress, v the modified strain rate (w + theta ||w_T|| e_N for
   * Drucker-Prager, w for Von Mises), multiplied by rho, and the natural-map residual
   * F = r - P_K(r - rho v). */
  double* rv = (double*)malloc((size_t)(6 * nc + nc) * sizeof(double));
  if (!rv) {
    free(rho);
    return numerics_error("plasticity_2d_compute_error", "bad alloc");
  }
  double* r_all = rv;          /* 3 nc */
  double* v_all = rv + 3 * nc; /* 3 nc */
  double* f2 = rv + 6 * nc;    /* ||F_e||^2, nc */
  double* eta_all = NULL;      /* cone coefficients for the CONE_ETA norm */
  PlasticityStressCoordinates coordinates = problem->stress_coordinates;
  for (int ic = 0, ic3 = 0; ic < nc; ic++, ic3 += 3) {
    double* r = r_all + ic3;
    double* v = v_all + ic3;
    double P[3];
    plasticity_2d_stress_to_cone(coordinates, 1, z + ic3, r);
    plasticity_2d_strainrate_to_cone(coordinates, 1, w + ic3, v);
    if (drucker_prager) {
      /* with a cohesion c, the error is computed for sigma + c/eta e_N (cone without
       * cohesion, see Plasticity_DruckerPrager_model) */
      const double* c = problem->model.drucker_prager->cohesion;
      double eta = problem->model.drucker_prager->eta[ic];
      if (c && c[ic] != 0.0) r[0] += c[ic] / eta;
      v[0] += problem->model.drucker_prager->theta[ic] * sqrt(v[1] * v[1] + v[2] * v[2]);
      if (rho)
        for (int i = 0; i < 3; i++) v[i] *= rho[ic];
      for (int i = 0; i < 3; i++) P[i] = r[i] - v[i];
      projectionOnCone(P, eta);
    } else {
      double sigma_y = problem->model.von_mises->sigma_y[ic];
      PlasticityVonMisesHypothesis hypothesis = problem->model.von_mises->hypothesis;
      if (hypothesis == PLASTICITY_VONMISES_PLANE_STRESS) {
        /* error in the scaled variables y = D^{-1} sigma, w_y = D w,
         * D = diag(sqrt2, sqrt(2/3), sqrt(2/3)): the ellipsoid becomes the ball */
        const double cc = sqrt(2.0 / 3.0);
        r[0] *= M_SQRT1_2;
        r[1] /= cc;
        r[2] /= cc;
        v[0] *= M_SQRT2;
        v[1] *= cc;
        v[2] *= cc;
      }
      if (rho)
        for (int i = 0; i < 3; i++) v[i] *= rho[ic];
      for (int i = 0; i < 3; i++) P[i] = r[i] - v[i];
      if (hypothesis == PLASTICITY_VONMISES_PLANE_STRAIN)
        plasticity_2d_projectionOnVonMises(P, sigma_y);
      else
        plasticity_2d_projectionOnVonMisesBall(P, sigma_y);
    }
    f2[ic] = (r[0] - P[0]) * (r[0] - P[0]) + (r[1] - P[1]) * (r[1] - P[1]) +
             (r[2] - P[2]) * (r[2] - P[2]);
  }

  if (norm_id == PLASTICITY_SCALING_CONE_ETA) {
    eta_all = (double*)malloc((size_t)nc * sizeof(double));
    for (int ic = 0; ic < nc; ic++)
      eta_all[ic] = drucker_prager ? problem->model.drucker_prager->eta[ic] : 1.0;
  }
  *error = 0.;
  if (norm_id == PLASTICITY_SCALING_NO) {
    for (int ic = 0; ic < nc; ic++) *error += f2[ic];
    *error = sqrt(*error);
  } else if (scope == PLASTICITY_SCALING_SCOPE_PER_VECTOR) {
    double s = plasticity_2d_scaling(norm_id, reduction, r_all, v_all, nc, eta_all, alpha);
    for (int ic = 0; ic < nc; ic++) *error += f2[ic];
    *error = sqrt(*error) / s;
  } else {
    for (int ic = 0; ic < nc; ic++) {
      double s = plasticity_2d_scaling(norm_id, reduction, r_all + 3 * ic, v_all + 3 * ic, 1,
                                       eta_all ? eta_all + ic : NULL, alpha);
      *error += f2[ic] / (s * s);
    }
    *error = sqrt(*error);
  }
  free(eta_all);
  free(rho);
  free(rv);
  DEBUG_PRINTF("scaled error in complementarity = %12.8e\n", *error);
  DEBUG_END("plasticity_2d_compute_error(...)\n");
  if (*error > tolerance) return 1;

  return 0;
}
