// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*!\file plasticity_2d_nsn.c
 * \brief Global semi-smooth Newton method for plasticity problems (PLASTICITY_2D_NSN).
 *
 * Problem (cone coordinates, as given by plasticity_2d_driver after its change of
 * variables): find sigma with, for each element e,
 *
 *   Phi_e(sigma) = sigma_e - P_{K_e}(sigma_e - rho_e v_e(w)) = 0,   w = M sigma + q,
 *
 * v_e = w_e (Von Mises), v_e = w_e + theta_e ||w_eT|| e_N (Drucker-Prager, dilatancy),
 * rho_e = c ||M_ee^{-1}||_2 (default; c = 10 for the full Newton steps, 1 for the merit
 * phase) or c/||M_ee||_2. An element of the generalized
 * Jacobian of Phi is
 *
 *   H = I - J_P (I - R G M) = (I - J_P) + J_P R G M,
 *
 * J_P = blockdiag(J_P_e) (Jacobian of the projection at y_e = sigma_e - rho_e v_e),
 * R = blockdiag(rho_e I), G = blockdiag(dv_e/dw_e). H has the sparsity of M. The Newton
 * direction solves H d = -Phi (sparse LU).
 *
 * Globalization (PLASTICITY_NSN_GLOBALIZATION): full Newton steps (the plastic elements are
 * corrected all at once, like a primal-dual active set method; the residual is not
 * monotone), with a watchdog: after PLASTICITY_NSN_WATCHDOG_ITERATIONS steps without
 * convergence, restart from the best iterate with a merit function, the QP objective
 * f = 1/2 s^T M s + q^T s (feasible iterates, projected arc search, projected gradient
 * fallback) for convex QP problems, 1/2 ||Phi||^2 otherwise (Armijo, Levenberg-Marquardt
 * direction (H^T H + lambda I) d = -H^T Phi, lambda = c ||Phi||, when the iteration stalls
 * or H is singular, nonsmooth_newton_levenberg_marquardt.h). The stopping test is
 * plasticity_2d_compute_error
 * (the error of the NSGS solvers).
 * See docs/devel/PLASTICITY.md.
 */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "NumericsMatrix.h"
#include "PlasticityProblem.h"
#include "Plasticity_options.h"
#include "SiconosBlas.h"
#include "SolverOptions.h"
#include "naming_conventions.h"
#include "nonsmooth_newton_levenberg_marquardt.h"
#include "numerics_errors.h"
#include "numerics_verbose.h"
#include "plasticity_2d_compute_error.h"
#include "plasticity_2d_druckerprager.h"
#include "plasticity_2d_local_newton.h"
#include "plasticity_2d_solvers.h"
#include "plasticity_2d_vonmises.h"
#include "projection_soc.h"
#include "solver_registry.h"

/* data of the global natural map */
typedef struct {
  PlasticityProblem* problem;
  NumericsMatrix* M; /* sparse copy (or problem->M if already sparse) */
  int n, nc;
  double* rho; /* per element */
  int ball;    /* Von Mises: 1 ball, 0 cylinder */
  double* w;   /* work: M sigma + q */
  double* y;   /* work: sigma - rho v */
} NSNData;

/* Phi(sigma); if JP, RG are not NULL, the blocks J_P_e and J_P_e rho_e G_e (column-major,
 * 9 values per element) and returns ||Phi|| */
static double nsn_natural_map(NSNData* d, const double* sigma, double* Phi, double* JP,
                              double* JPRG) {
  PlasticityProblem* problem = d->problem;
  cblas_dcopy(d->n, problem->q, 1, d->w, 1);
  NM_gemv(1.0, d->M, (double*)sigma, 1.0, d->w);
  double norm2 = 0.0;
  for (int e = 0; e < d->nc; e++) {
    const double* s = &sigma[3 * e];
    const double* w = &d->w[3 * e];
    double* y = &d->y[3 * e];
    double rho = d->rho[e];
    double G[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1}; /* dv/dw, column-major */
    double P[3], J[9];
    if (problem->model_type == PLASTICITY_MODEL_DRUCKER_PRAGER) {
      double eta = problem->model.drucker_prager->eta[e];
      double theta = problem->model.drucker_prager->theta[e];
      double nwT = sqrt(w[1] * w[1] + w[2] * w[2]);
      y[0] = s[0] - rho * (w[0] + theta * nwT);
      y[1] = s[1] - rho * w[1];
      y[2] = s[2] - rho * w[2];
      if (nwT > 0.0) {
        G[3] = theta * w[1] / nwT; /* row 0, column 1 */
        G[6] = theta * w[2] / nwT; /* row 0, column 2 */
      }
      P[0] = y[0];
      P[1] = y[1];
      P[2] = y[2];
      soc_projection_3d(P, eta);
      if (JP) plasticity_2d_DruckerPrager_projection_jacobian(y, eta, J);
    } else {
      double sigma_y = problem->model.von_mises->sigma_y[e];
      for (int i = 0; i < 3; i++) P[i] = y[i] = s[i] - rho * w[i];
      if (d->ball)
        plasticity_2d_projectionOnVonMisesBall(P, sigma_y);
      else
        plasticity_2d_projectionOnVonMises(P, sigma_y);
      if (JP) plasticity_2d_vonMises_projection_jacobian(d->ball, y, sigma_y, J);
    }
    for (int i = 0; i < 3; i++) {
      Phi[3 * e + i] = s[i] - P[i];
      norm2 += Phi[3 * e + i] * Phi[3 * e + i];
    }
    if (JP) {
      memcpy(&JP[9 * e], J, 9 * sizeof(double));
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
          double sum = 0.0;
          for (int k = 0; k < 3; k++) sum += J[i + 3 * k] * G[k + 3 * j];
          JPRG[9 * e + i + 3 * j] = rho * sum;
        }
    }
  }
  return sqrt(norm2);
}

/* H = (I - J_P) + blockdiag(J_P R G) M (new NM_SPARSE matrix) */
static NumericsMatrix* nsn_jacobian(NSNData* d, const double* JP, const double* JPRG) {
  int n = d->n;
  NumericsMatrix* A = NM_create(NM_SPARSE, n, n); /* I - J_P */
  NumericsMatrix* B = NM_create(NM_SPARSE, n, n); /* J_P R G */
  NM_triplet_alloc(A, 9 * d->nc);
  NM_triplet_alloc(B, 9 * d->nc);
  for (int e = 0; e < d->nc; e++)
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
        double a = (i == j ? 1.0 : 0.0) - JP[9 * e + i + 3 * j];
        double b = JPRG[9 * e + i + 3 * j];
        if (a != 0.0) NM_zentry(A, 3 * e + i, 3 * e + j, a, 0.0);
        if (b != 0.0) NM_zentry(B, 3 * e + i, 3 * e + j, b, 0.0);
      }
  NumericsMatrix* BM = NM_multiply(B, d->M);
  NumericsMatrix* H = NM_add(1.0, A, 1.0, BM);
  NM_free(A);
  NM_free(B);
  NM_free(BM);
  return H;
}

/* projection of each element of x on its yield set */
static void nsn_project(NSNData* d, double* x) {
  PlasticityProblem* problem = d->problem;
  for (int e = 0; e < d->nc; e++) {
    double* s = &x[3 * e];
    if (problem->model_type == PLASTICITY_MODEL_DRUCKER_PRAGER)
      soc_projection_3d(s, problem->model.drucker_prager->eta[e]);
    else if (d->ball)
      plasticity_2d_projectionOnVonMisesBall(s, problem->model.von_mises->sigma_y[e]);
    else
      plasticity_2d_projectionOnVonMises(s, problem->model.von_mises->sigma_y[e]);
  }
}

/* f = 1/2 s^T M s + q^T s = 1/2 s^T (w + q), w = M s + q (left in d->w) */
static double nsn_qp_objective(NSNData* d, const double* s) {
  cblas_dcopy(d->n, d->problem->q, 1, d->w, 1);
  NM_gemv(1.0, d->M, (double*)s, 1.0, d->w);
  double f = 0.0;
  for (int i = 0; i < d->n; i++) f += 0.5 * s[i] * (d->w[i] + d->problem->q[i]);
  return f;
}

/* Armijo search along the projected arc s(t) = P_K(s + t dir), t = 1, 1/2, ...:
   f(s(t)) <= f0 + armijo w^T (s(t) - s) + 16 eps |f0|, w the gradient at s (the last term
   accepts the steps whose decrease is below the round-off of f, close to the solution).
   Returns 1 and s(t) in s_t if accepted. */
static int nsn_arc_search(NSNData* d, const double* s, const double* dir, const double* w,
                          double f0, double armijo, int max_iter, double* s_t) {
  double t = 1.0;
  for (int k = 0; k < max_iter; k++, t *= 0.5) {
    for (int i = 0; i < d->n; i++) s_t[i] = s[i] + t * dir[i];
    nsn_project(d, s_t);
    double slope = 0.0;
    for (int i = 0; i < d->n; i++) slope += w[i] * (s_t[i] - s[i]);
    if (slope >= 0.0) continue; /* not a descent point of the arc */
    if (nsn_qp_objective(d, s_t) <= f0 + armijo * slope + 16.0 * DBL_EPSILON * fabs(f0)) {
      numerics_printf_verbose(2, "---- PLASTICITY_2D_NSN -   arc step t = %.3e", t);
      return 1;
    }
  }
  return 0;
}

/* the problem is a convex QP: associated flow (Von Mises, Drucker-Prager with theta = 0)
   and M symmetric */
static int nsn_is_qp(NSNData* d, SolverOptions* options) {
  (void)options;
  PlasticityProblem* problem = d->problem;
  if (problem->model_type == PLASTICITY_MODEL_DRUCKER_PRAGER)
    for (int e = 0; e < d->nc; e++)
      if (problem->model.drucker_prager->theta[e] != 0.0) return 0;
  NumericsMatrix* Mt = NM_transpose(d->M);
  NumericsMatrix* D = NM_add(1.0, d->M, -1.0, Mt);
  double asym = NM_norm_inf(D), norm = NM_norm_inf(d->M);
  NM_free(Mt);
  NM_free(D);
  return asym <= 1e-12 * norm;
}

/* rho_e of the natural map of an element (column-major 3x3 diagonal block M_ee), before the
   factor dparam[PLASTICITY_NSN_RHO_FACTOR] of the full Newton steps: INVERSE_SPECTRAL_NORM:
   ||M_ee^{-1}||_2 = 1/lambda_min (default; 1/||M_ee||_2 if M_ee is singular), SPECTRAL_NORM:
   1/||M_ee||_2. rho decides which elements a full Newton step treats as plastic (y = sigma -
   rho v); on the FEM problems (strip with a hole) the number of full steps decreases when rho
   increases up to ~30 ||M_ee^{-1}|| and the iteration diverges beyond ~100 (footing), see
   docs/devel/PLASTICITY.md. */
static double nsn_element_rho(const double M[9], int strategy) {
  double n2 = plasticity_2d_local_spectral_norm(M);
  double rho_max = (n2 > 0.0) ? 1.0 / n2 : 1.0;
  if (strategy == PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_SPECTRAL_NORM) return rho_max;
  double det = M[0] * (M[4] * M[8] - M[7] * M[5]) - M[3] * (M[1] * M[8] - M[7] * M[2]) +
               M[6] * (M[1] * M[5] - M[4] * M[2]);
  if (!(fabs(det) > DBL_EPSILON * n2 * n2 * n2)) return rho_max;
  double inv[9];
  inv[0] = (M[4] * M[8] - M[7] * M[5]) / det;
  inv[3] = -(M[3] * M[8] - M[6] * M[5]) / det;
  inv[6] = (M[3] * M[7] - M[6] * M[4]) / det;
  inv[1] = -(M[1] * M[8] - M[7] * M[2]) / det;
  inv[4] = (M[0] * M[8] - M[6] * M[2]) / det;
  inv[7] = -(M[0] * M[7] - M[6] * M[1]) / det;
  inv[2] = (M[1] * M[5] - M[4] * M[2]) / det;
  inv[5] = -(M[0] * M[5] - M[3] * M[2]) / det;
  inv[8] = (M[0] * M[4] - M[3] * M[1]) / det;
  return plasticity_2d_local_spectral_norm(inv);
}

void plasticity_2d_nsn(PlasticityProblem* problem, double* stress, double* strainrate,
                       int* info, SolverOptions* options) {
  *info = 1;
  int nc = problem->numberOfCones;
  int n = 3 * nc;
  if (problem->dimension != 3) {
    *info = numerics_error("plasticity_2d_nsn", "dimension must be 3");
    return;
  }
  if (problem->stress_coordinates != PLASTICITY_STRESS_COORDINATES_CONE ||
      (problem->model_type == PLASTICITY_MODEL_VON_MISES &&
       problem->model.von_mises->hypothesis == PLASTICITY_VONMISES_PLANE_STRESS) ||
      (problem->model_type == PLASTICITY_MODEL_DRUCKER_PRAGER &&
       problem->model.drucker_prager->cohesion)) {
    *info = numerics_error("plasticity_2d_nsn",
                           "Cartesian, plane stress or cohesive problem: call "
                           "plasticity_2d_driver (change of variables)");
    return;
  }
  if (problem->model_type != PLASTICITY_MODEL_VON_MISES &&
      problem->model_type != PLASTICITY_MODEL_DRUCKER_PRAGER) {
    *info = numerics_error("plasticity_2d_nsn", "unknown plasticity model");
    return;
  }

  NSNData d;
  d.problem = problem;
  d.n = n;
  d.nc = nc;
  d.ball = (problem->model_type == PLASTICITY_MODEL_VON_MISES &&
            problem->model.von_mises->hypothesis == PLASTICITY_VONMISES_BALL);
  d.M =
      (problem->M->storageType == NM_SPARSE) ? problem->M : NM_copy_to_sparse(problem->M, 0.0);
  d.rho = (double*)malloc((size_t)nc * sizeof(double));
  d.w = (double*)malloc((size_t)n * sizeof(double));
  d.y = (double*)malloc((size_t)n * sizeof(double));
  double block[9];
  double* pblock = block;
  int rho_strategy = options->iparam[PLASTICITY_NSN_RHO_STRATEGY];
  double rho_factor = options->dparam[PLASTICITY_NSN_RHO_FACTOR] > 0.0
                          ? options->dparam[PLASTICITY_NSN_RHO_FACTOR]
                          : 1.0;
  for (int e = 0; e < nc; e++) {
    NM_copy_diag_block3(problem->M, e, &pblock);
    d.rho[e] = nsn_element_rho(block, rho_strategy);
  }

  double* Phi = (double*)malloc((size_t)n * sizeof(double));
  double* Phi_t = (double*)malloc((size_t)n * sizeof(double));
  double* dir = (double*)malloc((size_t)n * sizeof(double));
  double* sigma_t = (double*)malloc((size_t)n * sizeof(double));
  double* JP = (double*)malloc((size_t)(9 * nc) * sizeof(double));
  double* JPRG = (double*)malloc((size_t)(9 * nc) * sizeof(double));

  int max_iter = SOLVER_MAX_ITER(options);
  double tolerance = SOLVER_TOL(options);
  int ls_max_iter = options->iparam[PLASTICITY_NSN_LINESEARCH_MAX_ITER] > 0
                        ? options->iparam[PLASTICITY_NSN_LINESEARCH_MAX_ITER]
                        : 30;
  int use_lm = options->iparam[PLASTICITY_NSN_LM_FALLBACK];
  double lm_coef = options->dparam[PLASTICITY_NSN_LM_COEF] > 0.0
                       ? options->dparam[PLASTICITY_NSN_LM_COEF]
                       : 1.0;
  double norm_q = cblas_dnrm2(n, problem->q, 1);
  const double armijo = 1e-4;

  /* Globalization (iparam[PLASTICITY_NSN_GLOBALIZATION]):
     - full Newton steps (WATCHDOG, NONE): the active set is corrected for all the elements at
       once (primal-dual active set behaviour), the residual is not monotone;
     - merit function (MERIT, and WATCHDOG after PLASTICITY_NSN_WATCHDOG_ITERATIONS full
       steps without convergence, from the best iterate): the QP objective
       f = 1/2 s^T M s + q^T s (feasible iterates, projected arc search, projected gradient
       fallback) when the problem is a convex QP (associated flow: Von Mises, or
       Drucker-Prager without dilatancy, and M symmetric); 1/2 ||Phi||^2 otherwise (Armijo,
       Levenberg-Marquardt fallback). */
  int globalization = options->iparam[PLASTICITY_NSN_GLOBALIZATION];
  int full_step = (globalization != PLASTICITY_NSN_GLOBALIZATION_MERIT);
  int qp_problem = nsn_is_qp(&d, options);
  int qp = 0; /* QP merit in use */
  double* best_stress = (double*)malloc((size_t)n * sizeof(double));
  double best_error = INFINITY;
  memcpy(best_stress, stress, (size_t)n * sizeof(double));

  double error = INFINITY, previous_norm = INFINITY;
  int stalled = 0, iter = 0, n_lm = 0, n_gradient = 0, n_full = 0;
  /* the factor c of rho applies to the full Newton steps only (the merit phase converges
     faster with c = 1) */
  if (full_step)
    for (int e = 0; e < nc; e++) d.rho[e] *= rho_factor;
  if (!full_step && qp_problem) {
    qp = 1;
    nsn_project(&d, stress);
  }
  double normPhi = nsn_natural_map(&d, stress, Phi, JP, JPRG);
  numerics_printf_verbose(1, "---- PLASTICITY_2D_NSN - start, %d elements, %s", nc,
                          full_step ? "full Newton steps"
                          : qp      ? "merit: QP objective"
                                    : "merit: natural map");
  while (iter < max_iter) {
    plasticity_2d_compute_error(problem, stress, strainrate, tolerance, options, norm_q,
                                &error);
    numerics_printf_verbose(
        1, "---- PLASTICITY_2D_NSN - iteration %d, |Phi| = %.3e, error = %.3e", iter, normPhi,
        error);
    if (verbose >= 2) {
      int n_active = 0; /* elements whose projection is not the identity */
      for (int e = 0; e < nc; e++)
        if (fabs(JP[9 * e] - 1.0) + fabs(JP[9 * e + 4] - 1.0) + fabs(JP[9 * e + 8] - 1.0) >
            1e-14)
          n_active++;
      numerics_printf_verbose(2, "---- PLASTICITY_2D_NSN -   active elements %d", n_active);
    }
    if (error <= tolerance) {
      *info = 0;
      break;
    }
    if (full_step && isfinite(error) && error < best_error) {
      best_error = error;
      memcpy(best_stress, stress, (size_t)n * sizeof(double));
    }
    /* watchdog: switch to the merit function from the best iterate */
    if (full_step && globalization == PLASTICITY_NSN_GLOBALIZATION_WATCHDOG &&
        (n_full >= PLASTICITY_NSN_WATCHDOG_ITERATIONS || !isfinite(error))) {
      full_step = 0;
      for (int e = 0; e < nc; e++) d.rho[e] /= rho_factor;
      memcpy(stress, best_stress, (size_t)n * sizeof(double));
      if (qp_problem) {
        qp = 1;
        nsn_project(&d, stress);
      }
      normPhi = nsn_natural_map(&d, stress, Phi, JP, JPRG);
      previous_norm = INFINITY;
      numerics_printf_verbose(1, "---- PLASTICITY_2D_NSN - watchdog: switch to the %s merit",
                              qp ? "QP objective" : "natural map");
      continue;
    }
    iter++;

    /* Newton direction H d = -Phi (Levenberg-Marquardt when stalled / singular) */
    NumericsMatrix* H = nsn_jacobian(&d, JP, JPRG);
    int lm_now = !qp && !full_step && use_lm &&
                 nonsmooth_newton_lm_is_stalled(normPhi, previous_norm, &stalled);
    int lu_info = 1;
    for (int i = 0; i < n; i++) dir[i] = -Phi[i];
    if (!lm_now) {
      lu_info = NM_LU_solve(NM_preserve(H), dir, 1);
      for (int i = 0; i < n && !lu_info; i++)
        if (!isfinite(dir[i])) lu_info = 1;
    }
    if (lm_now || (lu_info && use_lm)) {
      n_lm++;
      for (int i = 0; i < n; i++) dir[i] = 0.0;
      lu_info = nonsmooth_newton_lm_direction(NM_unpreserve(H), Phi, lm_coef * normPhi, dir);
    }
    NM_free(H);

    int accepted = 0;
    if (full_step) {
      n_full++;
      if (lu_info) {
        /* singular system: let the watchdog switch to the merit function */
        if (globalization == PLASTICITY_NSN_GLOBALIZATION_NONE) break;
        n_full = PLASTICITY_NSN_WATCHDOG_ITERATIONS;
        continue;
      }
      for (int i = 0; i < n; i++) sigma_t[i] = stress[i] + dir[i];
    } else if (qp) {
      /* projected arc s(t) = P_K(s + t d), Armijo on f: f(s(t)) <= f(s) + a w^T (s(t) - s) */
      double f0 = nsn_qp_objective(&d, stress);       /* also w = M s + q in d.w */
      memcpy(Phi_t, d.w, (size_t)n * sizeof(double)); /* gradient w at s */
      if (!lu_info)
        accepted = nsn_arc_search(&d, stress, dir, Phi_t, f0, armijo, ls_max_iter, sigma_t);
      if (!accepted) {
        /* projected gradient: d = -R w (descent for any t small enough) */
        n_gradient++;
        for (int e = 0; e < nc; e++)
          for (int i = 0; i < 3; i++) dir[3 * e + i] = -d.rho[e] * Phi_t[3 * e + i];
        accepted =
            nsn_arc_search(&d, stress, dir, Phi_t, f0, armijo, 2 * ls_max_iter, sigma_t);
      }
      if (!accepted) {
        numerics_printf_verbose(1, "---- PLASTICITY_2D_NSN - no descent step");
        break;
      }
    } else {
      if (lu_info) {
        numerics_printf_verbose(1, "---- PLASTICITY_2D_NSN - singular Newton system");
        break;
      }
      /* Armijo line search on psi = 1/2 |Phi|^2 */
      double psi = 0.5 * normPhi * normPhi, t = 1.0;
      for (int ls = 0; ls < ls_max_iter; ls++) {
        for (int i = 0; i < n; i++) sigma_t[i] = stress[i] + t * dir[i];
        double normPhi_t = nsn_natural_map(&d, sigma_t, Phi_t, NULL, NULL);
        if (0.5 * normPhi_t * normPhi_t <= (1.0 - 2.0 * armijo * t) * psi) {
          accepted = 1;
          break;
        }
        t *= 0.5;
      }
      if (!accepted)
        /* take the smallest step anyway (the stall test switches to Levenberg-Marquardt) */
        numerics_printf_verbose(1, "---- PLASTICITY_2D_NSN - line search failed (t = %.1e)",
                                t);
    }
    memcpy(stress, sigma_t, (size_t)n * sizeof(double));
    previous_norm = normPhi;
    normPhi = nsn_natural_map(&d, stress, Phi, JP, JPRG);
  }
  /* full steps only: return the best iterate */
  if (*info && globalization == PLASTICITY_NSN_GLOBALIZATION_NONE && best_error < error)
    memcpy(stress, best_stress, (size_t)n * sizeof(double));
  free(best_stress);
  if (*info) {
    plasticity_2d_compute_error(problem, stress, strainrate, tolerance, options, norm_q,
                                &error);
    *info = (error <= tolerance) ? 0 : 1;
  }
  numerics_printf_verbose(
      1,
      "---- PLASTICITY_2D_NSN - end: info %d, %d iterations (%d full "
      "steps, %d Levenberg-Marquardt, %d projected gradient), error = %.3e",
      *info, iter, n_full, n_lm, n_gradient, error);
  SET_SOLVER_ITER_DONE(options, iter);
  SET_SOLVER_RESIDUAL(options, error);

  if (d.M != problem->M) NM_free(d.M);
  free(d.rho);
  free(d.w);
  free(d.y);
  free(Phi);
  free(Phi_t);
  free(dir);
  free(sigma_t);
  free(JP);
  free(JPRG);
}

void plasticity_2d_nsn_set_default(SolverOptions* options) {
  options->iparam[PLASTICITY_NSN_LINESEARCH_MAX_ITER] = 30;
  options->iparam[PLASTICITY_NSN_LM_FALLBACK] = 1;
  options->iparam[PLASTICITY_NSN_GLOBALIZATION] = PLASTICITY_NSN_GLOBALIZATION_WATCHDOG;
  options->iparam[PLASTICITY_NSN_RHO_STRATEGY] =
      PLASTICITY_NSN_FORMULATION_RHO_STRATEGY_INVERSE_SPECTRAL_NORM;
  options->dparam[PLASTICITY_NSN_RHO_FACTOR] = 10.0;
  options->dparam[PLASTICITY_NSN_LM_COEF] = 1.0;
}

static int plasticity_2d_nsn_init_wrap(void* problem, SolverOptions* options) {
  (void)problem;
  (void)options;
  return NUMERICS_OK;
}

static int plasticity_2d_nsn_solve_wrap(void* problem, double* reaction, double* velocity,
                                        SolverOptions* options) {
  int info = 1;
  plasticity_2d_nsn((PlasticityProblem*)problem, reaction, velocity, &info, options);
  return info;
}

REGISTER_SOLVER(PLASTICITY_2D_NSN, "PLASTICITY_2D_NSN",
                "Global semi-smooth Newton on the natural map (plasticity)",
                plasticity_2d_nsn_init_wrap, plasticity_2d_nsn_solve_wrap, NULL, NULL,
                plasticity_2d_nsn_set_default, 100, /* default_max_iter */
                1e-10,                              /* default_tol */
                0 /* is_local_solver */);
