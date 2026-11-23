// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "plasticity_2d_vonmises_projection.h"

#include <float.h>
#include <math.h>
#include <stdio.h>

/* #define DEBUG_NOCOLOR */
/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "siconos_debug.h"

double plasticity_2d_vonMises_radius(double sigma_y) { return sqrt(2.0 / 3.0) * sigma_y; }

double plasticity_2d_vonMises_equivalent_stress(const double stress[3]) {
  return sqrt(1.5) * hypot(stress[1], stress[2]);
}

int plasticity_2d_vonMises_check_yield(const double stress[3], double sigma_y) {
  return (hypot(stress[1], stress[2]) <= plasticity_2d_vonMises_radius(sigma_y)) ? 1 : 0;
}

int plasticity_2d_projectionOnVonMises(double stress[3], double sigma_y) {
  double R = plasticity_2d_vonMises_radius(sigma_y);
  double nT = hypot(stress[1], stress[2]);
  DEBUG_PRINTF("Von Mises projection: |sigma_T| = %e, R = %e\n", nT, R);
  if (nT <= R) return 1;
  /* radial return of the deviatoric part, the hydrostatic part is unchanged */
  double scale = R / nT;
  stress[1] *= scale;
  stress[2] *= scale;
  return 0;
}

void plasticity_2d_vonmises_unitary_compute_and_add_error(const double r[3], const double u[3],
                                                          double sigma_y, double* error,
                                                          double worktmp[3]) {
  for (int i = 0; i < 3; i++) worktmp[i] = r[i] - u[i];
  plasticity_2d_projectionOnVonMises(worktmp, sigma_y);
  for (int i = 0; i < 3; i++) *error += (r[i] - worktmp[i]) * (r[i] - worktmp[i]);
}

/* x = -(S + mu I)^{-1} b for the 2x2 matrix S. Returns 0 if (S + mu I) is invertible. */
static int vonmises_solve2(const double S[4], double mu, const double b[2], double x[2]) {
  double a = S[0] + mu, d = S[3] + mu; /* S = [S[0] S[2]; S[1] S[3]] (column-major) */
  double det = a * d - S[2] * S[1];
  if (!(fabs(det) > DBL_MIN)) return 1;
  x[0] = -(d * b[0] - S[2] * b[1]) / det;
  x[1] = -(-S[1] * b[0] + a * b[1]) / det;
  return 0;
}

int plasticity_2d_vonMises_local_solve(const double M[9], const double q[3], double sigma_y,
                                       double stress[3]) {
  double m = M[0];
  if (!(m > 0.0)) return 1;
  /* M = [m, M_NT; M_TN, M_TT] (column-major) */
  const double mNT[2] = {M[3], M[6]};
  const double mTN[2] = {M[1], M[2]};
  /* Schur complement S = M_TT - M_TN M_NT / m and q~ = q_T - M_TN q_N / m */
  const double S[4] = {M[4] - mTN[0] * mNT[0] / m, M[5] - mTN[1] * mNT[0] / m,
                       M[7] - mTN[0] * mNT[1] / m, M[8] - mTN[1] * mNT[1] / m};
  const double qt[2] = {q[1] - mTN[0] * q[0] / m, q[2] - mTN[1] * q[0] / m};
  double R = plasticity_2d_vonMises_radius(sigma_y);
  double x[2];
  int status = 0;

  /* elastic trial: S x = -q~ */
  if (vonmises_solve2(S, 0.0, qt, x)) return 1;
  double nx = hypot(x[0], x[1]);
  DEBUG_PRINTF("Von Mises local solve: elastic trial |sigma_T| = %e, R = %e\n", nx, R);

  if (nx > R) {
    if (R <= 0.0) {
      x[0] = x[1] = 0.0;
    } else {
      /* plastic: find mu > 0 with ||x(mu)|| = R, x(mu) = -(S + mu I)^{-1} q~.
       * g(mu) = 1/||x(mu)|| - 1/R is increasing (and concave for a symmetric S): Newton from
       * the left of the root, safeguarded by bisection on the bracket [lo, hi]. */
      double lo = 0.0, hi = hypot(qt[0], qt[1]) / R;
      double xt[2];
      int k = 0;
      /* ||x(mu)|| <= ||q~|| / mu when the symmetric part of S is positive semi-definite; in
       * general, enlarge hi until ||x(hi)|| <= R */
      while (k++ < 100 && (vonmises_solve2(S, hi, qt, xt) || hypot(xt[0], xt[1]) > R)) {
        lo = hi;
        hi *= 2.0;
      }
      double mu = lo;
      double f = nx - R;
      for (int it = 0; it < 100; it++) {
        if (vonmises_solve2(S, mu, qt, x)) {
          status = 1;
          break;
        }
        nx = hypot(x[0], x[1]);
        f = nx - R;
        DEBUG_PRINTF("  it %d, mu = %.17e, ||x|| - R = %e\n", it, mu, f);
        if (fabs(f) <= 4.0 * DBL_EPSILON * R) break;
        if (f > 0.0)
          lo = mu;
        else
          hi = mu;
        if (hi - lo <= 2.0 * DBL_EPSILON * hi) break;
        /* Newton step on g: dx/dmu = -(S + mu I)^{-1} x */
        double dx[2];
        vonmises_solve2(S, mu, x, dx);
        double dnx = (x[0] * dx[0] + x[1] * dx[1]) / nx; /* d||x||/dmu */
        double g = 1.0 / nx - 1.0 / R, dg = -dnx / (nx * nx);
        double mu_new = (dg != 0.0) ? mu - g / dg : 0.5 * (lo + hi);
        if (!(mu_new > lo && mu_new < hi)) mu_new = 0.5 * (lo + hi);
        mu = mu_new;
      }
      if (status == 0 && fabs(f) > 1e3 * DBL_EPSILON * R) status = 2;
      /* on the yield surface */
      nx = hypot(x[0], x[1]);
      if (nx > 0.0) {
        x[0] *= R / nx;
        x[1] *= R / nx;
      }
    }
  }
  stress[1] = x[0];
  stress[2] = x[1];
  /* w_N = 0 */
  stress[0] = -(q[0] + mNT[0] * x[0] + mNT[1] * x[1]) / m;
  return status;
}
