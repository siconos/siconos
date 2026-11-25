/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
/* Tools of the 2D rolling friction tests: random problems W = H M^-1 H^T
   of a chain of disks (contact c touches the disks c and c - 1, rows
   normal, tangent, rolling), and the check of the law independent of
   rolling_friction_2d_compute_error. */
#ifndef RFC2D_TEST_TOOLS_H
#define RFC2D_TEST_TOOLS_H

#include <math.h>
#include <stdlib.h>

static inline double rfc2d_urand(double a, double b) {
  return a + (b - a) * rand() / (double)RAND_MAX;
}

/* W (column major, n x n, n = 3 nc), q, mu, mu_r of a random chain;
   with_mur = 0: mu_r = 0 */
static inline void rfc2d_random_chain(int nc, int with_mur, double* W, double* q, double* mu,
                                      double* mur) {
  int n = 3 * nc, nb = 3 * nc;
  double* H = (double*)calloc((size_t)n * nb, sizeof(double));
  double* Minv = (double*)calloc(nb, sizeof(double));
  for (int b = 0; b < nc; b++) {
    double m = rfc2d_urand(0.1, 10.), R = rfc2d_urand(0.1, 3.);
    Minv[3 * b] = Minv[3 * b + 1] = 1. / m;
    Minv[3 * b + 2] = 1. / (0.5 * m * R * R);
  }
  for (int c = 0; c < nc; c++) {
    double a = rfc2d_urand(0, 2 * M_PI), nx = cos(a), ny = sin(a), tx = -ny, ty = nx;
    int r0 = 3 * c * nb, r1 = (3 * c + 1) * nb, r2 = (3 * c + 2) * nb;
    double R1 = rfc2d_urand(0.1, 3.);
    H[r0 + 3 * c] = nx;
    H[r0 + 3 * c + 1] = ny;
    H[r1 + 3 * c] = tx;
    H[r1 + 3 * c + 1] = ty;
    H[r1 + 3 * c + 2] = -R1;
    H[r2 + 3 * c + 2] = 1.;
    if (c > 0) {
      int b2 = 3 * (c - 1);
      double R2 = rfc2d_urand(0.1, 3.);
      H[r0 + b2] = -nx;
      H[r0 + b2 + 1] = -ny;
      H[r1 + b2] = -tx;
      H[r1 + b2 + 1] = -ty;
      H[r1 + b2 + 2] = -R2;
      H[r2 + b2 + 2] = -1.;
    }
  }
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      double s = 0.;
      for (int l = 0; l < nb; l++) s += H[i * nb + l] * Minv[l] * H[j * nb + l];
      W[i + j * n] = s;
    }
  for (int i = 0; i < n; i++) q[i] = rfc2d_urand(-1., 1.);
  for (int c = 0; c < nc; c++) {
    mu[c] = rfc2d_urand(0.05, 1.);
    mur[c] = with_mur ? rfc2d_urand(0.001, 0.5) : 0.;
  }
  free(H);
  free(Minv);
}

/* the largest violation of the law by r: r in K, u = W r + q with
   u_n + mu |u_t| + mu_r |u_r| >= 0 (u tilde in the dual cone), and the
   complementarity r . u tilde = 0; relative to 1 + |r_n| */
static inline double rfc2d_law_violation(int nc, const double* W, const double* q,
                                         const double* mu, const double* mur,
                                         const double* r) {
  int n = 3 * nc;
  double violation = 0.;
  for (int c = 0; c < nc; c++) {
    const double* rc = r + 3 * c;
    double uc[3];
    for (int i = 0; i < 3; i++) {
      double s = q[3 * c + i];
      for (int j = 0; j < n; j++) s += W[(3 * c + i) + j * n] * r[j];
      uc[i] = s;
    }
    double scale = 1. + fabs(rc[0]);
    double ut0 = uc[0] + mu[c] * fabs(uc[1]) + mur[c] * fabs(uc[2]);
    violation = fmax(violation, (fabs(rc[1]) - mu[c] * rc[0]) / scale);
    violation = fmax(violation, (fabs(rc[2]) - mur[c] * rc[0]) / scale);
    violation = fmax(violation, -rc[0] / scale);
    violation = fmax(violation, -uc[0]);
    violation = fmax(violation, fabs(rc[0] * ut0 + rc[1] * uc[1] + rc[2] * uc[2]) / scale);
  }
  return violation;
}

#endif
