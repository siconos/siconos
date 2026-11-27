// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "convex_set.h"

static double max_abs_diff(const double* a, const double* b, int n) {
  double d = 0.0;
  for (int i = 0; i < n; ++i) {
    double di = fabs(a[i] - b[i]);
    if (di > d) d = di;
  }
  return d;
}

int main(void) {
  int info = 0;
  const int nc = 6;
  const int n = 3 * nc;
  const double rho = 1.5;

  /* reaction and velocity stacked by contact */
  double r[] = {/* contact 0: w inside cone */
                2.0, 0.01, 0.0,
                /* contact 1: w in dual cone */
                0.0, 0.5, 0.0,
                /* contact 2: w on cone boundary */
                1.2, 1.0, 0.0,
                /* contact 3: w inside, nonzero tangential u */
                3.0, 0.2, -0.1,
                /* contact 4: w on boundary, nonzero normal u */
                1.5, 1.0, 0.0,
                /* contact 5: w in dual, nonzero u */
                1.0, 0.5, 0.5};
  double u[] = {1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                1.0, 0.1, 0.0, 0.5, 0.0, 0.0, 2.0, 0.0, 0.0};
  double mu[] = {0.1, 0.3, 0.5, 0.7, 0.5, 0.2};

  double F_loop[3 * 6];
  double F_vec[3 * 6];
  double v[3 * 6];
  for (int i = 0; i < n; ++i) v[i] = rho * u[i];
  ConvexSet K = convex_set(CONVEX_SET_SOC, 3, nc, mu, NULL);

  /* reference: cone-wise natural map r - P_K(r - rho u) */
  double F2_loop = 0.0;
  for (int i = 0; i < nc; ++i) {
    F2_loop += convex_set_natural_map_block(&K, i, &r[3 * i], &v[3 * i], &F_loop[3 * i]);
  }

  /* natural map on the whole stack (vectorized projection) */
  memset(F_vec, 0, sizeof(F_vec));
  double F2_vec = convex_set_natural_map(&K, r, v, F_vec);
  double err = max_abs_diff(F_loop, F_vec, n);
  printf("convex_set_natural_map max diff = %e\n", err);
  if (err > 1e-14 || fabs(F2_loop - F2_vec) > 1e-14) {
    fprintf(stderr, "Stacked natural map differs from the cone-wise version!\n");
    info += 1;
  }

  /* contact 1: r - rho u in the polar cone, F = r */
  if (max_abs_diff(&F_loop[3], &r[3], 3) > 0.0) {
    fprintf(stderr, "natural map in the polar cone should be r!\n");
    info += 1;
  }

  return info;
}
