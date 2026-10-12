/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "projectionOnCone.h"

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
  double r[] = {
      /* contact 0: w inside cone */
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
  double u[] = {
      1.0, 0.0, 0.0,
      1.0, 0.0, 0.0,
      1.0, 0.0, 0.0,
      1.0, 0.1, 0.0,
      0.5, 0.0, 0.0,
      2.0, 0.0, 0.0};
  double mu[] = {0.1, 0.3, 0.5, 0.7, 0.5, 0.2};

  double F_loop[3 * 6];
  double F_vec[3 * 6];
  double work[6];

  /* reference: cone-wise natural map */
  for (int i = 0; i < nc; ++i) {
    naturalMapOnCone(&r[3 * i], &u[3 * i], rho, mu[i], &F_loop[3 * i]);
  }

  /* vectorized version with caller-provided workspace */
  naturalMapOnConeVector(nc, r, u, rho, mu, F_vec, work);
  double err = max_abs_diff(F_loop, F_vec, n);
  printf("naturalMapOnConeVector (with work) max diff = %e\n", err);
  if (err > 1e-14) {
    fprintf(stderr, "Vectorized natural map differs from cone-wise version!\n");
    info += 1;
  }

  /* vectorized version with internal allocation */
  memset(F_vec, 0, sizeof(F_vec));
  naturalMapOnConeVector(nc, r, u, rho, mu, F_vec, NULL);
  err = max_abs_diff(F_loop, F_vec, n);
  printf("naturalMapOnConeVector (work=NULL) max diff = %e\n", err);
  if (err > 1e-14) {
    fprintf(stderr, "Vectorized natural map (NULL work) differs from cone-wise version!\n");
    info += 1;
  }

  return info;
}
