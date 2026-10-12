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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fc3d_compute_error.h"

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

  /* Local velocities covering various tangential/normal magnitudes. */
  double u[] = {
      1.0, 0.0, 0.0,    /* pure normal */
      -0.5, 0.3, 0.4,   /* negative normal with tangential part */
      0.2, 1.0, 0.0,    /* boundary-like */
      2.0, 0.2, -0.3,   /* inside-like */
      -1.0, 0.5, 0.5,   /* dual-like */
      0.0, 0.1, 0.1};
  double mu[] = {0.1, 0.3, 0.5, 0.7, 0.5, 0.2};

  double u_mod[3 * 6];
  double expected[3 * 6];
  double work[6];

  /* Expected modified local velocity computed cone by cone. */
  for (int i = 0; i < nc; ++i) {
    const int pos = 3 * i;
    const double normUT = sqrt(u[pos + 1] * u[pos + 1] + u[pos + 2] * u[pos + 2]);
    expected[pos] = u[pos] + mu[i] * normUT;
    expected[pos + 1] = u[pos + 1];
    expected[pos + 2] = u[pos + 2];
  }

  /* Vectorized version with caller-provided workspace. */
  fc3d_compute_modified_local_velocity(nc, u, mu, u_mod, work);
  double err = max_abs_diff(expected, u_mod, n);
  printf("fc3d_compute_modified_local_velocity (with work) max diff = %e\n", err);
  if (err > 1e-14) {
    fprintf(stderr, "Vectorized modified velocity differs from cone-wise reference!\n");
    info += 1;
  }

  /* Vectorized version with internal allocation. */
  memset(u_mod, 0, sizeof(u_mod));
  fc3d_compute_modified_local_velocity(nc, u, mu, u_mod, NULL);
  err = max_abs_diff(expected, u_mod, n);
  printf("fc3d_compute_modified_local_velocity (work=NULL) max diff = %e\n", err);
  if (err > 1e-14) {
    fprintf(stderr, "Vectorized modified velocity (NULL work) differs from cone-wise reference!\n");
    info += 1;
  }

  /* Consistency check: r - u_modified must equal the worktmp used in
     fc3d_unitary_compute_and_add_error before projection. */
  double r[] = {
      0.5, 0.1, 0.0,
      0.2, -0.1, 0.2,
      1.0, 0.0, 0.0,
      -0.3, 0.2, 0.1,
      0.0, 0.5, -0.5,
      0.7, 0.0, 0.0};
  double diff_expected[3 * 6];
  for (int i = 0; i < nc; ++i) {
    const int pos = 3 * i;
    const double normUT = sqrt(u[pos + 1] * u[pos + 1] + u[pos + 2] * u[pos + 2]);
    diff_expected[pos] = r[pos] - (u[pos] + mu[i] * normUT);
    diff_expected[pos + 1] = r[pos + 1] - u[pos + 1];
    diff_expected[pos + 2] = r[pos + 2] - u[pos + 2];
  }

  double diff[3 * 6];
  fc3d_compute_modified_local_velocity(nc, u, mu, u_mod, work);
  for (int i = 0; i < n; ++i) {
    diff[i] = r[i] - u_mod[i];
  }
  err = max_abs_diff(diff_expected, diff, n);
  printf("r - u_modified consistency max diff = %e\n", err);
  if (err > 1e-14) {
    fprintf(stderr, "r - u_modified differs from expected complementarity argument!\n");
    info += 1;
  }

  return info;
}
