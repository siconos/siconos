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
#include <sys/time.h>

#include "projectionOnCone.h"

static double now(void) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return (double)tv.tv_sec + (double)tv.tv_usec * 1e-6;
}

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
  const int nc = 200000;
  const int n = 3 * nc;
  const int nrep = 20;
  const double rho = 1.5;

  double* r = (double*)malloc(n * sizeof(double));
  double* u = (double*)malloc(n * sizeof(double));
  double* mu = (double*)malloc(nc * sizeof(double));
  double* work = (double*)malloc(nc * sizeof(double));

  srand48(123456);
  for (int i = 0; i < nc; ++i) {
    mu[i] = 0.1 + 0.9 * drand48();
    r[3 * i] = 2.0 * drand48() - 1.0;
    r[3 * i + 1] = 2.0 * drand48() - 1.0;
    r[3 * i + 2] = 2.0 * drand48() - 1.0;
    u[3 * i] = 2.0 * drand48() - 1.0;
    u[3 * i + 1] = 2.0 * drand48() - 1.0;
    u[3 * i + 2] = 2.0 * drand48() - 1.0;
  }

  double* r_loop = (double*)malloc(n * sizeof(double));
  double* r_vec = (double*)malloc(n * sizeof(double));
  double* F_loop = (double*)malloc(n * sizeof(double));
  double* F_vec = (double*)malloc(n * sizeof(double));

  /* ---------- projectionOnConeVector ---------- */
  double t0 = now();
  for (int rep = 0; rep < nrep; ++rep) {
    memcpy(r_loop, r, n * sizeof(double));
    for (int i = 0; i < nc; ++i) {
      projectionOnCone(&r_loop[3 * i], mu[i]);
    }
  }
  double t_loop_proj = (now() - t0) / nrep;

  t0 = now();
  for (int rep = 0; rep < nrep; ++rep) {
    memcpy(r_vec, r, n * sizeof(double));
    projectionOnConeVector(nc, r_vec, mu, work);
  }
  double t_vec_proj = (now() - t0) / nrep;

  double err_proj = max_abs_diff(r_loop, r_vec, n);
  printf(
      "projectionOnConeVector: nc=%d  cone-wise=%.6f s  vectorized=%.6f s  speedup=%.2fx  "
      "err=%e\n",
      nc, t_loop_proj, t_vec_proj, t_loop_proj / fmax(t_vec_proj, 1e-12), err_proj);
  if (err_proj > 1e-12) {
    fprintf(stderr, "projectionOnConeVector correctness check failed!\n");
    info += 1;
  }

  /* ---------- naturalMapOnConeVector ---------- */
  t0 = now();
  for (int rep = 0; rep < nrep; ++rep) {
    for (int i = 0; i < nc; ++i) {
      naturalMapOnCone(&r[3 * i], &u[3 * i], rho, mu[i], &F_loop[3 * i]);
    }
  }
  double t_loop_nat = (now() - t0) / nrep;

  t0 = now();
  for (int rep = 0; rep < nrep; ++rep) {
    naturalMapOnConeVector(nc, r, u, rho, mu, F_vec, work);
  }
  double t_vec_nat = (now() - t0) / nrep;

  double err_nat = max_abs_diff(F_loop, F_vec, n);
  printf(
      "naturalMapOnConeVector: nc=%d  cone-wise=%.6f s  vectorized=%.6f s  speedup=%.2fx  "
      "err=%e\n",
      nc, t_loop_nat, t_vec_nat, t_loop_nat / fmax(t_vec_nat, 1e-12), err_nat);
  if (err_nat > 1e-12) {
    fprintf(stderr, "naturalMapOnConeVector correctness check failed!\n");
    info += 1;
  }

  free(r);
  free(u);
  free(mu);
  free(work);
  free(r_loop);
  free(r_vec);
  free(F_loop);
  free(F_vec);

  return info;
}
