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
#include "projectionOnCone.h"

#include "op3x3.h"   // for SET3, eig_3x3
                     //
#include "SiconosBlas.h"
#include <math.h>    // for sqrt
#include <stdio.h>   // for fprintf, stderr
#include <stdlib.h>  // for exit, EXIT_FAILURE, malloc, free

/* #define DEBUG_NOCOLOR */
/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "numerics_errors.h"
#include "siconos_debug.h"  // for DEBUG_PRINTF

unsigned projectionOnCone(double* r, double mu) {
  double normT = sqrt(r[1] * r[1] + r[2] * r[2]);
  /* hypot of libm is sure but really slow */
  /* double normT = hypot(r[1], r[2]); */
  if (mu * normT <= -r[0]) {
    r[0] = 0.0;
    r[1] = 0.0;
    r[2] = 0.0;
    return PROJCONE_DUAL;
  } else if (normT <= mu * r[0]) {
    return PROJCONE_INSIDE;
  } else {
    double mu2 = mu * mu;
    r[0] = (mu * normT + r[0]) / (mu2 + 1.0);
    r[1] = mu * r[0] * r[1] / normT;
    r[2] = mu * r[0] * r[2] / normT;
    return PROJCONE_BOUNDARY;
  }
}
void naturalMapOnCone(const double r[3], const double u[3], double rho, double mu, double F[3]) {
  double w[3];
  w[0] = r[0] - rho * u[0];
  w[1] = r[1] - rho * u[1];
  w[2] = r[2] - rho * u[2];

  projectionOnCone(w, mu);

  F[0] = r[0] - w[0];
  F[1] = r[1] - w[1];
  F[2] = r[2] - w[2];
}

void projectionOnConeVector(int nc, double* r, const double* mu, double* work) {
  int alloc_work = 0;
  if (work == NULL) {
    work = (double*)malloc(nc * sizeof(double));
    alloc_work = 1;
  }

  /* Tangential norms of all contacts at once.
     A plain loop is faster here than per-contact BLAS dnrm2 calls because
     the blocks are short (length 2) and contiguous in memory. */
  for (int i = 0; i < nc; ++i) {
    const double t1 = r[3 * i + 1];
    const double t2 = r[3 * i + 2];
    work[i] = sqrt(t1 * t1 + t2 * t2);
  }

  for (int i = 0; i < nc; ++i) {
    const int pos = 3 * i;
    const double rN = r[pos];
    const double normT = work[i];
    const double mui = mu[i];

    if (mui * normT <= -rN) {
      r[pos] = 0.0;
      r[pos + 1] = 0.0;
      r[pos + 2] = 0.0;
    } else if (normT <= mui * rN) {
      /* projection is identity, leave r unchanged */
    } else {
      const double mu2 = mui * mui;
      const double rN_proj = (mui * normT + rN) / (mu2 + 1.0);
      const double scale = mui * rN_proj / normT;
      r[pos] = rN_proj;
      r[pos + 1] = scale * r[pos + 1];
      r[pos + 2] = scale * r[pos + 2];
    }
  }

  if (alloc_work) free(work);
}

void naturalMapOnConeVector(int nc, const double* r, const double* u, double rho,
                            const double* mu, double* F, double* work) {
  const int n = 3 * nc;
  int alloc_work = 0;
  if (work == NULL) {
    work = (double*)malloc(nc * sizeof(double));
    alloc_work = 1;
  }

  /* F = r - rho * u (then used as the argument of the projection) */
  cblas_dcopy(n, r, 1, F, 1);
  cblas_daxpy(n, -rho, u, 1, F, 1);

  /* Reuse the vectorized projection on w = r - rho*u */
  projectionOnConeVector(nc, F, mu, work);

  /* Natural map residual: r - projection(w) */
  for (int i = 0; i < n; ++i) {
    F[i] = r[i] - F[i];
  }

  if (alloc_work) free(work);
}

unsigned projectionOnDualCone(double* u, double mu) {
  double normT = sqrt(u[1] * u[1] + u[2] * u[2]);
  /* hypot of libm is sure but really slow */
  /* double normT = hypot(u[1], u[2]); */

  if (normT <= -mu * u[0]) {
    u[0] = 0.0;
    u[1] = 0.0;
    u[2] = 0.0;
    return PROJCONE_DUAL;
  } else if (mu * normT <= u[0]) {
    return PROJCONE_INSIDE;
  } else {
    double mu2 = mu * mu;
    u[0] = (normT + mu * u[0]) / (mu2 + 1.0);
    u[1] = u[0] * u[1] / normT;
    u[2] = u[0] * u[2] / normT;
    u[0] = mu * u[0];

    return PROJCONE_BOUNDARY;
  }
  /* return projectionOnCone(u, 1.0/mu); */
}

void projectionOnSecondOrderCone(double* r, double mu, int size) {
  if (size == 3) {
    projectionOnCone(r, mu);
  } else {
    fprintf(stderr,
            "Numerics, projectionOnSecondOrderCone not yet implemented for size != 3\n");
    exit(EXIT_FAILURE);
  }
}

unsigned subdifferentialProjectionOnCone(double* H, double* r, double mu) {
  double normT = sqrt(r[1] * r[1] + r[2] * r[2]);
  /* hypot of libm is sure but really slow */
  /* double normT = hypot(r[1], r[2]); */
  SET3X3(H);
  zero3x3(H00);

  if (mu * normT <= -r[0]) {
    // printf("We are in the polar cone\n");
    return PROJCONE_DUAL;
  } else if (normT <= mu * r[0]) {
    *H00 = 1.;
    *H11 = 1.;
    *H22 = 1.;
    // printf("We are in the cone\n");
    return PROJCONE_INSIDE;
  } else {
    double oneoveroneplusmu2 = 1. / (1. + mu * mu);

    double s1 = r[1] / normT;
    double s2 = r[2] / normT;
    *H00 = 1. * oneoveroneplusmu2;

    double muoveroneplusmu2 = mu * oneoveroneplusmu2;
    *H10 = muoveroneplusmu2 * s1;
    *H20 = muoveroneplusmu2 * s2;

    *H01 = *H10;
    *H02 = *H20;

    double muoveroneplusmu2overnormT = muoveroneplusmu2 / normT;
    double r0plusmunormT = r[0] + mu * normT;

    *H11 = muoveroneplusmu2overnormT * (r0plusmunormT - r[0] * s1 * s1);
    *H12 = muoveroneplusmu2overnormT * (-r[0] * s1 * s2);
    *H21 = muoveroneplusmu2overnormT * (-r[0] * s1 * s2);
    *H22 = muoveroneplusmu2overnormT * (r0plusmunormT - r[0] * s2 * s2);
    // printf("We are outside the cone and its polar\n");
    return PROJCONE_BOUNDARY;
  }
}
