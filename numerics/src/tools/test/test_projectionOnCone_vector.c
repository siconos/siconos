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
  const int nc = 5;
  const int n = 3 * nc;

  /* Stacked reactions covering the three possible projection statuses.
     Contacts are chosen so that the vectorized code exercises inside,
     dual and boundary branches with different mu values. */
  double r_ori[] = {
      /* contact 0: inside */
      1.0, 0.05, 0.0,
      /* contact 1: dual */
      -1.0, 0.5, 0.0,
      /* contact 2: boundary */
      0.2, 1.0, 0.0,
      /* contact 3: inside, nonzero tangential component */
      2.0, 0.2, 0.3,
      /* contact 4: boundary, large normal */
      5.0, 3.0, 4.0};
  double mu[] = {0.1, 0.3, 0.5, 0.7, 0.8};

  double r_loop[3 * 5];
  double r_vec[3 * 5];
  double work[5];

  memcpy(r_loop, r_ori, sizeof(r_loop));
  memcpy(r_vec, r_ori, sizeof(r_vec));

  /* reference: cone-wise projection */
  for (int i = 0; i < nc; ++i) {
    projectionOnCone(&r_loop[3 * i], mu[i]);
  }

  /* vectorized projection with caller-provided workspace */
  projectionOnConeVector(nc, r_vec, mu, work);
  double err = max_abs_diff(r_loop, r_vec, n);
  printf("projectionOnConeVector (with work) max diff = %e\n", err);
  if (err > 1e-14) {
    fprintf(stderr, "Vectorized projection differs from cone-wise version!\n");
    info += 1;
  }

  /* vectorized projection with internal allocation */
  memcpy(r_vec, r_ori, sizeof(r_vec));
  projectionOnConeVector(nc, r_vec, mu, NULL);
  err = max_abs_diff(r_loop, r_vec, n);
  printf("projectionOnConeVector (work=NULL) max diff = %e\n", err);
  if (err > 1e-14) {
    fprintf(stderr, "Vectorized projection (NULL work) differs from cone-wise version!\n");
    info += 1;
  }

  return info;
}
