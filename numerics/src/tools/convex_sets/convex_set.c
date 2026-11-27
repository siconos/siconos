// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "convex_set.h"

#include <stdlib.h>  // for malloc, free

/* blocks up to this size use a buffer on the stack */
#define CONVEX_SET_STACK_DIM 16

double convex_set_natural_map_block(const ConvexSet* K, int i, const double* r,
                                    const double* v, double* F) {
  const int dim = K->dim;
  double buffer[CONVEX_SET_STACK_DIM];
  double* P = (dim <= CONVEX_SET_STACK_DIM) ? buffer : (double*)malloc(dim * sizeof(double));
  for (int k = 0; k < dim; k++) P[k] = r[k] - v[k];
  convex_set_project_block(K, i, P);
  double F2 = 0.0;
  for (int k = 0; k < dim; k++) {
    const double f = r[k] - P[k];
    if (F) F[k] = f;
    F2 += f * f;
  }
  if (P != buffer) free(P);
  return F2;
}

void convex_set_project(const ConvexSet* K, double* x) {
  const int dim = K->dim;
  switch (K->type) {
    case CONVEX_SET_FREE:
      return;
    case CONVEX_SET_ORTHANT:
      orthant_projection(x, dim * K->nblocks);
      return;
    case CONVEX_SET_BOX:
      box_projection(x, dim * K->nblocks, K->p1, K->p2);
      return;
    case CONVEX_SET_SOC:
      if (dim == 3) {
        soc_projection_3d_stack(K->nblocks, x, K->p1, NULL);
        return;
      }
      /* fall through */
    default:
      for (int i = 0; i < K->nblocks; i++) convex_set_project_block(K, i, x + i * dim);
  }
}

int convex_set_projection_jacobian_block(const ConvexSet* K, int i, double* J,
                                         const double* x) {
  const int dim = K->dim;
  switch (K->type) {
    case CONVEX_SET_SOC:
      if (dim != 3) return 1;
      soc_projection_jacobian_3d(J, x, K->p1[i]);
      return 0;
    case CONVEX_SET_BALL:
      ball_projection_jacobian(J, x, dim, K->p1[i]);
      return 0;
    case CONVEX_SET_CYLINDER:
      cylinder_projection_jacobian(J, x, dim, K->p1[i]);
      return 0;
    case CONVEX_SET_ROLLING_CONE:
      if (dim != 5) return 1;
      rolling_cone_projection_jacobian(J, (double*)x, K->p1[i], K->p2[i]);
      return 0;
    default:
      break;
  }
  /* diagonal Jacobians */
  double buffer[CONVEX_SET_STACK_DIM];
  double* d = (dim <= CONVEX_SET_STACK_DIM) ? buffer : (double*)malloc(dim * sizeof(double));
  int info = 0;
  switch (K->type) {
    case CONVEX_SET_FREE:
      for (int k = 0; k < dim; k++) d[k] = 1.0;
      break;
    case CONVEX_SET_ORTHANT:
      orthant_projection_jacobian_diag(d, x, dim);
      break;
    case CONVEX_SET_BOX:
      box_projection_jacobian_diag(d, x, dim, K->p1 + i * dim, K->p2 + i * dim);
      break;
    case CONVEX_SET_POSITIVE_CYLINDER: {
      /* product R+ x B: diag(1_{x_N > 0}, Jacobian of the ball projection) */
      double* JT = (double*)malloc((dim - 1) * (dim - 1) * sizeof(double));
      ball_projection_jacobian(JT, x + 1, dim - 1, K->p1[i]);
      for (int k = 0; k < dim * dim; k++) J[k] = 0.0;
      J[0] = (x[0] > 0.0) ? 1.0 : 0.0;
      for (int a = 1; a < dim; a++)
        for (int b = 1; b < dim; b++) J[a + dim * b] = JT[(a - 1) + (dim - 1) * (b - 1)];
      free(JT);
      if (d != buffer) free(d);
      return 0;
    }
    default:
      info = 1;
  }
  if (!info) {
    for (int k = 0; k < dim * dim; k++) J[k] = 0.0;
    for (int k = 0; k < dim; k++) J[k + dim * k] = d[k];
  }
  if (d != buffer) free(d);
  return info;
}

double convex_set_natural_map(const ConvexSet* K, const double* r, const double* v,
                              double* F) {
  const int n = K->dim * K->nblocks;
  for (int k = 0; k < n; k++) F[k] = r[k] - v[k];
  convex_set_project(K, F);
  double F2 = 0.0;
  for (int k = 0; k < n; k++) {
    F[k] = r[k] - F[k];
    F2 += F[k] * F[k];
  }
  return F2;
}
