// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "projection_ball.h"

#include <math.h>  // for sqrt

static inline double squared_norm(const double* x, int n) {
  double s = 0.0;
  for (int i = 0; i < n; i++) s += x[i] * x[i];
  return s;
}

unsigned ball_projection(double* x, int dim, double R) {
  double n = sqrt(squared_norm(x, dim));
  if (n <= R) return BALL_PROJECTION_INSIDE;
  double scale = (R > 0.0) ? R / n : 0.0;
  for (int i = 0; i < dim; i++) x[i] *= scale;
  return BALL_PROJECTION_BOUNDARY;
}

unsigned ball_projection_jacobian(double* J, const double* x, int dim, double R) {
  double n = sqrt(squared_norm(x, dim));
  for (int i = 0; i < dim * dim; i++) J[i] = 0.0;
  if (n <= R) {
    for (int i = 0; i < dim; i++) J[i + dim * i] = 1.0;
    return BALL_PROJECTION_INSIDE;
  }
  /* (R/||x||) (I - u u^T), u = x/||x|| */
  double f = R / n;
  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++)
      J[i + dim * j] = f * ((i == j ? 1.0 : 0.0) - x[i] * x[j] / (n * n));
  return BALL_PROJECTION_BOUNDARY;
}

unsigned cylinder_projection(double* x, int dim, double R) {
  return ball_projection(x + 1, dim - 1, R);
}

unsigned cylinder_projection_jacobian(double* J, const double* x, int dim, double R) {
  for (int i = 0; i < dim * dim; i++) J[i] = 0.0;
  J[0] = 1.0;
  int m = dim - 1;
  double nT = sqrt(squared_norm(x + 1, m));
  if (nT <= R) {
    for (int i = 1; i < dim; i++) J[i + dim * i] = 1.0;
    return BALL_PROJECTION_INSIDE;
  }
  double f = R / nT;
  for (int i = 1; i < dim; i++)
    for (int j = 1; j < dim; j++)
      J[i + dim * j] = f * ((i == j ? 1.0 : 0.0) - (x[i] / nT) * (x[j] / nT));
  return BALL_PROJECTION_BOUNDARY;
}

unsigned positive_cylinder_projection(double* x, int dim, double R) {
  unsigned status = ball_projection(x + 1, dim - 1, R);
  if (x[0] < 0.0) {
    x[0] = 0.0;
    status = BALL_PROJECTION_BOUNDARY;
  }
  return status;
}
