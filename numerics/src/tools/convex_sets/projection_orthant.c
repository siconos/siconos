// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "projection_orthant.h"

void orthant_projection(double* x, int n) {
  for (int i = 0; i < n; i++)
    if (x[i] < 0.0) x[i] = 0.0;
}

void orthant_projection_jacobian_diag(double* d, const double* x, int n) {
  for (int i = 0; i < n; i++) d[i] = (x[i] > 0.0) ? 1.0 : 0.0;
}

void box_projection(double* x, int n, const double* lb, const double* ub) {
  for (int i = 0; i < n; ++i) {
    if (x[i] < lb[i])
      x[i] = lb[i];
    else if (x[i] > ub[i])
      x[i] = ub[i];
  }
}

void box_projection_jacobian_diag(double* d, const double* x, int n, const double* lb,
                                  const double* ub) {
  for (int i = 0; i < n; i++) d[i] = (x[i] > lb[i] && x[i] < ub[i]) ? 1.0 : 0.0;
}
