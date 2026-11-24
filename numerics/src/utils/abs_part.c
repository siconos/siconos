// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "NSSTools.h"  // for abs_part
#include "numerics_errors.h"

void abs_part(double x[], double sol[], int n)

{
  int i;
  double eps;

  eps = 1.e-12;

  for (i = 0; i < n; i++) {
    if (x[i] >= eps) {
      sol[i] = x[i];
    } else {
      sol[i] = -x[i];
    }
  }
}
