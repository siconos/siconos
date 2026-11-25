// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "NSSTools.h"  // for max_part
#include "numerics_errors.h"
void max_part(double* x, double* sol, int n)

{
  int i;
  double max;
  max = x[0];

  for (i = 1; i < n; i++) {
    if (max < x[i]) max = x[i];
  }

  *sol = max;
}
