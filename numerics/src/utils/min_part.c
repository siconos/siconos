/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#include "NSSTools.h"  // for min_part
#include "numerics_errors.h"

void min_part(double* x, double* sol, int n)

{
  int i;
  double minx;

  minx = x[0];

  for (i = 1; i < n; i++) {
    if (minx > x[i]) minx = x[i];
  }

  *sol = minx;
}
