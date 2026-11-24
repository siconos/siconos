// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "NSSTools.h"  // for pos_part
#include "numerics_errors.h"
void pos_part(unsigned n, double* restrict x, double* restrict x_plus) {
  for (unsigned i = 0; i < n; ++i) {
    if (x[i] > 0.0) {
      x_plus[i] = x[i];
    } else {
      x_plus[i] = 0.0;
    }
  }
}
