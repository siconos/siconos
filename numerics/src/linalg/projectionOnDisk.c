// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "projectionOnDisk.h"

#include <math.h>

#include "numerics_errors.h"

void projectionOnDisk(double* r, double R) {
  double normTsquare = r[0] * r[0] + r[1] * r[1];
  if (normTsquare <= R * R) {
    return;
  } else {
    normTsquare = sqrt(normTsquare);
    r[0] = R * r[0] / normTsquare;
    r[1] = R * r[1] / normTsquare;
    return;
  }
}

void projectionOnGeneralDisk(double* r, double R, int dim) {}
