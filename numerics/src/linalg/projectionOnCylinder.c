/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#include "projectionOnCylinder.h"

#include <math.h>

#include "numerics_errors.h"

void projectionOnCylinder(double* r, double R) {
  double normTsquare = r[1] * r[1] + r[2] * r[2];

  if (r[0] >= 0) {
    if (normTsquare <= R * R) {
      return;
    } else {
      normTsquare = sqrt(normTsquare);
      r[1] = R * r[1] / normTsquare;
      r[2] = R * r[2] / normTsquare;
      return;
    }
  } else {
    r[0] = 0.0;
    if (0 < normTsquare) {
      normTsquare = sqrt(normTsquare);
      r[1] = R * r[1] / normTsquare;
      r[2] = R * r[2] / normTsquare;
      /*    r[1]=0.0; */
      /*    r[2]=0.0; */
      return;
    } else {
      r[1] = 0.0;
      r[2] = 0.0;
      return;
    }
  }
}
void projectionOnGeneralCylinder(double* r, double R, int dim) {}
