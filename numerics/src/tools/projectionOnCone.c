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
#include "projectionOnCone.h"

#include "op3x3.h"   // for SET3, eig_3x3
                     //
#include <math.h>    // for sqrt
#include <stdio.h>   // for fprintf, stderr
#include <stdlib.h>  // for exit, EXIT_FAILURE

/* #define DEBUG_NOCOLOR */
/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "siconos_debug.h"  // for DEBUG_PRINTF

unsigned projectionOnCone(double* r, double mu) {
  double normT = sqrt(r[1] * r[1] + r[2] * r[2]);
  /* hypot of libm is sure but really slow */
  /* double normT = hypot(r[1], r[2]); */
  if (mu * normT <= -r[0]) {
    r[0] = 0.0;
    r[1] = 0.0;
    r[2] = 0.0;
    return PROJCONE_DUAL;
  } else if (normT <= mu * r[0]) {
    return PROJCONE_INSIDE;
  } else {
    double mu2 = mu * mu;
    r[0] = (mu * normT + r[0]) / (mu2 + 1.0);
    r[1] = mu * r[0] * r[1] / normT;
    r[2] = mu * r[0] * r[2] / normT;
    return PROJCONE_BOUNDARY;
  }
}
unsigned projectionOnDualCone(double* u, double mu) {
  double normT = sqrt(u[1] * u[1] + u[2] * u[2]);
  /* hypot of libm is sure but really slow */
  /* double normT = hypot(u[1], u[2]); */

  if (normT <= -mu * u[0]) {
    u[0] = 0.0;
    u[1] = 0.0;
    u[2] = 0.0;
    return PROJCONE_DUAL;
  } else if (mu * normT <= u[0]) {
    return PROJCONE_INSIDE;
  } else {
    double mu2 = mu * mu;
    u[0] = (normT + mu * u[0]) / (mu2 + 1.0);
    u[1] = u[0] * u[1] / normT;
    u[2] = u[0] * u[2] / normT;
    u[0] = mu * u[0];

    return PROJCONE_BOUNDARY;
  }
  /* return projectionOnCone(u, 1.0/mu); */
}

void projectionOnSecondOrderCone(double* r, double mu, int size) {
  if (size == 3) {
    projectionOnCone(r, mu);
  } else {
    fprintf(stderr,
            "Numerics, projectionOnSecondOrderCone f not yet implemented for size != 3 \n");
    exit(EXIT_FAILURE);
  }
}

unsigned subdifferentialProjectionOnCone(double* H, double* r, double mu) {
  double normT = sqrt(r[1] * r[1] + r[2] * r[2]);
  /* hypot of libm is sure but really slow */
  /* double normT = hypot(r[1], r[2]); */
  SET3X3(H);
  zero3x3(H00);

  if (mu * normT <= -r[0]) {
    // printf("We are in the polar cone\n");
    return PROJCONE_DUAL;
  } else if (normT <= mu * r[0]) {
    *H00 = 1.;
    *H11 = 1.;
    *H22 = 1.;
    // printf("We are in the cone\n");
    return PROJCONE_INSIDE;
  } else {
    double oneoveroneplusmu2 = 1. / (1. + mu * mu);

    *H00 = 1. * oneoveroneplusmu2;
    *H10 = mu * r[1] * oneoveroneplusmu2;
    *H20 = mu * r[2] * oneoveroneplusmu2;
    double s1 = r[1] / normT;
    double s2 = r[2] / normT;

    double muoveroneplusmu2 = mu * oneoveroneplusmu2;

    *H01 = muoveroneplusmu2 * s1;
    *H02 = muoveroneplusmu2 * s2;

    *H11 =
        muoveroneplusmu2 / normT * (r[0] + mu * normT + r[0] * r[1] * r[1] / (normT * normT));
    *H12 = muoveroneplusmu2 / normT * (r[0] * r[1] * r[2] / (normT * normT));
    *H21 = muoveroneplusmu2 / normT * (r[0] * r[1] * r[2] / (normT * normT));
    *H22 =
        muoveroneplusmu2 / normT * (r[0] + mu * normT + r[0] * r[2] * r[2] / (normT * normT));
    // printf("We are outside the cone and its polar\n");
    return PROJCONE_BOUNDARY;
  }
}
