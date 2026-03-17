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
#include "plasticity_2d_vonmises_projection.h"

#include <math.h>
#include <stdio.h>

/* #define DEBUG_NOCOLOR */
/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "siconos_debug.h"

double plasticity_2d_vonMises_equivalent_stress(const double stress[3]) {
  double sx = stress[0];
  double sy = stress[1];
  double txy = stress[2];
  
  /* Von Mises equivalent stress: q = sqrt(sx^2 + sy^2 - sx*sy + 3*txy^2) */
  double q = sqrt(sx * sx + sy * sy - sx * sy + 3.0 * txy * txy);
  
  return q;
}

int plasticity_2d_vonMises_check_yield(const double stress[3], double sigma_y) {
  double q = plasticity_2d_vonMises_equivalent_stress(stress);
  
  DEBUG_PRINTF("Von Mises check: q = %e, sigma_y = %e\n", q, sigma_y);
  
  return (q <= sigma_y) ? 1 : 0;
}

int plasticity_2d_projectionOnVonMises(double stress[3], double sigma_y) {
  double q = plasticity_2d_vonMises_equivalent_stress(stress);
  
  DEBUG_PRINTF("Projection: q = %e, sigma_y = %e\n", q, sigma_y);
  
  /* If stress is inside yield surface, no projection needed */
  if (q <= sigma_y) {
    DEBUG_PRINT("Von Mises: stress inside yield surface, no projection\n");
    return 1;
  }
  
  /* Radial return: scale stress to yield surface */
  double scale = sigma_y / q;
  
  DEBUG_PRINTF("Von Mises: scaling by %e\n", scale);
  
  stress[0] *= scale;
  stress[1] *= scale;
  stress[2] *= scale;
  
  return 0;
}
