// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef ROLLING_FRICTION_naturalmap_functions_H
#define ROLLING_FRICTION_naturalmap_functions_H

/*!\file rolling_naturalmap_functions.h

  Typedef and functions declarations related to naturalmap map
  formulation for rolling friction contact problems.


 */
#if defined(__cplusplus)
extern "C" {
#endif

void rolling_friction_3D_computeNaturalMap(double reaction[5], double velocity[5], double mu,
                                           double mur, double rho[1], double result[5],
                                           double A[25], double B[25]);

#if defined(__cplusplus)
}
#endif

#endif
