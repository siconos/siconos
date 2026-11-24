// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef PLASTICITY2D_naturalmap_functions_H
#define PLASTICITY2D_naturalmap_functions_H

/*!\file fc3d_AlartCurnier_functions.h

  Typedef and functions declarations related to naturalmap map
  formulation for Mohr Coulomb 2D contact problems.


 */

#if defined(__cplusplus)
extern "C" {
#endif

void plasticity_2d_computeNaturalMap(double reaction[3], double velocity[3], double eta,
                                     double theta, double rho[3], double result[3],
                                     double A[9], double B[9]);

#if defined(__cplusplus)
}
#endif

#endif
