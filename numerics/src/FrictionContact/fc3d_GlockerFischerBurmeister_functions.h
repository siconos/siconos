// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef FRICTIONCONTACT3D_GlockerFischerBurmeister_functions_H
#define FRICTIONCONTACT3D_GlockerFischerBurmeister_functions_H

/*!\file fc3d_GlockerFischerBurmeister_functions.h

  \brief Typedef and functions declarations related GlockerFischerBurmeister
 */
#if defined(__cplusplus)
extern "C" {
#endif
typedef void (*UpdateSolverPtr)(int, double*);

void F_GlockerFischerBurmeister(int sizeF, double* reaction, double* FVector, int up2Date);

void jacobianF_GlockerFischerBurmeister(int sizeF, double* reaction, double* jacobianFMatrix,
                                        int up2Date);

#if defined(__cplusplus)
}
#endif

#endif
