// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef NONSMOOTHNEWTONNEIGH_H
#define NONSMOOTHNEWTONNEIGH_H

#include "NonSmoothNewton.h"  // for NewtonFunctionPtr

/*!\file NonSmoothNewtonNeighbour.h
  Typedef and functions declarations related to non-smooth Newton solver

  Solve \f$ \phi(z) = 0 \f$ using a Newton method.

  The algorithm is alg 4.1 of the paper of Kanzow and Kleinmichel, "A new class of semismooth
  Newton-type methods for nonlinear complementarity problems", in Computational Optimization
  and Applications, 11, 227-251 (1998).

 */

#if defined(__cplusplus)
extern "C" {
#endif

double* nonSmoothNewtonNeighInitMemory(int n, double* dWork, int* iWork);

int nonSmoothNewtonNeigh(int, double*, NewtonFunctionPtr*, NewtonFunctionPtr*, int*, double*);

int nonSmoothNewtonNeigh_getNbIWork(int n, int m);
int nonSmoothNewtonNeigh_getNbDWork(int n, int m);

/*only for debug*/
void NSNN_thisIsTheSolution(int n, double* z);
void NSNN_reset(void);
#if defined(__cplusplus)
}
#endif

#endif
