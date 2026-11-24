// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*!\file pinv.h
 * \brief peseudo-inverse computations*/

#ifndef Pinv_H
#define Pinv_H

#if defined(__cplusplus)
extern "C" {
#endif

/** Compute the pseudo-inverse of dense matrix with column major storage
\param A the matrix to be inversed
\param n the number of rows of A
\param m the number of columns of A
\param tolerance threshold used to validate the computation: if the error is less than this
value, the computation is considered valid \return the conditioning
*/
double pinv(double* A, int n, int m, double tolerance);

#if defined(__cplusplus)
}
#endif
#endif
