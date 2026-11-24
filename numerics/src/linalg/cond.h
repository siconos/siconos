// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef Cond_H
#define Cond_H

/*!\file cond.h
 * \brief condition number computation
 */

#include <stddef.h>  // size_t

#ifdef __cplusplus
extern "C" {
#endif
/** Compute the conditioning of dense matrix with column major storage
\param A the matrix to be read
\param n the number of rows of A
\param m the number of columns of A
\return the condition number
*/
double cond(double* A, size_t n, size_t m);

#ifdef __cplusplus
}
#endif
#endif
