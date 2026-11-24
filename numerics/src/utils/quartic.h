/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef QUARTIC_H
#define QUARTIC_H

/*!\file quartic.h
 * \brief quartic algorithm to find roots of quartic polynomials*/

/** CACM Algorithm 326
   Roots of low order polynomials
   Author: Terence R.F.Nonweiler
   CACM  (Apr 1968) p269
   Translated into c and programmed by M.Dow
   ANUSF, Australian National University, Canberra, Australia
   m.dow@anu.edu.au

Suite of procedures for finding the (complex) roots of the
quadratic, cubic or quartic polynomials by explicit algebraic methods.
Each Returns
x=r[1][k] + i r[2][k]  k=1,...,n, where n={2,3,4}
as roots of
sum_{k=0:n} p[k] x^(n-k) =0
Assume p[0]<>0 (overflows otherwise)
**/

#if defined(__cplusplus)
extern "C" {
#endif
/** Suite of procedures for finding the (complex) roots of the quadratic,
\param p Coefficients of the polynomial
\param r root of the polynomial r[1][k] real part of the kth root r[2][k] imaginary part.
\return always 0
*/
int QUADROOTS(double p[5], double r[3][5]);
int CUBICROOTS(double p[5], double r[3][5]);
int BIQUADROOTS(double p[5], double r[3][5]);

#if defined(__cplusplus)
}
#endif
#endif
