/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.

 * Copyright 2016 INRIA.

 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at

 * http://www.apache.org/licenses/LICENSE-2.0

 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

#ifndef QUARTIC_H
#define QUARTIC_H

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

#include "SiconosConfig.h"

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif
  /** Suite of procedures for finding the (complex) roots of the quadratic,
  \param p Coefficients of the polynomial
  \param r root of the polynomial r[1][k] real part of the kth root r[2][k] imaginary part.
  */
  int QUADROOTS(double  p[5], double r[3][5]);
  int CUBICROOTS(double p[5], double r[3][5]);
  int BIQUADROOTS(double p[5], double r[3][5]);


#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif
#endif
