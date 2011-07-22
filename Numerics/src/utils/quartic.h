/* Siconos-Numerics, Copyright INRIA 2005-2011.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
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

#ifdef __cplusplus
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


#ifdef __cplusplus
}
#endif
#endif
