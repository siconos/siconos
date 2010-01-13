/* Siconos-Numerics, Copyright INRIA 2005-2010.
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

#ifndef NSSTOOLS_H
#define NSSTOOLS_H

/*!\file NSSTools.h
  Header to collect basic tools, structures definition or any usefull things for NSSpack

*/

#ifdef __cplusplus
extern "C" {
#endif

  /** Search for the max. element of a vector
      \param[in] x, the vector
      \param[in-out] solution, value of the greatest element of x
      \param[in] n, size of x
  */
  void max_part(double*, double*, int);

  /** compare two double a and b, and return the max.
   *  \param a, double*
   *  \param b, double*
   *  \param c, double*, the max
   */
  void maxf(double*, double*, double*);

  /** Search for the min. element of a vector
      \param[in] x, the vector
      \param[in-out] solution, value of the smallest element of x
      \param[in] n, size of x
  */
  void min_part(double*, double*, int);

  /** compare two double a and b, and return the min.
   *  \param a, double*
   *  \param b, double*
   *  \param c, double*, the min
   */
  void minf(double*, double*, double*);

  /** Positive part values of the components of a vector
      \param[in] x, the vector
      \param[in-out] solution, vector of positive part values of x components
      \param[in] n, size of x
  */
  void pos_part(double*, double*, int);

  /** Absolute values of the components of a vector
      \param[in] x, the vector
      \param[in-out] solution, vector of absolute values of x components
      \param[in] n, size of x
  */
  void abs_part(double*, double*, int);

  /**
      Input na, a, nb, b
      Output nc, c
      a and b: interger vectors in increasing order
      c : vector of integers of a that are not in b.
      \author Nineb Sheherazade & Dureisseix David.
  */
  void diffns(int *na, int *a, int *nb, int * b, int *nc, int *c);

  /** */
  void sortsn_(int*ddl_i, int *sort, int *n);

#ifdef __cplusplus
}
#endif

#endif
