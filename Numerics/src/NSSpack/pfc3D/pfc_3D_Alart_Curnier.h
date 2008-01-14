/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2007.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */
/*!\file pfc_3D_Alart_Curnier.h
 *
 *
 *
 *
 *
 * Computation of the equation to be solved (We try to solve G = 0) with Newton method. \n
 *
 *            |          w - M z - q             |\n
 *            |                                  |\n
 *       G =  | 1/rn*[wn - pos_part(wn - rn*zn)] |;\n
 *            |                                  |\n
 *            |   1/rt*[wT - proj(wT - rt*zT)]   |\n
 *
 *
 *
 * where M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional  vector and w an  * n-dimensional vector.\n
 *
 * Inside, we compute G, JacG and the inverse of JacG. \n
 *
 * \fn  void Compute_G
 *
 * \fn  void Compute_JacG
 *
 * \fn  void matrix_inv3
 *
 * \fn void Linesearch_AC
 *
 *
 * \author Houari Khenous  last modification (08/10/2007)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "LA.h"
#include <time.h>

#ifdef __cplusplus
extern "C" {
#endif

  /** Compute function G */
  void Compute_G_AC(int m, double *G, double *x , double *C, double *y, double *b, double *param1, double *param2, double *param3, double rn, double rt, double coef);

  /** Compute Jacobian of function G */
  void Compute_JacG_AC(int m, double *JacG , double *x ,  double *C , double *y, double *b , double *param1, double *param2, double *param3, double rn, double rt, double coef);

  /**  Inverse Matrix 3x3  */
  void matrix_inv3(double *a, double *b);

  /** Linesearch */
  void Linesearch_AC(int n, double *zz, double *ww, double *www, double *b, double *C, double *param1, double *param2, double *param3, double an, double at, double mu, double err1);

#ifdef __cplusplus
}
#endif
