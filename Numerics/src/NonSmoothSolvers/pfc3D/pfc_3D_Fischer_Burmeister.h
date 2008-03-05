/* Siconos-Numerics version 3.0.0, Copyright INRIA 2005-2008.
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
/*!\file pfc_3D_Fischer_Burmeister.h
 *
 *
 *   void Compute_G_FB
 *
 *   void Compute_JacG_FB
 *
 *   void matrix_inv2
 *
 *   void Linesearch_FB
 *
 *
 * \author Houari Khenous last modification (08/10/2007)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "LA.h"
#include <time.h>

/** Compute function G */
void Compute_G_FB(int m, double *G , double *y , double *C , double *x , double *b, double *Ip , double *IP, double *I3, double coef3 , double coef21, double coef);

/** Compute Jacobian of function G */
void Compute_JacG_FB(int m, double *JacG , double *y , double *C , double *x , double *b, double *Ip , double *IP, double *I3, double coef3 , double coef21, double coef);

/**  Inverse Matrix 2x2  */
void matrix_inv2(double *a, double *b);

/** Lineserach */
void Linesearch_FB(int n, double *zz, double *ww, double *www, double *b, double *C, double *Ip, double *IP, double *I3, double coef21, double coef3, double mu, double err1);
