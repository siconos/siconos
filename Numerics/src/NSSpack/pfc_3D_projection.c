
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
/*!\file pfc_3D_nlgs.c
 *
 * This subroutine allows the primal resolution of contact problems with friction.\n
 *
 *   Try \f$(z,w)\f$ such that:\n
 *   \f$
 *    \left\lbrace
 *     \begin{array}{l}
 *      M z + q = w\\
 *      0 \le z_n \perp w_n \ge 0\\
 *      -w_t \in \partial\psi_{[-\mu z_n, \mu z_n]}(z_t)\\
 *     \end{array}
 *    \right.
 *   \f$
 *
 * here M is an n by n  matrix, q an n-dimensional vector, z an n-dimensional  vector and w an n-dimensional vector.
 *
 * \fn  pfc_3D_projection( int *nn , double *vec , double *q , double *z , double *w , double coef )
 *
 * Generic pfc_3D parameters:\n
 *
 * \param nn      Unchanged parameter which represents the dimension of the system.
 * \param vec     Unchanged parameter which contains the components of the matrix with a fortran storage.
 * \param q       Unchanged parameter which contains the components of the right hand side vector.
 * \param z       Modified parameter which contains the initial solution and returns the solution of the problem.
 * \param w       Modified parameter which returns the solution of the problem.
 * \param coef    Unchanged parameter which represents the friction coefficient
 *
 *
 * \author houari khenous 14/09/2007 .
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "LA.h"
#include <time.h>


void pfc_3D_projection(int n , double *C , double *b , double *zz , double *ww , double coef , void (*Compute_G), void (*Compute_JacG), int *iparam_local , double *dparam_local)
{

  double mrn, num, coef2;

  coef2 = coef * coef;

  if (b[0] > 0.)
  {
    zz[0] = 0.;
    zz[1] = 0.;
    zz[2] = 0.;
  }
  else
  {
    zz[0] = -b[0] / C[0 * n + 0];
    zz[1] = -b[1] / C[1 * n + 1];
    zz[2] = -b[2] / C[2 * n + 2];

    mrn = zz[1] * zz[1] + zz[2] * zz[2];

    if (mrn > coef2 * zz[0]*zz[0])
    {
      num = coef * zz[0] / sqrt(mrn);
      zz[1] = zz[1] * num;
      zz[2] = zz[2] * num;
    }
  }
}

