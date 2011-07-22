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

#include "SparseBlockMatrix.h"
#include "LA.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void NCP_block_compute_error(int n, SparseBlockStructuredMatrix *M , double *q , double *z , int verbose, double *w, double *err)
{
  double error = 0;
  int incx = 1, incy = 1;
  double normq = DNRM2(n , q , incx);

  /* w is initialized with q */
  DCOPY(n , q , incx , w , incy);

  /* Computes w += Mz */
  prodSBM(n, n, 1.0, M, z, 1.0, w);

  /* Computes error */
  for (int i = 0 ; i < n ; i++)
    error += abs(z[i] + w[i]) - (z[i] + w[i]);

  *err = error / normq;

  if (verbose > 0) printf("Siconos/Numerics: NCP_compute_error: Error evaluation = %g \n", *err);
}
