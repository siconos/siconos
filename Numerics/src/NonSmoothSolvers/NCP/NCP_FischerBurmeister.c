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
#include "LA.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* Computation of  Fischer-Burmeister function, phi(z,F(z)) = sqrt(z*z + F(z)*F(z)) - z - F(z) */
void phi_FB(int size, double* z, double* F, double* phiVector)
{
  if (z == NULL || F == NULL || phiVector == NULL)
  {
    fprintf(stderr, "FisherBurmeister::phi_FB failed, null input vector(s).\n");
    exit(EXIT_FAILURE);
  }

  int i;
  for (i = 0 ; i < size ; ++i)
  {
    phiVector[i] = sqrt(z[i] * z[i] + F[i] * F[i]) - z[i] - F[i];
  }
}

/* Compute Jacobian of function G */
void jacobianPhi_FB(int size, double* z, double* F, double* jacobianF, double* jacobianPhiMatrix)
{
  if (z == NULL || F == NULL || jacobianF == NULL || jacobianPhiMatrix == NULL)
  {
    fprintf(stderr, "FisherBurmeister::jacobianPhi_FB failed, null input vector(s) or matrices.\n");
    exit(EXIT_FAILURE);
  }

  /* jacobianPhiMatrix is initialized with jacobianF */
  DCOPY(size * size, jacobianF, 1, jacobianPhiMatrix, 1);

  double ri, ai, bi;
  int i;
  for (i = 0; i < size; i++)
  {
    ri = sqrt(z[i] * z[i] +  F[i] * F[i]);
    if (ri > 0.0)
    {
      ai = z[i] / ri - 1.0;
      bi = F[i] / ri - 1.0;
    }
    else
    {
      ai = -1.0;
      bi = -1.0;
    }
    /* jacobianPhiMatrix_ij = delta_ij*ai + bi * jacobianF_ij
       delta_ij being the Kronecker symbol
    */
    DSCAL(size, bi, &jacobianPhiMatrix[i], size);
    jacobianPhiMatrix[i * size + i] += ai;
  }
}

