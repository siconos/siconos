/* Siconos-Numerics, Copyright INRIA 2005-2014
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

#include "SiconosBlas.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "assert.h"

void F_min(int n1, int n2, double* restrict z, double* restrict F, double* restrict Fmin)
{
  assert(z != NULL);
  assert(F != NULL);
  assert(Fmin != NULL);

  for (int i = 0; i < n1; ++i)
  {
    Fmin[i] = F[i];
  }
  for (int i = n1 ; i < n1+n2 ; ++i)
  {
    Fmin[i] = z[i] <= F[i] ? z[i] : F[i];
  }
}

void Jac_F_min(int n1, int n2, double* restrict z, double* restrict F, double* restrict nabla_F, double* restrict H)
{

  int n = n1 + n2;
  if (n1 > 0)
  {
    //printf("Jac_F_min: mixed case has to validated -- xhub\n");
    cblas_dcopy(n*n, nabla_F, 1, H, 1);
  }
  // Facchinei--Pang p. 660 and 661
  // i \in alpha if F_i > z_i
  // i \in beta if F_i == z_i
  // i \in gamma if F_i < z_i
  // We made the choice here to have the trivial case z_i + d_i = 0 for alpha U
  // beta. See Facchinei--Pang for more details
  //
  // TODO implement sparse version where the trivial case can be fully exploited
  for (int i = n1; i < n; ++i)
  {
    if (z[i] <= F[i]) // i in beta U alpha
    {
      for (int j = 0; j < n; ++j)
      {
        H[j * n + i] = 0.0;
      }
      H[i * n + i] = 1.0;

    }
    else // i in gamma
    {
      for (int j = 0; j < n; ++j)
      {
        H[j * n + i] = nabla_F[j * n + i];
      }
    }
  }
}
