/* Siconos-Numerics, Copyright INRIA 2005-2012.
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
#include <math.h>
#include "projectionOnCone.h"

unsigned projectionOnCone(double* r, double  mu)
{
  double normT = hypot(r[1], r[2]);
  if (mu * normT <= - r[0])
  {
    r[0] = 0.0;
    r[1] = 0.0;
    r[2] = 0.0;
    return PROJCONE_DUAL;
  }
  else if (normT <= mu * r[0])
  {
    return PROJCONE_INSIDE;
  }
  else
  {
    double mu2 = mu * mu;
    r[0] = (mu * normT + r[0]) / (mu2 + 1.0);
    r[1] = mu * r[0] * r[1] / normT;
    r[2] = mu * r[0] * r[2] / normT;
    return PROJCONE_BOUNDARY;
  }
}

void projectionOnSecondOrderCone(double* r, double  mu, int size)
{
  if (size ==3)
  {
    projectionOnCone(r, mu);
  }
  else
  {
    fprintf(stderr, "Numerics, projectionOnSecondOrderCone f not yet implementes for size != 3 \n");
    exit(EXIT_FAILURE);
  }

}
