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

#include "LA.h"
#include "FrictionContact3D_Solvers.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#define ACFun frictionContact3D_GlobalAlartCurnierFunction
#define DUnACFunT1 frictionContact3D_GlobalAlartCurnierChooseDUnACFunT1
#define DRnACFunT1 frictionContact3D_GlobalAlartCurnierChooseDRnACFunT1
#define DRnACFunT2 frictionContact3D_GlobalAlartCurnierChooseDRnACFunT2
#define DUtACFunT1 frictionContact3D_GlobalAlartCurnierChooseDUtACFunT1
#define DUtACFunT2 frictionContact3D_GlobalAlartCurnierChooseDUtACFunT2

void frictionContact3D_GlobalAlartCurnierFunction(
  unsigned int dim,
  double *reaction,
  double *velocity,
  double *rho,
  double *mu,
  double *result)
{
  assert(reaction);
  assert(velocity);
  assert(rho);
  assert(mu);
  assert(result);
  assert(dim / 3 > 0);
  assert(dim % 3 == 0);

  unsigned int i, iN, iT1, iT2;

  for (i = 0, iN = 0, iT1 = 1, iT2 = 2; iN < dim; iN += 3, iT1 += 3, iT2 += 3)
  {
    result[iN] = reaction[iN] - rho[iN] * velocity[iN];
    result[iT1] = reaction[iT1] - rho[iT1] * velocity[iT1];
    result[iT2] = reaction[iT2] - rho[iT2] * velocity[iT2];

    if (result[iN] > 0)
      result[iN] -= reaction[iN];
    else
      result[iN] = -reaction[iN];

    double rmu = reaction[iN] * mu[i];

    if (rmu == 0.)
    {
      result[iT1] = 0.;
      result[iT2] = 0.;
    }
    else
    {
      double p = hypot(result[iT1], result[iT2]);
      if ((p > rmu) && (p > 0.))
      {
        result[iT1] *= rmu / p;
        result[iT2] *= rmu / p;
      }
    }

    result[iT1] -= reaction[iT1];
    result[iT2] -= reaction[iT2];

  }

}

double frictionContact3D_GlobalAlartCurnierChooseDUnACFunT1(
  unsigned int frictionContact3D_problemSize,
  double *reaction,
  double *velocity,
  double *rho)
{
  return 0.;
};

double frictionContact3D_GlobalAlartCurnierChooseDRnACFunT1(
  unsigned int frictionContact3D_problemSize,
  double *reaction,
  double *velocity,
  double *rhoN)
{
  return 0.;
};

double frictionContact3D_GlobalAlartCurnierChooseDRnACFunT2(
  unsigned int frictionContact3D_problemSize,
  double *reaction,
  double *velocity,
  double *rho)
{
  return 0.;
};

void frictionContact3D_GlobalAlartCurnierChooseDUtACFunT2(
  unsigned int frictionContact3D_problemSize,
  double *reaction,
  double *velocity,
  double *rho,
  double *result2x2) {};

void frictionContact3D_GlobalAlartCurnierChooseDRtACFunT2(
  unsigned int frictionContact3D_problemSize,
  double *reaction,
  double *velocity,
  double *rho,
  double *result2x2) {};



#undef ACFun      // frictionContact3D_GlobalAlartCurnierFunction
#undef DUnACFunT1 // frictionContact3D_GlobalAlartCurnierChooseDUnACFunT1
#undef DRnACFunT1 // frictionContact3D_GlobalAlartCurnierChooseDRnACFunT1
#undef DRnACFunT2 // frictionContact3D_GlobalAlartCurnierChooseDRnACFunT2
#undef DUtACFunT1 // frictionContact3D_GlobalAlartCurnierChooseDUtACFunT1
#undef DUtACFunT2 // frictionContact3D_GlobalAlartCurnierChooseDUtACFunT2a
