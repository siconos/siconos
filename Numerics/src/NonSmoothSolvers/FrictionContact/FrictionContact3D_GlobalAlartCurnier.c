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

#define GAMMA frictionContact3D_gamma

void frictionContact3D_gamma(double* x, int i1, int i2, double *r, int i11, int i12, int i21, int i22, double *p)
{

  assert(x);
  assert(r);
  // p optional

  double invp_, p_, p3;

  if (p)
  {
    p_ = p[0];
  }
  else
  {
    p_ = hypot(x[i1], x[i2]);
  }

  p3 = p_ * p_ * p_;
  invp_ = 1. / p_;
  r[i11] = invp_ - x[i1] * x[i1] / p3;
  r[i12] = - x[i1] * x[i2] / p3;
  r[i21] = r[i12];
  r[i22] = invp_ - x[i2] * x[i2] / p3;

}



void frictionContact3D_GlobalAlartCurnierFunction(
  unsigned int problemSize,
  double *reaction,
  double *velocity,
  double *rho,
  double *mu,
  double *result,
  double *A,
  double *B)
{
  assert(reaction);
  assert(velocity);
  assert(rho);
  assert(mu);
  assert(result);
  assert(problemSize / 3 > 0);
  assert(problemSize % 3 == 0);

  unsigned int ip1, ip2, i, i0, i1, i2, i00, i01, i02, i10, i11, i12, i20, i21, i22;

  for (i = 0, ip1 = 0, ip2 = 0; ip1 < problemSize; ++i)
  {

    i0 = ip1++;
    i1 = ip1++;
    i2 = ip1++;

    if (A && B)
    {
      i00 = ip2++;
      i01 = ip2++;
      i02 = ip2++;
      i10 = ip2++;
      i11 = ip2++;
      i12 = ip2++;
      i20 = ip2++;
      i21 = ip2++;
      i22 = ip2++;

      A[i01] = 0.;
      A[i02] = 0.;
      A[i10] = 0.;
      A[i20] = 0.;

      B[i01] = 0.;
      B[i02] = 0.;

    }


    result[i0] = reaction[i0] - rho[i0] * velocity[i0];
    result[i1] = reaction[i1] - rho[i1] * velocity[i1];
    result[i2] = reaction[i2] - rho[i2] * velocity[i2];

    if (result[i0] > 0)
    {
      result[i0] -= reaction[i0]; // note : this is -PHI p425

      if (A && B)
      {
        // DUnPHI2
        A[i00] = rho[i0];
        B[i00] = 0.;
      }


    }

    else
    {
      result[i0] = -reaction[i0];

      if (A && B)
      {
        A[i00] = 0.;
        B[i00] = 1.;
      }
    }


    double rmu = reaction[i0] * mu[i];
    double p = hypot(result[i1], result[i2]);

    if (p > rmu)
    {

      // outside disk
      if (A && B)
      {
        double mureact, rho1mureact, rho2mureact;
        mureact = mu[i] * reaction[i0];
        rho1mureact = rho[i1] * mureact;
        rho2mureact = rho[i2] * mureact;


        GAMMA(result, i1, i2, A, i11, i12, i21, i22, &p);
        A[i11] *= rho1mureact;
        A[i12] *= rho1mureact;
        A[i21] *= rho2mureact;
        A[i22] *= rho2mureact;

        B[i10] = - mureact * result[i1] / p;
        B[i20] = - mureact * result[i2] / p;

        B[i11] = 1. - A[i11];
        B[i12] = - A[i12];
        B[i21] = - A[i21];
        B[i22] = 1. - A[i22];

      }

      if (rmu <= 0.)
      {
        result[i1] = 0.;
        result[i2] = 0.;
      }
      else
      {
        assert(p > 0.);

        result[i1] *= rmu / p;
        result[i2] *= rmu / p;
      }

    }

    else
    {
      // inside disk
      if (A && B)
      {
        A[i11] = rho[i1];
        A[i12] = 0.;
        A[i21] = 0.;
        A[i22] = rho[i2];

        B[i10] = 0.;
        B[i20] = 0.;

        B[i11] = 0.;
        B[i12] = 0.;
        B[i21] = 0.;
        B[i22] = 0.;

      }

    }

    result[i1] -= reaction[i1];
    result[i2] -= reaction[i2];


    assert(! isnan(result[i0]));
    assert(! isnan(result[i1]));
    assert(! isnan(result[i2]));

    assert(! isnan(A[i00]));
    assert(! isnan(A[i01]));
    assert(! isnan(A[i02]));
    assert(! isnan(A[i10]));
    assert(! isnan(A[i11]));
    assert(! isnan(A[i12]));
    assert(! isnan(A[i20]));
    assert(! isnan(A[i21]));
    assert(! isnan(A[i22]));


    assert(! isnan(B[i00]));
    assert(! isnan(B[i01]));
    assert(! isnan(B[i02]));
    assert(! isnan(B[i10]));
    assert(! isnan(B[i11]));
    assert(! isnan(B[i12]));
    assert(! isnan(B[i20]));
    assert(! isnan(B[i21]));
    assert(! isnan(B[i22]));


  }



}


