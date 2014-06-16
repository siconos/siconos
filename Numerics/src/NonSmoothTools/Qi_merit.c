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
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

#include "Qi_merit.h"
#include <math.h>
#include <assert.h>
#include <float.h>

void phi_Qi(int n, double* restrict x, double* restrict F, double* restrict Fbox, double* restrict lb, double* restrict ub)
{
  assert(x);
  assert(Fbox);
  assert(lb);
  assert(ub);

  double val;
  double val2;
  double Fi2;

  for (int i = 0; i < n; ++i)
  {
    if (F[i] >= 0.0)
    {
      val = x[i] - lb[i];
      if (x[i] >= lb[i])
      {
        Fi2 = F[i]*F[i];
        val -= sqrt(val*val + Fi2);
        if (x[i] <= ub[i])
        {
          val += F[i]; // done
        }
        else
        {
          val2 = x[i] - ub[i];
          val += sqrt(val2*val2 + Fi2); // done
        }
      }
    }
    else
    {
      val = x[i] - ub[i];
      if (x[i] <= ub[i])
      {
        Fi2 = F[i]*F[i];
        val += sqrt(val*val + Fi2);
        if (x[i] >= lb[i])
        {
          val += F[i]; // done
        }
        else
        {
          val2 = x[i] - lb[i];
          val -= sqrt(val2*val2 + Fi2); // done
        }
      }
    }
    Fbox[i] = val;
  }
}

void Jac_F_Qi(int n, double* restrict x, double* restrict F, double* restrict workV1, double* restrict workV2, double* restrict nabla_F, double* restrict lb, double* restrict ub, double* restrict H)
{
  // function based on the formula given p. 873 in Facchinei--Pang (2003)
  // when the gradient fails to exists, we set a_i = 0 and b_i = 1, as in
  // L. Qi (1999) doi:10.1287/moor.24.2.440

  double normi;
  double normii;
  double diff_l;
  double diff_u;
  double a;
  double b;
  double aa = 1.0-sqrt(2.0)/2.0;
  double bb = 1.0-sqrt(2.0)/2.0;

  // constructing the set beta
  // Introduce a tolerance ? -- xhub
  for (int i = 0; i < n; ++i)
  {
    diff_l = x[i] - lb[i];
    diff_u = x[i] - ub[i];
    if ((fabs(F[i]) < DBL_EPSILON) && ((fabs(diff_l) < DBL_EPSILON) || (fabs(diff_u) < DBL_EPSILON)))
    { // a = 0.0, b = 1.0 ; other choices are possible
      for (int j = 0; j < n; ++j)
      {
        H[j * n + i] = bb*nabla_F[j * n + i];
      }
      H[i * n + i] = aa;
    }
    else // now the rest.... Easy things first
    {
      if (((diff_l <= 0.0) && (F[i] >= 0.0)) || ((diff_u >= 0.0) && (F[i] <= 0.0)))
      {
        a = 1.0;
        b = 0.0;
      }
      // x in the box
      else if ((diff_l >= 0.0) && (diff_u <= 0.0))
      {
        if (F[i] >= 0.0)
        {
          normi = sqrt(diff_l*diff_l + F[i]*F[i]);
          a = 1.0 - diff_l/normi;
          b = 1.0 - F[i]/normi;
        }
        else // F[i] < 0.0
        {
          normi = sqrt(diff_u*diff_u + F[i]*F[i]);
          a = 1.0 + diff_u/normi;
          b = 1.0 + F[i]/normi;
        }
      }
      else // most complex case : x not in the box
      {
        assert(((diff_l < 0.0) && F[i] < 0.0) || ((diff_u > 0.0) && (F[i] > 0.0)));
        normi = sqrt(diff_l*diff_l + F[i]*F[i]);
        normii = sqrt(diff_u*diff_u + F[i]*F[i]);
        a = 1.0 - diff_l/normi + diff_u/normii;
        b = -F[i]/normi + F[i]/normii;
      }

      // now fill H
      for (int j = 0; j < n; ++j)
      {
        H[j * n + i] = b*nabla_F[j * n + i];
      }
      H[i * n + i] += a;
    }
  }
}
