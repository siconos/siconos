/* Siconos-Example version 3.0.0, Copyright INRIA 2005-2010.
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
 * Foundation, Inc., 51 Franklin St, Fifth FLOOR, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 *
 */

#include <math.h>
#include "DiskMovingPlanR.hpp"

DiskMovingPlanR::DiskMovingPlanR(FTime FA, FTime FB, FTime FC,
                                 FTime FAD, FTime FBD, FTime FCD,
                                 double radius) : LagrangianRheonomousR()
{
  setComputeAFunction(FA);
  setComputeBFunction(FB);
  setComputeCFunction(FC);
  setComputeADotFunction(FAD);
  setComputeBDotFunction(FBD);
  setComputeCDotFunction(FCD);
  _r = radius;
}


void DiskMovingPlanR::init(double time)
{
  if (time != _time)
  {
    _time = time;
    computeA(time);
    computeB(time);
    computeC(time);
    computeADot(time);
    computeBDot(time);
    computeCDot(time);

    _sqrA2pB2 = hypot(_A, _B);
    _AADot = _A * _ADot;
    _BBDot = _B * _BDot;
    _cubsqrA2pB2 = _sqrA2pB2 * _sqrA2pB2 * _sqrA2pB2;
  }
}

double DiskMovingPlanR::distance(double x, double y, double rad)
{
  return (fabs(_A * x + _B * y + _C) / _sqrA2pB2 - rad);
}

/* Called compute h, but only the gap function is needed! */
void DiskMovingPlanR::computeh(double time)
{
  init(time);

  SiconosVector *y = interaction()->y(0).get();
  double q_0 = (*data[q0])(0);
  double q_1 = (*data[q0])(1);

  (*y)(0) = distance(q_0, q_1, _r);

}

void DiskMovingPlanR::computeJachq(double time)
{
  init(time);

  SimpleMatrix *g = (SimpleMatrix *) _jachq.get();

  double x = (*data[q0])(0);
  double y = (*data[q0])(1);

  double D1 = _A * x + _B * y + _C;
  double signD1 = copysign(1, D1);

  (*g)(0, 0) = _A * signD1 / _sqrA2pB2;
  (*g)(1, 0) = -_B * signD1 / _sqrA2pB2;
  (*g)(0, 1) = _B * signD1 / _sqrA2pB2;
  (*g)(1, 1) = _A * signD1 / _sqrA2pB2;
  (*g)(0, 2) = 0;
  (*g)(1, 2) = -_r;
}

void DiskMovingPlanR::computehDot(double time)
{
  init(time);

  double x = (*data[q0])(0);
  double y = (*data[q0])(1);

  double D1 = _A * x + _B * y + _C;
  double signD1 = copysign(1, D1);
  (*_hDot)(0) = (-_AADot - _BBDot) * fabs(D1) / _cubsqrA2pB2 + (_ADot * x + _BDot * y + _CDot) * signD1 / _sqrA2pB2;
}

bool DiskMovingPlanR::equal(FTime pA, FTime pB, FTime pC, double pr) const
{
  return (_AFunction->fPtr == pA && _BFunction->fPtr == pB &&
          _CFunction->fPtr == pC && _r == pr);
}

