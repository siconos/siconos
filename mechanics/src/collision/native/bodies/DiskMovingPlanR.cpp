/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

#include <cmath>
#include "DiskMovingPlanR.hpp"
#include <BlockVector.hpp>
#include "SimpleMatrix.hpp"

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
void DiskMovingPlanR::computeh(double time, SiconosVector& q, SiconosVector& z, SiconosVector& y)
{
  init(time);

  double q_0 = q(0);
  double q_1 = q(1);

  y(0) = distance(q_0, q_1, _r);

}

void DiskMovingPlanR::computeJachq(double time, SiconosVector& q, SiconosVector& z)
{
  init(time);

  SimpleMatrix *g = (SimpleMatrix *) _jachq.get();

  double x = q(0);
  double y = q(1);

  double D1 = _A * x + _B * y + _C;
  double signD1 = copysign(1, D1);

  (*g)(0, 0) = _A * signD1 / _sqrA2pB2;
  (*g)(1, 0) = -_B * signD1 / _sqrA2pB2;
  (*g)(0, 1) = _B * signD1 / _sqrA2pB2;
  (*g)(1, 1) = _A * signD1 / _sqrA2pB2;
  (*g)(0, 2) = 0;
  (*g)(1, 2) = -_r;
}

void DiskMovingPlanR::computehDot(double time, SiconosVector& q, SiconosVector& z)
{
  init(time);

  double x = q(0);
  double y = q(1);

  double D1 = _A * x + _B * y + _C;
  double signD1 = copysign(1, D1);
  (*_hDot)(0) = (-_AADot - _BBDot) * fabs(D1) / _cubsqrA2pB2 + (_ADot * x + _BDot * y + _CDot) * signD1 / _sqrA2pB2;
}

bool DiskMovingPlanR::equal(FTime pA, FTime pB, FTime pC, double pr) const
{
  return ((FTime)_AFunction->fPtr == pA && (FTime)_BFunction->fPtr == pB &&
          (FTime)_CFunction->fPtr == pC && _r == pr);
}

