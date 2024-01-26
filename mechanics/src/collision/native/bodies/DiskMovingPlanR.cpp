/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

#include "DiskMovingPlanR.hpp"

#include <cmath>

#include "BlockVector.hpp"
#include "PluggedObject.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"

siconos::collision::native::bodies::DiskMovingPlanR::DiskMovingPlanR(
    siconos::plugins::FTime FA, siconos::plugins::FTime FB, siconos::plugins::FTime FC,
    siconos::plugins::FTime FAD, siconos::plugins::FTime FBD, siconos::plugins::FTime FCD,
    double radius)
    : siconos::modeling::LagrangianRheonomousR{}
{
  setComputeAFunction(FA);
  setComputeBFunction(FB);
  setComputeCFunction(FC);
  setComputeADotFunction(FAD);
  setComputeBDotFunction(FBD);
  setComputeCDotFunction(FCD);
  _r = radius;
}

void siconos::collision::native::bodies::DiskMovingPlanR::init(double time)
{
  if (time != _time) {
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

double siconos::collision::native::bodies::DiskMovingPlanR::distance(double x, double y,
                                                                     double rad)
{
  return (fabs(_A * x + _B * y + _C) / _sqrA2pB2 - rad);
}

/* Called compute h, but only the gap function is needed! */
void siconos::collision::native::bodies::DiskMovingPlanR::computeh(
    double time, const siconos::algebra::BlockVector& q, siconos::algebra::BlockVector& z,
    siconos::algebra::SiconosVector& y)
{
  init(time);

  double q_0 = q(0);
  double q_1 = q(1);

  y(0) = distance(q_0, q_1, _r);
}

void siconos::collision::native::bodies::DiskMovingPlanR::computeJachq(
    double time, const siconos::algebra::BlockVector& q, siconos::algebra::BlockVector& z)
{
  init(time);

  double x = q(0);
  double y = q(1);

  double D1 = _A * x + _B * y + _C;
  double signD1 = copysign(1, D1);

  _jachq->setValue(0, 0, _A * signD1 / _sqrA2pB2);
  _jachq->setValue(1, 0, -_B * signD1 / _sqrA2pB2);
  _jachq->setValue(0, 1, _B * signD1 / _sqrA2pB2);
  _jachq->setValue(1, 1, _A * signD1 / _sqrA2pB2);
  _jachq->setValue(0, 2, 0);
  _jachq->setValue(1, 2, -_r);
}

void siconos::collision::native::bodies::DiskMovingPlanR::computehDot(
    double time, const siconos::algebra::BlockVector& q, siconos::algebra::BlockVector& z)
{
  init(time);

  double x = q(0);
  double y = q(1);

  double D1 = _A * x + _B * y + _C;
  double signD1 = copysign(1, D1);
  (*_hDot)(0) = (-_AADot - _BBDot) * fabs(D1) / _cubsqrA2pB2 +
                (_ADot * x + _BDot * y + _CDot) * signD1 / _sqrA2pB2;
}

bool siconos::collision::native::bodies::DiskMovingPlanR::equal(siconos::plugins::FTime pA,
                                                                siconos::plugins::FTime pB,
                                                                siconos::plugins::FTime pC,
                                                                double pr) const
{
  return ((siconos::plugins::FTime)_AFunction->fPtr == pA &&
          (siconos::plugins::FTime)_BFunction->fPtr == pB &&
          (siconos::plugins::FTime)_CFunction->fPtr == pC && _r == pr);
}

void siconos::collision::native::bodies::DiskMovingPlanR::setComputeAFunction(
    siconos::plugins::FTime f)
{
  if (!_AFunction) _AFunction = std::make_shared<siconos::plugins::PluggedObject>();
  _AFunction->setComputeFunction((void*)f);
}

void siconos::collision::native::bodies::DiskMovingPlanR::setComputeBFunction(
    siconos::plugins::FTime f)
{
  if (!_BFunction) _BFunction = std::make_shared<siconos::plugins::PluggedObject>();
  _BFunction->setComputeFunction((void*)f);
}

void siconos::collision::native::bodies::DiskMovingPlanR::setComputeCFunction(
    siconos::plugins::FTime f)
{
  if (!_CFunction) _CFunction = std::make_shared<siconos::plugins::PluggedObject>();
  _CFunction->setComputeFunction((void*)f);
}

void siconos::collision::native::bodies::DiskMovingPlanR::setComputeADotFunction(
    siconos::plugins::FTime f)
{
  if (!_ADotFunction) _ADotFunction = std::make_shared<siconos::plugins::PluggedObject>();
  _ADotFunction->setComputeFunction((void*)f);
}

void siconos::collision::native::bodies::DiskMovingPlanR::setComputeBDotFunction(
    siconos::plugins::FTime f)
{
  if (!_BDotFunction) _BDotFunction = std::make_shared<siconos::plugins::PluggedObject>();
  _BDotFunction->setComputeFunction((void*)f);
}

void siconos::collision::native::bodies::DiskMovingPlanR::setComputeCDotFunction(
    siconos::plugins::FTime f)
{
  if (!_CDotFunction) _CDotFunction = std::make_shared<siconos::plugins::PluggedObject>();
  _CDotFunction->setComputeFunction((void*)f);
}

void siconos::collision::native::bodies::DiskMovingPlanR::computeA(double t)
{
  if (_AFunction->fPtr)
    _A = ((siconos::plugins::FTime)(_AFunction->fPtr))(t);
  else
    _A = 0.;
}

void siconos::collision::native::bodies::DiskMovingPlanR::computeB(double t)
{
  if (_BFunction->fPtr)
    _B = ((siconos::plugins::FTime)(_BFunction->fPtr))(t);
  else
    _B = 0.;
}

void siconos::collision::native::bodies::DiskMovingPlanR::computeC(double t)
{
  if (_CFunction->fPtr)
    _C = ((siconos::plugins::FTime)(_CFunction->fPtr))(t);
  else
    _C = 0.;
}

void siconos::collision::native::bodies::DiskMovingPlanR::computeADot(double t)
{
  if (_ADotFunction->fPtr)
    _ADot = ((siconos::plugins::FTime)(_ADotFunction->fPtr))(t);
  else
    _ADot = 0.;
}

void siconos::collision::native::bodies::DiskMovingPlanR::computeBDot(double t)
{
  if (_BDotFunction->fPtr)
    _BDot = ((siconos::plugins::FTime)(_BDotFunction->fPtr))(t);
  else
    _BDot = 0.;
}

void siconos::collision::native::bodies::DiskMovingPlanR::computeCDot(double t)
{
  if (_CDotFunction->fPtr)
    _CDot = ((siconos::plugins::FTime)(_CDotFunction->fPtr))(t);
  else
    _CDot = 0.;
}
