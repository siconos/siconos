/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#include "SphereLDS.hpp"

#include <cmath>

#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"

constexpr double TWOPI = 2.0 * M_PI;

static void normalize(siconos::algebra::SiconosVector& q, unsigned int i)
{
  q.setValue(i, fmod(q.getValue(i), TWOPI));

  assert(fabs(q.getValue(i)) - std::numeric_limits<double>::epsilon() >= 0.);
  assert(fabs(q.getValue(i)) < TWOPI);
}

siconos::collision::native::bodies::SphereLDS::SphereLDS(
    double r, double m, std::shared_ptr<siconos::algebra::SiconosVector> qinit,
    std::shared_ptr<siconos::algebra::SiconosVector> vinit)
  : siconos::modeling::LagrangianDS{qinit, vinit}, radius{r}, massValue{m}
{
  normalize(*q(), 3);
  normalize(*q(), 4);
  normalize(*q(), 5);
  _ndof = 6;

  assert(qinit->size() == _ndof);
  assert(vinit->size() == _ndof);

  _mass = std::make_shared<siconos::algebra::SimpleMatrix>(_ndof, _ndof);
  _mass->zero();
  I = massValue * radius * radius * 2. / 5.;
  (*_mass)(0, 0) = (*_mass)(1, 1) = (*_mass)(2, 2) = massValue;
  ;
  (*_mass)(3, 3) = (*_mass)(4, 4) = (*_mass)(5, 5) = I;

  computeMass();

  _jacobianFGyrq = std::make_shared<siconos::algebra::SimpleMatrix>(_ndof, _ndof);
  _jacobianFGyrqDot = std::make_shared<siconos::algebra::SimpleMatrix>(_ndof, _ndof);

  _fGyr = std::make_shared<siconos::algebra::SiconosVector>(_ndof);
  _fGyr->zero();

  computeFGyr();
}

void siconos::collision::native::bodies::SphereLDS::computeMass()
{
  normalize(*q(), 3);
  normalize(*q(), 4);
  normalize(*q(), 5);

  // // SS: Forcing modification of qold, is this necessary?
  // if (qMemory().nbVectorsInMemory() >= 1)
  // {
  //   SiconosVector& qold = qMemory().getsiconos::algebra::SiconosVector(0);
  //   normalize(qold, 3);
  //   normalize(qold, 4);
  //   normalize(qold, 5);
  // }

  double theta = q()->getValue(3);

  assert(fabs(theta) - std::numeric_limits<double>::epsilon() >= 0.);
  // assert (fabs(theta) - TWOPI < 0.);

  (*_mass)(4, 5) = (*_mass)(5, 4) = I * cos(theta);
}

void siconos::collision::native::bodies::SphereLDS::computeFGyr()
{
  siconos::collision::native::bodies::SphereLDS::computeFGyr(q(), velocity());
}

void siconos::collision::native::bodies::SphereLDS::computeFGyr(
    std::shared_ptr<siconos::algebra::SiconosVector> q,
    std::shared_ptr<siconos::algebra::SiconosVector> v)
{
  assert(q->size() == 6);
  assert(v->size() == 6);

  //  normalize(q,3);
  // normalize(q,4);
  // normalize(q,5);

  double theta = q->getValue(3);

  double thetadot = v->getValue(3);
  double phidot = v->getValue(4);
  double psidot = v->getValue(5);

  double sintheta = sin(theta);

  (*_fGyr)(0) = (*_fGyr)(1) = (*_fGyr)(2) = 0;

  (*_fGyr)(3) = I * psidot * phidot * sintheta;
  (*_fGyr)(4) = -I * psidot * thetadot * sintheta;
  (*_fGyr)(5) = -I * phidot * thetadot * sintheta;
}

void siconos::collision::native::bodies::SphereLDS::computeJacobianFGyrq()
{
  siconos::collision::native::bodies::SphereLDS::computeJacobianFGyrq(q(), velocity());
}
void siconos::collision::native::bodies::SphereLDS::computeJacobianFGyrqDot()
{
  siconos::collision::native::bodies::SphereLDS::computeJacobianFGyrqDot(q(), velocity());
}

void siconos::collision::native::bodies::SphereLDS::computeJacobianFGyrq(
    std::shared_ptr<siconos::algebra::SiconosVector> q,
    std::shared_ptr<siconos::algebra::SiconosVector> v)
{
  double theta = q->getValue(3);

  double thetadot = v->getValue(3);
  double phidot = v->getValue(4);
  double psidot = v->getValue(5);

  double costheta = cos(theta);

  _jacobianFGyrq->zero();

  (*_jacobianFGyrq)(3, 3) = -I * psidot * phidot * costheta;
  (*_jacobianFGyrq)(4, 3) = I * psidot * thetadot * costheta;
  (*_jacobianFGyrq)(5, 3) = I * psidot * thetadot * costheta;
}
void siconos::collision::native::bodies::SphereLDS::computeJacobianFGyrqDot(
    std::shared_ptr<siconos::algebra::SiconosVector> q,
    std::shared_ptr<siconos::algebra::SiconosVector> v)
{
  double theta = q->getValue(3);

  double thetadot = v->getValue(3);
  double phidot = v->getValue(4);
  double psidot = v->getValue(5);

  double sintheta = sin(theta);

  _jacobianFGyrqDot->zero();

  (*_jacobianFGyrqDot)(3, 3) = 0;
  (*_jacobianFGyrqDot)(3, 4) = I * psidot * sintheta;
  (*_jacobianFGyrqDot)(3, 5) = I * phidot * sintheta;

  (*_jacobianFGyrqDot)(4, 3) = -I * psidot * sintheta;
  (*_jacobianFGyrqDot)(4, 4) = 0;
  (*_jacobianFGyrqDot)(4, 5) = -I * thetadot * sintheta;

  (*_jacobianFGyrqDot)(5, 3) = -I * phidot * sintheta;
  (*_jacobianFGyrqDot)(5, 4) = -I * thetadot * sintheta;
  (*_jacobianFGyrqDot)(5, 5) = 0;
}

double siconos::collision::native::bodies::SphereLDS::getQ(unsigned int pos)
{
  assert(pos < _ndof);
  return (*_q[0])(pos);
};

double siconos::collision::native::bodies::SphereLDS::getVelocity(unsigned int pos)
{
  assert(pos < _ndof);
  return (*_q[1])(pos);
};
