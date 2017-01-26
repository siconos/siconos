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

#include "BulletWeightedShape.hpp"

#include <BulletCollision/CollisionShapes/btCollisionShape.h>

#include <SimpleMatrix.hpp>

#include <CxxStd.hpp>

BulletWeightedShape::BulletWeightedShape(SP::btCollisionShape shape, const double& mass) :
  _mass(mass),
  _shape(shape)
{
  btVector3 localinertia;
  _shape->calculateLocalInertia(_mass, localinertia);

  assert( ! ((localinertia.x() == 0.0
              && localinertia.y() == 0.0
              && localinertia.z() == 0.0)
             || isinf(localinertia.x())
             || isinf(localinertia.y())
             || isinf(localinertia.z()))
          && "calculateLocalInertia() returned garbage" );

  _inertia.reset(new SimpleMatrix(3, 3));

  _inertia->zero();

  (*_inertia)(0, 0) = localinertia[0];
  (*_inertia)(1, 1) = localinertia[1];
  (*_inertia)(2, 2) = localinertia[2];

}


void BulletWeightedShape::setInertia(double ix, double iy, double iz)
{
  _inertia->zero();

  (*_inertia)(0, 0) = ix;
  (*_inertia)(1, 1) = iy;
  (*_inertia)(2, 2) = iz;
}
