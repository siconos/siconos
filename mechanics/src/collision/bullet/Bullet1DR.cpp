/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

#include "Bullet1DR.hpp"

#include <Interaction.hpp>
#include <BulletCollision/NarrowPhaseCollision/btManifoldPoint.h>
#include <BulletCollision/CollisionDispatch/btCollisionObject.h>

Bullet1DR::Bullet1DR(SP::btManifoldPoint point) : NewtonEuler1DR(), _contactPoints(point)
{
}

void Bullet1DR::computeh(double time, BlockVector& q0, SiconosVector& y)
{
  y.setValue(0, _contactPoints->getDistance());
  btVector3 posa = _contactPoints->getPositionWorldOnA();
  btVector3 posb = _contactPoints->getPositionWorldOnB();
  (*pc1())(0) = posa[0];
  (*pc1())(1) = posa[1];
  (*pc1())(2) = posa[2];
  (*pc2())(0) = posb[0];
  (*pc2())(1) = posb[1];
  (*pc2())(2) = posb[2];

  (*nc())(0) = _contactPoints->m_normalWorldOnB[0];
  (*nc())(1) = _contactPoints->m_normalWorldOnB[1];
  (*nc())(2) = _contactPoints->m_normalWorldOnB[2];;

}
