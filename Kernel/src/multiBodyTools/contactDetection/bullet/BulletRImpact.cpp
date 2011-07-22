/* Siconos-Kernel, Copyright INRIA 2005-2011.
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

#include "BulletRImpact.hpp"

#include <bullet/BulletCollision/NarrowPhaseCollision/btManifoldPoint.h>
#include <BulletCollision/CollisionDispatch/btCollisionObject.h>

BulletRImpact::BulletRImpact(SP::btManifoldPoint point) : NewtonEulerRImpact(), _contactPoints(point)
{
};

void BulletRImpact::computeh(double)
{
  interaction()->y(0)->setValue(0, _contactPoints->getDistance());
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

};
