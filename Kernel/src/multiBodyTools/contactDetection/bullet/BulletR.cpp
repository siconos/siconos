/* Siconos-Kernel, Copyright INRIA 2005-2012.
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

//#define DEBUG_MESSAGES 1
#include <debug.h>

#include "BulletR.hpp"

#include <bullet/BulletCollision/NarrowPhaseCollision/btManifoldPoint.h>
#include <BulletCollision/CollisionDispatch/btCollisionObject.h>

#include <btBulletCollisionCommon.h>

BulletR::BulletR(SP::btManifoldPoint point, SP::btPersistentManifold contactManifold) :
  NewtonEulerFrom3DLocalFrameR(),
  _contactPoints(point),
  _contactManifold(contactManifold)
{
}

void BulletR::computeh(const double time, Interaction& inter)
{
  DEBUG_PRINT("start of computeh\n");

  btCollisionObject* obA =
    static_cast<btCollisionObject*>(_contactManifold->getBody0());
  btCollisionObject* obB =
    static_cast<btCollisionObject*>(_contactManifold->getBody1());

  _contactManifold->refreshContactPoints(obA->getWorldTransform(), obB->getWorldTransform());

  btVector3 posa = _contactPoints->getPositionWorldOnA();
  btVector3 posb = _contactPoints->getPositionWorldOnB();

  (*pc1())(0) = posa[0];
  (*pc1())(1) = posa[1];
  (*pc1())(2) = posa[2];
  (*pc2())(0) = posb[0];
  (*pc2())(1) = posb[1];
  (*pc2())(2) = posb[2];

  btScalar d = sqrt((posb - posa).dot(posb - posa));

  if (fabs(fabs(_contactPoints->getDistance()) - d) <= 1e-5)
  {
    inter.y(0)->setValue(0, _contactPoints->getDistance());

    (*nc())(0) = _contactPoints->m_normalWorldOnB[0];
    (*nc())(1) = _contactPoints->m_normalWorldOnB[1];
    (*nc())(2) = _contactPoints->m_normalWorldOnB[2];
  }
  else
  {
    if (_contactPoints->getDistance() < 0.) d = -d;
    //    printf("WARNING! distance = %g != %g\n",  _contactPoints->getDistance(), d);
    inter.y(0)->setValue(0, d);
    btVector3 n = (posb - posa).normalize();
    (*nc())(0) = n[0];
    (*nc())(1) = n[1];
    (*nc())(2) = n[2];
  }

  DEBUG_PRINTF("distance : %g\n",  inter.y(0)->getValue(0));


  DEBUG_PRINTF("position on A : %g,%g,%g\n", posa[0], posa[1], posa[2]);
  DEBUG_PRINTF("position on B : %g,%g,%g\n", posb[0], posb[1], posb[2]);
  DEBUG_PRINTF("normal on B   : %g,%g,%g\n", (*nc())(0), (*nc())(1), (*nc())(2));

  DEBUG_PRINT("end of computeh\n");


}
