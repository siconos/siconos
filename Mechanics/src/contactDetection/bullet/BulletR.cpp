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

BulletR::BulletR(unsigned int contact_num, SP::btPersistentManifold contactManifold) :
  NewtonEulerFrom3DLocalFrameR(),
  _contact_num(contact_num),
  _contactManifold(contactManifold)
{
}

void BulletR::computeh(const double time, Interaction& inter)
{
  DEBUG_PRINT("start of computeh\n");

  const btCollisionObject* obA =
    static_cast<const btCollisionObject*>(_contactManifold->getBody0());
  const btCollisionObject* obB =
    static_cast<const btCollisionObject*>(_contactManifold->getBody1());
      
  _contactManifold->refreshContactPoints(obA->getWorldTransform(), 
                                         obB->getWorldTransform());

  unsigned int numContacts = _contactManifold->getNumContacts();

  DEBUG_PRINTF("contact_num : %d\n", _contact_num);

  DEBUG_PRINTF("number of contacts : %d\n", numContacts);

  if (_contact_num > numContacts)
  {
    inter.y(0)->setValue(0, 1e30);
  }
  else
  {

    btManifoldPoint& cpoint = _contactManifold->getContactPoint(_contact_num);

    if ((cpoint.m_normalWorldOnB[0]*cpoint.m_normalWorldOnB[0] +
         cpoint.m_normalWorldOnB[1]*cpoint.m_normalWorldOnB[1] +
         cpoint.m_normalWorldOnB[2]*cpoint.m_normalWorldOnB[2]) > 0.)
      
    {
 
      
      btVector3 posa = cpoint.getPositionWorldOnA();
      btVector3 posb = cpoint.getPositionWorldOnB();
      
      (*pc1())(0) = posa[0];
      (*pc1())(1) = posa[1];
      (*pc1())(2) = posa[2];
      (*pc2())(0) = posb[0];
      (*pc2())(1) = posb[1];
      (*pc2())(2) = posb[2];
      
      inter.y(0)->setValue(0, cpoint.getDistance());
      
      (*nc())(0) = cpoint.m_normalWorldOnB[0];
      (*nc())(1) = cpoint.m_normalWorldOnB[1];
      (*nc())(2) = cpoint.m_normalWorldOnB[2];
      
      
      DEBUG_PRINTF("%d, distance : %g\n",  _contact_num, inter.y(0)->getValue(0));
      DEBUG_PRINTF("%d, position on A : %g,%g,%g\n", _contact_num, posa[0], posa[1], posa[2]);
      DEBUG_PRINTF("%d, position on B : %g,%g,%g\n", _contact_num, posb[0], posb[1], posb[2]);
      DEBUG_PRINTF("%d, normal on B   : %g,%g,%g\n", _contact_num, (*nc())(0), (*nc())(1), (*nc())(2));
      
    }
    else
    {
      btVector3 posa = cpoint.getPositionWorldOnA();
      btVector3 posb = cpoint.getPositionWorldOnB();
      DEBUG_PRINTF("-> %d, distance : %g\n",  _contact_num, cpoint.getDistance());

      DEBUG_PRINTF("-> %d, position on A : %g,%g,%g\n", _contact_num, posa[0], posa[1], posa[2]);
      DEBUG_PRINTF("-> %d, position on B : %g,%g,%g\n", _contact_num, posb[0], posb[1], posb[2]);
      DEBUG_PRINTF("-> %d, normal on B   : %g,%g,%g\n", _contact_num, (*nc())(0), (*nc())(1), (*nc())(2));

      inter.y(0)->setValue(0, fmax(0.1, cpoint.getDistance()));
    }
  }

  DEBUG_PRINT("end of computeh\n");
      
      
}
