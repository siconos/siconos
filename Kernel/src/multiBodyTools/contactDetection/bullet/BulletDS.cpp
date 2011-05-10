/* Siconos-Kernel, Copyright INRIA 2005-2010.
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

#include "BulletDS.hpp"
#include "BulletWeightedShape.hpp"

#ifdef DEBUG_BULLETDS
#define DEBUG_MESSAGES 1
#endif
#include <debug.h>

#include <BulletCollision/CollisionDispatch/btCollisionObject.h>




BulletDS::BulletDS(SP::BulletWeightedShape weightedShape,
                   SP::SiconosVector position,
                   SP::SiconosVector velocity) :
  NewtonEulerDS(position, velocity, weightedShape->mass(), weightedShape->inertiaMatrix()),
  _weightedShape(weightedShape)
{

  _collisionObject.reset(new btCollisionObject());

  _collisionObject->setUserPointer(this);
  _collisionObject->setCollisionFlags(_collisionObject->getCollisionFlags() | btCollisionObject::CF_KINEMATIC_OBJECT);

  _collisionObject->setCollisionShape(&*(weightedShape->collisionShape()));

  updateCollisionObject();
}

void BulletDS::updateCollisionObject() const
{

  DEBUG_PRINT("updateCollisionObject()");

  SimpleVector& q = *_q;

  DEBUG_EXPR(q.display());


  btMatrix3x3 basis;
  basis.setIdentity();
  _collisionObject->getWorldTransform().setBasis(basis);

  assert(fabs(sqrt(pow(q(3), 2) + pow(q(4), 2) +  pow(q(5), 2) +  pow(q(6), 2)) - 1.) < 1e-10);
  _collisionObject->getWorldTransform().getOrigin().setX(q(0));
  _collisionObject->getWorldTransform().getOrigin().setY(q(1));
  _collisionObject->getWorldTransform().getOrigin().setZ(q(2));

  _collisionObject->getWorldTransform().getBasis().setRotation(btQuaternion(q(4), q(5),
      q(6), q(3)));


  _collisionObject->setActivationState(ACTIVE_TAG);

}
