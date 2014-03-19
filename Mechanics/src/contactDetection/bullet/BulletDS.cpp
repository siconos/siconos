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

#include "BulletDS.hpp"
#include "BulletDS_impl.hpp"
#include "BulletWeightedShape.hpp"

#ifdef DEBUG_BULLETDS
#define DEBUG_MESSAGES 1
#endif
#include "debug.h"

#include <BulletCollision/CollisionDispatch/btCollisionObject.h>

BulletDS::BulletDS(SP::BulletWeightedShape weightedShape,
                   SP::SiconosVector position,
                   SP::SiconosVector velocity) :
  NewtonEulerDS(position, velocity, weightedShape->mass(),
                weightedShape->inertiaMatrix()),
  _weightedShape(weightedShape),
  _collisionObjects(new CollisionObjects())
{
  SiconosVector& q = *_q;

  if (fabs(sqrt(pow(q(3), 2) + pow(q(4), 2) +
                pow(q(5), 2) + pow(q(6), 2)) - 1.) >= 1e-10)
  {
    RuntimeException::selfThrow(
      "BulletDS: quaternion in position parameter is not a unit quaternion "
    );
  }

  /* initialisation is done with the weighted shape as the only one
   * collision object */
  SP::btCollisionObject collisionObject(new btCollisionObject());

  collisionObject->setUserPointer(this);
  collisionObject->setCollisionFlags(collisionObject->getCollisionFlags()|
                                     btCollisionObject::CF_KINEMATIC_OBJECT);

  collisionObject->setCollisionShape(&*(weightedShape->collisionShape()));

  boost::array<double, 7> centerOfMass = { 0,0,0,1,0,0,0 };

  (*_collisionObjects)[&*collisionObject] =
    boost::tuple<SP::btCollisionObject, OffSet , int>
    (collisionObject,centerOfMass,0);

  updateCollisionObjects();
}

unsigned int BulletDS::numberOfCollisionObjects() const
{
  return _collisionObjects->size();
}

SP::CollisionObjects BulletDS::collisionObjects() const
{
  return _collisionObjects;
}


void BulletDS::updateCollisionObjects() const
{

  for(CollisionObjects::iterator ico = _collisionObjects->begin();
      ico != _collisionObjects->end(); ++ ico)

  {

    SP::btCollisionObject collisionObject = boost::get<0>((*ico).second);
    OffSet offset = boost::get<1>((*ico).second);

    SiconosVector& q = *_q;

    DEBUG_EXPR(q.display());

    assert(fabs(sqrt(pow(q(3), 2) + pow(q(4), 2) +
                   pow(q(5), 2) +  pow(q(6), 2)) - 1.) < 1e-10);

    collisionObject->getWorldTransform().setOrigin(btVector3(q(0)+offset[0],
                                                             q(1)+offset[1],
                                                             q(2)+offset[2]));
    collisionObject->getWorldTransform().getBasis().
      setRotation(btQuaternion(offset[4], offset[5],
                               offset[6], offset[3]) *
                  btQuaternion(q(4), q(5),
                               q(6), q(3)));

    /* is this needed ? */
    collisionObject->setActivationState(ACTIVE_TAG);
    collisionObject->activate();

  //_collisionObject->setContactProcessingThreshold(BT_LARGE_FLOAT);

  }
}

void BulletDS::addCollisionObject(SP::btCollisionObject cobj,
                                  SP::SiconosVector pos,
                                  SP::SiconosVector ori,
                                  int group)
{
  boost::array<double, 7> xpos = { (*pos)(0), (*pos)(1), (*pos)(2),
                                   (*ori)(0), (*ori)(1), (*ori)(2), (*ori)(3)};

  (*_collisionObjects)[&*cobj] =  boost::tuple<SP::btCollisionObject,
                                               OffSet , int>
    (cobj, xpos, group);
}
