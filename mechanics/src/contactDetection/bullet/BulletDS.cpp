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

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunreachable-code"
#include <BulletCollision/CollisionDispatch/btCollisionObject.h>
#pragma clang diagnostic pop

BulletDS::BulletDS(SP::BulletWeightedShape weightedShape,
                   SP::SiconosVector position,
                   SP::SiconosVector velocity,
                   SP::SiconosVector relative_position,
                   SP::SiconosVector relative_orientation,
                   int group) :
  NewtonEulerDS(position, velocity, weightedShape->mass(),
                weightedShape->inertia()),
  _weightedShape(weightedShape),
  _collisionObjects(new CollisionObjects())
{
  SiconosVector& q = *_q;


  /* with 32bits input ... 1e-7 */
  if (fabs(sqrt(pow(q(3), 2) + pow(q(4), 2) +
                pow(q(5), 2) + pow(q(6), 2)) - 1.) >= 1e-7)
  {
    RuntimeException::selfThrow(
      "BulletDS: quaternion in position parameter is not a unit quaternion "
    );
  }

  /* initialisation is done with the weighted shape as the only one
   * collision object */
  if (! relative_position)
  {
    relative_position.reset(new SiconosVector(3));
    relative_position->zero();
  }
  if (! relative_orientation)
  {
    relative_orientation.reset(new SiconosVector(4));
    relative_orientation->zero();
    (*relative_orientation)(1) = 1;
  }

  addCollisionShape(weightedShape->collisionShape(), relative_position,
                      relative_orientation, group);

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

  SiconosVector& q = *_q;

  btQuaternion rbase = btQuaternion(q(4), q(5), q(6), q(3));


  for(CollisionObjects::iterator ico = _collisionObjects->begin();
      ico != _collisionObjects->end(); ++ ico)

  {

    SP::btCollisionObject collisionObject = boost::get<0>((*ico).second);
    OffSet offset = boost::get<1>((*ico).second);

    DEBUG_EXPR(q.display());

    /* with 32bits input ... 1e-7 */
    assert(fabs(sqrt(pow(q(3), 2) + pow(q(4), 2) +
                   pow(q(5), 2) +  pow(q(6), 2)) - 1.) < 1e-7);



/*    btVector3 boffset = btVector3(offset[0], offset[1], offset[2]);

    btVector3 rboffset = quatRotate(rbase, boffset);

    collisionObject->getWorldTransform().setOrigin(btVector3(q(0)+rboffset.x(),
                                                             q(1)+rboffset.y(),
                                                             q(2)+rboffset.z()));
    collisionObject->getWorldTransform().getBasis().
      setRotation(rbase * btQuaternion(offset[4], offset[5],
      offset[6], offset[3]));*/

    this->setRelativeTransform(collisionObject,
                               btVector3(q(0), q(1), q(2)),
                               rbase,
                               btVector3(offset[0], offset[1], offset[2]),
                               btQuaternion(offset[4], offset[5],
                                            offset[6], offset[3]));


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
  std11::array<double, 7> xpos = { (*pos)(0), (*pos)(1), (*pos)(2),
                                   (*ori)(0), (*ori)(1), (*ori)(2), (*ori)(3)};

  (*_collisionObjects)[&*cobj] =  boost::tuple<SP::btCollisionObject,
                                               OffSet , int>
    (cobj, xpos, group);

  updateCollisionObjects();
}

void BulletDS::addCollisionShape(SP::btCollisionShape shape,
                                 SP::SiconosVector pos,
                                 SP::SiconosVector ori,
                                 int group)
{
  SP::btCollisionObject collisionObject(new btCollisionObject());
  collisionObject->setUserPointer(this);
  collisionObject->setCollisionFlags(collisionObject->getCollisionFlags()|
                                     btCollisionObject::CF_KINEMATIC_OBJECT);
  collisionObject->setCollisionShape(&*shape);

  addCollisionObject(collisionObject, pos, ori, group);
}



void BulletDS::setRelativeTransform(SP::btCollisionObject cobj,
                                    const btVector3& base_translation,
                                    const btQuaternion& base_rotation,
                                    const btVector3& offset_translation,
                                    const btQuaternion& offset_rotation)
{
  btVector3 rboffset = quatRotate(base_rotation, offset_translation);

  cobj->getWorldTransform().setOrigin(btVector3(
                                        base_translation.x() + rboffset.x(),
                                        base_translation.y() + rboffset.y(),
                                        base_translation.z() + rboffset.z()));

  cobj->getWorldTransform().setRotation(
    base_rotation * offset_rotation);
}
