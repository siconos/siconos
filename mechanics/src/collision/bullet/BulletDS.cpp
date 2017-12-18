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

#include "BulletDS.hpp"
#include "BulletDS_impl.hpp"
#include "BulletWeightedShape.hpp"
//#define DEBUG_BULLETDS
#ifdef DEBUG_BULLETDS
#define DEBUG_NOCOLOR
#define DEBUG_STDOUT
#define DEBUG_MESSAGES 1
#endif
#include "debug.h"

// #pragma GCC diagnostic push
// #pragma GCC diagnostic ignored "-Wunknown-pragmas"
// #pragma clang diagnostic push
// #pragma clang diagnostic ignored "-Wunreachable-code"
#include <BulletCollision/CollisionDispatch/btCollisionObject.h>
// #pragma clang diagnostic pop
// #pragma GCC diagnostic push

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
  DEBUG_BEGIN("BulletDS::BulletDS(...)\n");
  DEBUG_EXPR(position->display());
  DEBUG_EXPR(velocity->display());
  DEBUG_EXPR(relative_position->display());
  DEBUG_EXPR(relative_orientation->display());

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
    (*relative_orientation)(0) = 1;
  }

  addCollisionShape(weightedShape->collisionShape(), relative_position,
                      relative_orientation, group);
  DEBUG_END("BulletDS::BulletDS(...)\n");

}

unsigned int BulletDS::numberOfCollisionObjects() const
{
  return _collisionObjects->size();
}

SP::CollisionObjects BulletDS::collisionObjects() const
{
  return _collisionObjects;
}

static
void display_btQuaternion(btQuaternion q)
{
  printf("q0 (w)= %e,q1(x)= %e,q2(y)= %e,q3(z)= %e\n ", q.w(), q.x(), q.y(), q.z());
}

static
void display_btVector3(btVector3 v)
{
  printf("v0 (x)= %e, v1(y)= %e, q2(z)= %e\n ", v.x(), v.y(), v.z());
}




void BulletDS::updateCollisionObjects() const
{
  DEBUG_BEGIN("BulletDS::updateCollisionObjects()\n");
  SiconosVector& q = *_q;

  btQuaternion rbase = btQuaternion(q(4), q(5), q(6), q(3));
  DEBUG_EXPR(q.display(););
  DEBUG_EXPR(display_btQuaternion(rbase););
  for(CollisionObjects::iterator ico = _collisionObjects->begin();
      ico != _collisionObjects->end(); ++ ico)

  {

    SP::btCollisionObject collisionObject = boost::get<0>((*ico).second);
    OffSet offset = boost::get<1>((*ico).second);

    /* with 32bits input ... 1e-7 */
    DEBUG_PRINTF("norm of q - 1.0 = %12.8e\n",fabs(sqrt(pow(q(3), 2) + pow(q(4), 2) +
                                                pow(q(5), 2) +  pow(q(6), 2)) -1. ));

    assert((fabs(sqrt(pow(q(3), 2) + pow(q(4), 2) + pow(q(5), 2) +  pow(q(6), 2)) - 1.) < 1e-7));



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
  DEBUG_END("BulletDS::updateCollisionObjects()\n");

}

void BulletDS::addCollisionObject(SP::btCollisionObject cobj,
                                  SP::SiconosVector pos,
                                  SP::SiconosVector ori,
                                  int group)
{
  std11::array<double, 7> xpos = { {(*pos)(0), (*pos)(1), (*pos)(2), (*ori)(0), (*ori)(1), (*ori)(2), (*ori)(3)}};

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
  DEBUG_BEGIN("BulletDS::setRelativeTransform(...)\n");
  btVector3 rboffset = quatRotate(base_rotation, offset_translation);
  DEBUG_EXPR(display_btQuaternion(base_rotation););
  DEBUG_EXPR(display_btQuaternion(offset_rotation););
  cobj->getWorldTransform().setOrigin(btVector3(
                                        base_translation.x() + rboffset.x(),
                                        base_translation.y() + rboffset.y(),
                                        base_translation.z() + rboffset.z()));

  cobj->getWorldTransform().setRotation(base_rotation * offset_rotation);

  DEBUG_END("BulletDS::setRelativeTransform(...)\n");

}
