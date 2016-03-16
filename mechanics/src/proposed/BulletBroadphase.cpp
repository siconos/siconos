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

/*! \file BulletBroadphase.cpp
  \brief Implementation of a Bullet-based broadphase algorithm.
*/

#include "BulletBroadphase.hpp"

#include <stdio.h>

#include <BulletCollision/CollisionShapes/btSphereShape.h>
#include <BulletCollision/CollisionShapes/btBoxShape.h>

DEFINE_SPTR(btSphereShape);
DEFINE_SPTR(btBoxShape);

class BulletShapeHandler : public SiconosShapeHandler
{
protected:
  std::vector<SP::SiconosSphere> dirtySpheres;
  std::vector<SP::SiconosBox> dirtyBoxes;

  friend class BulletBroadphase;

public:
  BulletShapeHandler(BulletBroadphase *_impl) {}

  virtual void onChanged(SP::SiconosSphere sphere);
  virtual void onChanged(SP::SiconosBox box);
};

void BulletShapeHandler::onChanged(SP::SiconosSphere sphere)
{
  dirtySpheres.push_back(sphere);
  printf("pushed sphere: %p(%ld)\n",
         &*sphere,sphere.use_count());
}

void BulletShapeHandler::onChanged(SP::SiconosBox box)
{
    dirtyBoxes.push_back(box);
}

BulletBroadphase::BulletBroadphase() {
  handler.reset(new BulletShapeHandler(this));
  _collisionConfiguration.reset(new btDefaultCollisionConfiguration());
  _dispatcher.reset(new btCollisionDispatcher(&*_collisionConfiguration));
  _broadphase.reset(new btDbvtBroadphase());
  _collisionWorld.reset(new btCollisionWorld(&*_dispatcher, &*_broadphase,
                                             &*_collisionConfiguration));
}

BulletBroadphase::~BulletBroadphase() {
  // must be the first de-allocated, otherwise segfault
  _collisionWorld.reset();
}

void BulletBroadphase::visit(SP::SiconosSphere sphere)
{
  printf("contactor: %p, ", currentContactor);
  printf("sphere: %p(%ld)\n",
         &*sphere,sphere.use_count());

  // create corresponding Bullet object and shape
  btCollisionObject *btobject = new btCollisionObject();
  btCollisionShape *btshape = new btSphereShape(sphere->radius());
  btobject->setCollisionShape(btshape);

  // put it in the world
  _collisionWorld->addCollisionObject(btobject);

  // install handler for updates
  sphere->setHandler(handler);
}

void BulletBroadphase::visit(SP::Contactor contactor)
{
  currentContactor = contactor;

  std::vector<SP::SiconosShape>::const_iterator it;
  for (it=contactor->shapes().begin();
       it!=contactor->shapes().end();
       it++)
  {
    (*it)->acceptSP(shared_from_this());
  }
}

void BulletBroadphase::buildGraph(SP::Contactor contactor)
{
  contactor->acceptSP(shared_from_this());
}

void BulletBroadphase::updateGraph()
{
  printf("dirty spheres:\n");
  std::vector<SP::SiconosSphere>::iterator it;
  for (it=handler->dirtySpheres.begin();
       it!=handler->dirtySpheres.end(); it++)
  {
    printf("  %p\n", &**it);
  }
  handler->dirtySpheres.clear();

  printf("dirty boxes:\n");
  std::vector<SP::SiconosBox>::iterator it2;
  for (it2=handler->dirtyBoxes.begin();
       it2!=handler->dirtyBoxes.end(); it2++)
  {
    printf("  %p\n", &**it2);
  }
  handler->dirtyBoxes.clear();
}

void BulletBroadphase::performBroadphase()
{
  _collisionWorld->performDiscreteCollisionDetection();
}
