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

/*! \file SiconosBroadphase.hpp
  \brief Definition of an abstract broadphase algorithm.
*/

#ifndef BulletBroadphase_h
#define BulletBroadphase_h

#include <SiconosBroadphase.hpp>
#include <MechanicsFwd.hpp>

#include <BulletCollision/CollisionDispatch/btCollisionWorld.h>
#include <BulletCollision/CollisionShapes/btSphereShape.h>
#include <BulletCollision/Gimpact/btGImpactCollisionAlgorithm.h>
#include <BulletCollision/CollisionDispatch/btDefaultCollisionConfiguration.h>
#include <BulletCollision/BroadphaseCollision/btDbvtBroadphase.h>

#include <SiconosPointers.hpp>
DEFINE_SPTR(btCollisionWorld);
DEFINE_SPTR(btCollisionObject);
DEFINE_SPTR(btSphereShape);
DEFINE_SPTR(btDefaultCollisionConfiguration);
DEFINE_SPTR(btCollisionDispatcher);
DEFINE_SPTR(btBroadphaseInterface);

#include <stdio.h>

class BulletBroadphase;

class BulletShapeHandler : public SiconosShapeHandler
{
protected:
  // Non-owning pointer back to implementation,
  // to avoid circular reference count.
  BulletBroadphase *impl;

public:
  BulletShapeHandler(BulletBroadphase *_impl)
    : impl(_impl) {}

  virtual void onChanged(SP::SiconosSphere sphere);
  virtual void onChanged(SP::SiconosBox box);
};

class BulletBroadphase : public SiconosBroadphase, public std11::enable_shared_from_this<BulletBroadphase>
{
protected:
  SP::btCollisionWorld _collisionWorld;
  SP::btDefaultCollisionConfiguration _collisionConfiguration;
  SP::btCollisionDispatcher _dispatcher;
  SP::btBroadphaseInterface _broadphase;

  std::vector<SP::SiconosSphere> dirtySpheres;
  std::vector<SP::SiconosBox> dirtyBoxes;

  SP::Contactor currentContactor;

  SP::SiconosShapeHandler handler;

  friend class BulletShapeHandler;

public:
  BulletBroadphase() {
    handler.reset(new BulletShapeHandler(this));

    _collisionConfiguration.reset(new btDefaultCollisionConfiguration());

    _dispatcher.reset(new btCollisionDispatcher(&*_collisionConfiguration));

    _broadphase.reset(new btDbvtBroadphase());

    _collisionWorld.reset(new btCollisionWorld(&*_dispatcher, &*_broadphase,
                                               &*_collisionConfiguration));
  }

  ~BulletBroadphase() {
    // must be the first de-allocated, otherwise segfault
    _collisionWorld.reset();
  }

protected:
  virtual void visit(SP::SiconosSphere sphere)
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

  virtual void visit(SP::Contactor contactor)
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

public:
  virtual void buildGraph(SP::Contactor contactor)
    { contactor->acceptSP(shared_from_this()); }

  virtual void updateGraph()
  {
    printf("dirty spheres:\n");
    std::vector<SP::SiconosSphere>::iterator it;
    for (it=dirtySpheres.begin(); it!=dirtySpheres.end(); it++) {
      printf("  %p\n", &**it);
    }
    dirtySpheres.clear();

    printf("dirty boxes:\n");
    std::vector<SP::SiconosBox>::iterator it2;
    for (it2=dirtyBoxes.begin(); it2!=dirtyBoxes.end(); it2++) {
      printf("  %p\n", &**it2);
    }
    dirtyBoxes.clear();
  }

  virtual void performBroadphase()
    { _collisionWorld->performDiscreteCollisionDetection(); }
};

inline void BulletShapeHandler::onChanged(SP::SiconosSphere sphere) {
  impl->dirtySpheres.push_back(sphere);
  printf("pushed sphere: %p(%ld)\n",
         &*sphere,sphere.use_count());
}
inline void BulletShapeHandler::onChanged(SP::SiconosBox box) {
    impl->dirtyBoxes.push_back(box);
}

#endif /* BulletBroadphase.hpp */
