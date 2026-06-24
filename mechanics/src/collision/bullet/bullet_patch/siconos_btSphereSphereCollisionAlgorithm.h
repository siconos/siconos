/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2006 Erwin Coumans  https://bulletphysics.org

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this
software. Permission is granted to anyone to use this software for any purpose, including
commercial applications, and to alter it and redistribute it freely, subject to the following
restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote
the original software. If you use this software in a product, an acknowledgment in the product
documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as
being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

#ifndef SICONOS_BT_SPHERE_SPHERE_COLLISION_ALGORITHM_H
#define SICONOS_BT_SPHERE_SPHERE_COLLISION_ALGORITHM_H

#include "BulletCollision/BroadphaseCollision/btBroadphaseProxy.h"
#include "BulletCollision/CollisionDispatch/btActivatingCollisionAlgorithm.h"
#include "BulletCollision/CollisionDispatch/btCollisionCreateFunc.h"
#include "BulletCollision/CollisionDispatch/btCollisionDispatcher.h"

class btPersistentManifold;

/// siconos_btSphereSphereCollisionAlgorithm  provides sphere-sphere collision detection.
/// Other features are frame-coherency (persistent data) and collision response.
/// Also provides the most basic sample for custom/user btCollisionAlgorithm
class siconos_btSphereSphereCollisionAlgorithm : public btActivatingCollisionAlgorithm {
  bool m_ownManifold;
  btPersistentManifold* m_manifoldPtr;

 public:
  siconos_btSphereSphereCollisionAlgorithm(btPersistentManifold* mf,
                                           const btCollisionAlgorithmConstructionInfo& ci,
                                           const btCollisionObjectWrapper* col0Wrap,
                                           const btCollisionObjectWrapper* col1Wrap);

  siconos_btSphereSphereCollisionAlgorithm(const btCollisionAlgorithmConstructionInfo& ci)
      : btActivatingCollisionAlgorithm(ci) {}

  virtual void processCollision(const btCollisionObjectWrapper* body0Wrap,
                                const btCollisionObjectWrapper* body1Wrap,
                                const btDispatcherInfo& dispatchInfo,
                                btManifoldResult* resultOut);

  virtual btScalar calculateTimeOfImpact(btCollisionObject* body0, btCollisionObject* body1,
                                         const btDispatcherInfo& dispatchInfo,
                                         btManifoldResult* resultOut);

  virtual void getAllContactManifolds(btManifoldArray& manifoldArray) {
    if (m_manifoldPtr && m_ownManifold) {
      manifoldArray.push_back(m_manifoldPtr);
    }
  }

  virtual ~siconos_btSphereSphereCollisionAlgorithm();

  struct CreateFunc : public btCollisionAlgorithmCreateFunc {
    virtual btCollisionAlgorithm* CreateCollisionAlgorithm(
        btCollisionAlgorithmConstructionInfo& ci, const btCollisionObjectWrapper* col0Wrap,
        const btCollisionObjectWrapper* col1Wrap) {
      void* mem = ci.m_dispatcher1->allocateCollisionAlgorithm(
          sizeof(siconos_btSphereSphereCollisionAlgorithm));
      return new (mem)
          siconos_btSphereSphereCollisionAlgorithm(nullptr, ci, col0Wrap, col1Wrap);
    }
  };
};

#endif  // BT_SPHERE_SPHERE_COLLISION_ALGORITHM_H
