/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

/*! \file SiconosBulletCollisionManager_impl.hpp
  \brief Internal details for the definition of a Bullet-based interaction handler for contact
  detection.
*/

#ifndef SICOBULLETCOLLISIONIMPL_H
#define SICOBULLETCOLLISIONIMPL_H

#include <BulletCollision/BroadphaseCollision/btBroadphaseInterface.h>
#include <BulletCollision/CollisionDispatch/btCollisionWorld.h>
#include <BulletCollision/CollisionDispatch/btDefaultCollisionConfiguration.h>
#include <BulletCollision/CollisionShapes/btCapsuleShape.h>
#include <BulletCollision/CollisionShapes/btConeShape.h>
#include <BulletCollision/CollisionShapes/btConvex2dShape.h>
#include <BulletCollision/CollisionShapes/btConvexHullShape.h>
#include <BulletCollision/CollisionShapes/btCylinderShape.h>
#include <BulletCollision/CollisionShapes/btSphereShape.h>

#include <map>
#include <variant>

#include "BodyBulletShapeRecord.hpp"
#include "SiconosBulletDefines.h"
#include "SiconosContactor.hpp"
#include "SiconosShape.hpp"

namespace siconos::modeling {
class SecondOrderDS;
class LagrangianDS;
class NewtonEulerDS;
}  // namespace siconos::modeling

namespace siconos::simulation {
class Simulation;
class InteractionManager;
}  // namespace siconos::simulation

namespace siconos::collision {
class RigidBodyDS;
class RigidBody2dDS;
class StaticBody;
}  // namespace siconos::collision

namespace siconos::collision::bullet {
class SiconosBulletOptions;
class SiconosBulletCollisionManager;
}  // namespace siconos::collision::bullet

namespace siconos::collision::bullet::internal {

class SiconosMeshData;
class SiconosHeightData;

using BodyBoxRecord =
    BodyRecord<siconos::collision::RigidBodyDS, siconos::collision::SiconosBox, BTBOXSHAPE>;

using BodySphereRecord = BodyRecord<siconos::collision::RigidBodyDS,
                                    siconos::collision::SiconosSphere, BTSPHERESHAPE>;
using BodyCHRecord = BodyRecord<siconos::collision::RigidBodyDS,
                                siconos::collision::SiconosConvexHull, btConvexHullShape>;
using BodyPlaneRecord = BodyRecord<siconos::collision::RigidBodyDS,
                                   siconos::collision::SiconosPlane, BTPLANESHAPE>;
using BodyCylinderRecord = BodyRecord<siconos::collision::RigidBodyDS,
                                      siconos::collision::SiconosCylinder, btCylinderShape>;
using BodyConeRecord =
    BodyRecord<siconos::collision::RigidBodyDS, siconos::collision::SiconosCone, btConeShape>;
using BodyCapsuleRecord = BodyRecord<siconos::collision::RigidBodyDS,
                                     siconos::collision::SiconosCapsule, btCapsuleShape>;
using BodyMeshRecord =
    BodyRecord<siconos::collision::RigidBodyDS, siconos::collision::SiconosMesh,
               siconos::collision::bullet::internal::SiconosMeshData>;
using BodyHeightRecord =
    BodyRecord<siconos::collision::RigidBodyDS, siconos::collision::SiconosHeightMap,
               siconos::collision::bullet::internal::SiconosHeightData>;

using BodyBox2dRecord = BodyRecord<siconos::collision::RigidBody2dDS,
                                   siconos::collision::SiconosBox2d, btConvex2dShape>;
using BodyCH2dRecord = BodyRecord<siconos::collision::RigidBody2dDS,
                                  siconos::collision::SiconosConvexHull2d, btConvex2dShape>;
using BodyDiskRecord = BodyRecord<siconos::collision::RigidBody2dDS,
                                  siconos::collision::SiconosDisk, btConvex2dShape>;

// A variant used to visit records during updateAllShapesForDS call
using RecordVariant =
    std::variant<std::shared_ptr<BodyPlaneRecord>, std::shared_ptr<BodySphereRecord>,
                 std::shared_ptr<BodyBoxRecord>, std::shared_ptr<BodyCylinderRecord>,
                 std::shared_ptr<BodyCapsuleRecord>, std::shared_ptr<BodyConeRecord>,
                 std::shared_ptr<BodyCH2dRecord>, std::shared_ptr<BodyCHRecord>,
                 std::shared_ptr<BodyMeshRecord>, std::shared_ptr<BodyHeightRecord>,
                 std::shared_ptr<BodyDiskRecord>, std::shared_ptr<BodyBox2dRecord>>;

using BodyShapeMap =
    std::map<const siconos::modeling::SecondOrderDS *, std::vector<RecordVariant>>;

using StaticBodyShapeMap = std::map<const siconos::collision::StaticBody *,
                                    std::vector<std::shared_ptr<BodyBulletShapeRecord>>>;

/* We derive a specific callback for filtering the broadphase of Bullet
 * based on collision group */
struct SiconosBulletFilterCallback : public btOverlapFilterCallback {
  siconos::simulation::InteractionManager *interactionManager{nullptr};
  // return true when pairs need collision
  virtual bool needBroadphaseCollision(btBroadphaseProxy *proxy0,
                                       btBroadphaseProxy *proxy1) const override;
};

class SiconosBulletCollisionManager_impl
    : public std::enable_shared_from_this<SiconosBulletCollisionManager_impl> {
 protected:
  std::shared_ptr<btCollisionWorld> _collisionWorld{nullptr};
  std::shared_ptr<btDefaultCollisionConfiguration> _collisionConfiguration{nullptr};
  std::shared_ptr<btCollisionDispatcher> _dispatcher{nullptr};
  std::shared_ptr<btBroadphaseInterface> _broadphase{nullptr};

  /* During iteration over DSs for position updates we need to access
   * btCollisionObject, so need a map DS->btXShape. */
  BodyShapeMap bodyShapeMap;

  StaticBodyShapeMap staticBodyShapeMap;

  std::weak_ptr<siconos::simulation::Simulation> _simulation{nullptr};

  std::shared_ptr<siconos::collision::bullet::SiconosBulletOptions> _options{nullptr};

  std::vector<std::pair<std::shared_ptr<btCollisionObject>, int>> _queuedCollisionObjects;

  // Rule of five
  SiconosBulletCollisionManager_impl() = delete;
  SiconosBulletCollisionManager_impl(const SiconosBulletCollisionManager_impl &) = delete;
  SiconosBulletCollisionManager_impl(SiconosBulletCollisionManager_impl &&) = delete;
  SiconosBulletCollisionManager_impl &operator=(const SiconosBulletCollisionManager_impl &) =
      delete;
  SiconosBulletCollisionManager_impl &operator=(SiconosBulletCollisionManager_impl &&) =
      delete;

  // Create collision objects for each shape type
  friend class CreateCollisionObjectShapeVisitor;

  // To visit RigidBody DS
  friend class RigidBodyDSVisitor;

  // To allow access to protected members ...
  friend class siconos::collision::bullet::SiconosBulletCollisionManager;

  // Generic version to ensure that all combination (DS / SHAPE) in the variants
  // are taken into account. This is required for std::visit
  // This function will return an exception.
  // A specialized function must be called (see createCollisionobject functions below, with
  // specific ds and shapes).
  template <typename DS, typename SHAPE>
  void createCollisionObject(const std::shared_ptr<siconos::algebra::SiconosVector> base,
                             const std::shared_ptr<DS> ds, std::shared_ptr<SHAPE> plane,
                             std::shared_ptr<SiconosContactor> contactor,
                             const std::shared_ptr<StaticBody> staticBody);

  /* Create collision objects for each shape type */
  void createCollisionObject(
      const std::shared_ptr<siconos::algebra::SiconosVector> base,
      const std::shared_ptr<siconos::collision::RigidBodyDS> ds,
      const std::shared_ptr<siconos::collision::SiconosPlane> plane,
      const std::shared_ptr<siconos::collision::SiconosContactor> contactor,
      const std::shared_ptr<siconos::collision::StaticBody> staticBody);
  void createCollisionObject(
      const std::shared_ptr<siconos::algebra::SiconosVector> base,
      const std::shared_ptr<siconos::collision::RigidBodyDS> ds,
      const std::shared_ptr<siconos::collision::SiconosSphere> sphere,
      const std::shared_ptr<siconos::collision::SiconosContactor> contactor,
      const std::shared_ptr<siconos::collision::StaticBody> staticBody);
  void createCollisionObject(
      const std::shared_ptr<siconos::algebra::SiconosVector> base,
      const std::shared_ptr<siconos::collision::RigidBodyDS> ds,
      const std::shared_ptr<siconos::collision::SiconosBox> box,
      const std::shared_ptr<siconos::collision::SiconosContactor> contactor,
      const std::shared_ptr<siconos::collision::StaticBody> staticBody);
  void createCollisionObject(
      const std::shared_ptr<siconos::algebra::SiconosVector> base,
      const std::shared_ptr<siconos::collision::RigidBodyDS> ds,
      const std::shared_ptr<siconos::collision::SiconosCylinder> cyl,
      const std::shared_ptr<siconos::collision::SiconosContactor> contactor,
      const std::shared_ptr<siconos::collision::StaticBody> staticBody);
  void createCollisionObject(
      const std::shared_ptr<siconos::algebra::SiconosVector> base,
      const std::shared_ptr<siconos::collision::RigidBodyDS> ds,
      const std::shared_ptr<siconos::collision::SiconosCone> cone,
      const std::shared_ptr<siconos::collision::SiconosContactor> contactor,
      const std::shared_ptr<siconos::collision::StaticBody> staticBody);
  void createCollisionObject(
      const std::shared_ptr<siconos::algebra::SiconosVector> base,
      const std::shared_ptr<siconos::collision::RigidBodyDS> ds,
      const std::shared_ptr<siconos::collision::SiconosCapsule> capsule,
      const std::shared_ptr<siconos::collision::SiconosContactor> contactor,
      const std::shared_ptr<siconos::collision::StaticBody> staticBody);
  void createCollisionObject(
      const std::shared_ptr<siconos::algebra::SiconosVector> base,
      const std::shared_ptr<siconos::collision::RigidBodyDS> ds,
      const std::shared_ptr<siconos::collision::SiconosConvexHull> ch,
      const std::shared_ptr<siconos::collision::SiconosContactor> contactor,
      const std::shared_ptr<siconos::collision::StaticBody> staticBody);
  void createCollisionObject(
      const std::shared_ptr<siconos::algebra::SiconosVector> base,
      const std::shared_ptr<siconos::collision::RigidBodyDS> ds,
      const std::shared_ptr<siconos::collision::SiconosMesh> mesh,
      const std::shared_ptr<siconos::collision::SiconosContactor> contactor,
      const std::shared_ptr<siconos::collision::StaticBody> staticBody);
  void createCollisionObject(
      const std::shared_ptr<siconos::algebra::SiconosVector> base,
      const std::shared_ptr<siconos::collision::RigidBodyDS> ds,
      const std::shared_ptr<siconos::collision::SiconosHeightMap> height,
      const std::shared_ptr<siconos::collision::SiconosContactor> contactor,
      const std::shared_ptr<siconos::collision::StaticBody> staticBody);

  void createCollisionObject(
      const std::shared_ptr<siconos::algebra::SiconosVector> base,
      const std::shared_ptr<siconos::collision::RigidBody2dDS> ds,
      const std::shared_ptr<siconos::collision::SiconosDisk> disk,
      const std::shared_ptr<siconos::collision::SiconosContactor> contactor,
      const std::shared_ptr<siconos::collision::StaticBody> staticBody);
  void createCollisionObject(
      const std::shared_ptr<siconos::algebra::SiconosVector> base,
      const std::shared_ptr<siconos::collision::RigidBody2dDS> ds,
      const std::shared_ptr<siconos::collision::SiconosBox2d> box2d,
      const std::shared_ptr<siconos::collision::SiconosContactor> contactor,
      const std::shared_ptr<siconos::collision::StaticBody> staticBody);
  void createCollisionObject(
      const std::shared_ptr<siconos::algebra::SiconosVector> base,
      const std::shared_ptr<siconos::collision::RigidBody2dDS> ds,
      const std::shared_ptr<siconos::collision::SiconosConvexHull2d> ch2d,
      const std::shared_ptr<siconos::collision::SiconosContactor> contactor,
      const std::shared_ptr<siconos::collision::StaticBody> staticBody);

  /* Call the above functions for each shape associated with a body or contactor. */
  void createCollisionObjectsForStaticBodyContactorSet(
      const std::shared_ptr<StaticBody> staticBody,
      const std::shared_ptr<const siconos::collision::SiconosContactorSet> contactor);

  void createCollisionObjectsForBodyContactorSetFromDS(
      const std::shared_ptr<siconos::modeling::SecondOrderDS> ds);

  /* A helper function used to initialise new shapes, generic to the
   * shape type */
  template <typename ST, typename BT, typename DST, typename BR>
  std::shared_ptr<btCollisionObject> createCollisionObjectHelper(
      std::shared_ptr<siconos::algebra::SiconosVector> base, const std::shared_ptr<DST> ds,
      std::shared_ptr<ST> shape, std::shared_ptr<BT> btshape, BodyShapeMap &bodyShapeMap,
      std::shared_ptr<SiconosContactor> contactor, StaticBodyShapeMap &StaticBodyShapeMap,
      std::shared_ptr<StaticBody> staticBody);
  void updateShape(std::shared_ptr<BodySphereRecord> record);
  void updateShape(std::shared_ptr<BodyPlaneRecord> record);
  void updateShape(std::shared_ptr<BodyBoxRecord> record);
  void updateShape(std::shared_ptr<BodyCylinderRecord> record);
  void updateShape(std::shared_ptr<BodyConeRecord> record);
  void updateShape(std::shared_ptr<BodyCapsuleRecord> record);
  void updateShape(std::shared_ptr<BodyCHRecord> record);
  void updateShape(std::shared_ptr<BodyMeshRecord> record);
  void updateShape(std::shared_ptr<BodyHeightRecord> record);

  void updateShape(std::shared_ptr<BodyDiskRecord> record);
  void updateShape(std::shared_ptr<BodyBox2dRecord> record);
  void updateShape(std::shared_ptr<BodyCH2dRecord> record);

  void updateAllShapesForDS(const siconos::modeling::SecondOrderDS &bds);
  void updateShapePosition(std::shared_ptr<BodyBulletShapeRecord> record);

  /* Helper to apply an offset transform to a position and return as a
   * btTransform */
  btTransform offsetTransform(const siconos::algebra::SiconosVector &position,
                              const siconos::algebra::SiconosVector &offset);

  btTransform offsetTransform(const siconos::algebra::SiconosVector &position);

  /** Helper to set the inertia of a NewtonEulerDS based on a
   * btCollisionShape */
  void updateContactorInertia(std::shared_ptr<siconos::modeling::NewtonEulerDS> ds,
                              std::shared_ptr<btCollisionShape> btshape);

  /** Helper to set the inertia of a LagrangianDS based on a
   * btCollisionShape */
  void update2DContactorInertia(std::shared_ptr<siconos::modeling::LagrangianDS> ds,
                                std::shared_ptr<btCollisionShape> btshape);

 public:
  SiconosBulletCollisionManager_impl(std::shared_ptr<SiconosBulletOptions> op)
      : _options{op} {};

  ~SiconosBulletCollisionManager_impl() noexcept = default;

  template <class DS>
  void updateShapes(std::shared_ptr<DS> ds) {
    auto bds = std::dynamic_pointer_cast<DS>(ds);
    if (bds->contactors()) {
      if (bodyShapeMap.find(&*bds) == bodyShapeMap.end()) {
        createCollisionObjectsForBodyContactorSetFromDS(bds);
      }
      updateAllShapesForDS(*bds);
    }
  }
};

}  // namespace siconos::collision::bullet::internal

#endif
