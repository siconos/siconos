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

#include "SiconosBulletCollisionManager_impl.hpp"

#include <BulletCollision/BroadphaseCollision/btAxisSweep3.h>
#include <BulletCollision/BroadphaseCollision/btDbvtBroadphase.h>
#include <BulletCollision/CollisionShapes/btConvexHullShape.h>
#include <BulletCollision/CollisionShapes/btTriangleIndexVertexArray.h>
#include <LinearMath/btConvexHullComputer.h>

#include <vector>

#include "BodyBulletShapeRecord.hpp"
#include "RigidBody2dDS.hpp"
#include "RigidBodyDS.hpp"
#include "SecondOrderDS.hpp"
#include "SiconosBulletDefines.h"
#include "SiconosBulletOptions.hpp"
#include "SiconosBulletShape.hpp"
#include "SiconosBulletVisitors.hpp"
#include "SiconosCollisionManager.hpp"
#include "SiconosContactor.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "StaticBody.hpp"
#include "siconos_debug.h"

namespace {  // anonymous because only for local use
             // ======================= Helper functions
             // ========================== Only used and available in this file
             // --> do not declare them in user interface

// helper to enable polyhedral contact clipping for shape types
// derived from btPolyhedralConvexShape
void initPolyhedralFeatures(btPolyhedralConvexShape& btshape) {
  btshape.initializePolyhedralFeatures();
}

void initPolyhedralFeatures(btCollisionShape& btshape) {}

int find_index_closest_point_btConvexHullShape(btVector3& pointA, btConvexHullShape& btch) {
  int numPoints = btch.getNumPoints();
  const btVector3* points = btch.getPoints();
  btScalar min_dist = 1e30;
  int p_idx = -1;
  for (int p = 0; p < numPoints; p++) {
    btScalar l2 = (points[p] - pointA).length2();
    if (l2 < min_dist) {
      min_dist = l2;
      p_idx = p;
    }
  }
  return p_idx;
}

// If type of SiconosMatrix is the same as btScalar, we can avoid a copy
template <typename SCALAR>
std::pair<std::shared_ptr<btTriangleIndexVertexArray>, SCALAR*> make_bt_vertex_array(
    std::shared_ptr<siconos::collision::SiconosMesh> mesh, SCALAR _s1, SCALAR _s2) {
  assert(mesh->vertices()->rows() == 3);
  auto bttris = std::make_shared<btTriangleIndexVertexArray>(
      mesh->indexes()->size() / 3, (int*)mesh->indexes()->data(), sizeof(int) * 3,
      mesh->vertices()->cols(), mesh->vertices()->data(), sizeof(btScalar) * 3);

  return std::make_pair(bttris, (btScalar*)nullptr);
}

// If type of SiconosMatrix is not the same as btScalar, we must copy
template <typename SCALAR1, typename SCALAR2>
std::pair<std::shared_ptr<btTriangleIndexVertexArray>, btScalar*> make_bt_vertex_array(
    std::shared_ptr<siconos::collision::SiconosMesh> mesh, SCALAR1 _s1, SCALAR2 _s2) {
  assert(mesh->vertices()->rows() == 3);
  auto numIndices = mesh->indexes()->size();
  auto numVertices = mesh->vertices()->cols();
  btScalar* vertices = new btScalar[numVertices * 3];
  for (siconos::algebra::Index i = 0; i < numVertices; i++) {
    vertices[i * 3 + 0] = (*mesh->vertices())(0, i);
    vertices[i * 3 + 1] = (*mesh->vertices())(1, i);
    vertices[i * 3 + 2] = (*mesh->vertices())(2, i);
  }
  auto bttris = std::make_shared<btTriangleIndexVertexArray>(
      numIndices / 3, (int*)mesh->indexes()->data(), sizeof(int) * 3, numVertices, vertices,
      sizeof(btScalar) * 3);
  return std::make_pair(bttris, vertices);
}
}  // namespace

bool siconos::collision::bullet::internal::SiconosBulletFilterCallback::
    needBroadphaseCollision(btBroadphaseProxy* proxy0, btBroadphaseProxy* proxy1) const {
  DEBUG_BEGIN("SiconosBulletFilterCallback :: needBroadphaseCollision\n");

  /* standard filter in Bullet */
  // bool collides = (proxy0->m_collisionFilterGroup &
  // proxy1->m_collisionFilterMask) != 0; collides = collides &&
  // (proxy1->m_collisionFilterGroup & proxy0->m_collisionFilterMask);

  // add some additional logic here that modified 'collides'
  auto nslaw = interactionManager->nonSmoothLaw(proxy0->m_collisionFilterGroup,
                                                proxy1->m_collisionFilterGroup);

  bool collides = (bool)nslaw;

  DEBUG_END("SiconosBulletFilterCallback :: needBroadphaseCollision\n");
  return collides;
}

void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::
    updateAllShapesForDS(const siconos::modeling::SecondOrderDS& bds) {
  for (auto& it : bodyShapeMap[&bds])
    std::visit([this](auto rec) { this->updateShape(rec); }, it);
}

template <typename ST, typename BT, typename DST, typename BR>
std::shared_ptr<btCollisionObject> siconos::collision::bullet::internal::
    SiconosBulletCollisionManager_impl::createCollisionObjectHelper(
        std::shared_ptr<siconos::algebra::SiconosVector> base, const std::shared_ptr<DST> ds,
        std::shared_ptr<ST> shape, std::shared_ptr<BT> btshape, BodyShapeMap& bodyShapeMap,
        std::shared_ptr<SiconosContactor> contactor, StaticBodyShapeMap& StaticBodyShapeMap,
        std::shared_ptr<StaticBody> staticBody) {
  assert(base && "Collision objects must have a base position.");

  // create corresponding Bullet object and shape
  auto btobject = std::make_shared<btCollisionObject>();

  // default parameters
  btobject->setContactProcessingThreshold(_options->contactProcessingThreshold);

  // associate the shape with the object
  btobject->setCollisionShape(&*btshape);

  // enable contact clipping for SAT
  if (_options->enablePolyhedralContactClipping) initPolyhedralFeatures(*btshape);

  if (!ds)
    btobject->setCollisionFlags(btCollisionObject::CF_STATIC_OBJECT);
  else
    btobject->setCollisionFlags(btCollisionObject::CF_KINEMATIC_OBJECT);

  // put it in the world
  int collisionFilterGroup = contactor->collision_group;
  int collisionFilterMask = 1;

#ifdef QUEUE_STATIC_CONTACTORS
  if (!ds) {
    _queuedCollisionObjects.push_back(std::make_pair(btobject, collisionFilterGroup));
  } else
    _collisionWorld->addCollisionObject(&*btobject, collisionFilterGroup, collisionFilterMask);
#else
  _collisionWorld->addCollisionObject(&*btobject, collisionFilterGroup, collisionFilterMask);
#endif

  // create a record to keep track of things
  // (for static contactor, ds=nil)
  auto record =
      std::make_shared<BR>(base, ds, shape, btshape, btobject, contactor, staticBody);

  bodyShapeMap[ds ? &*ds : nullptr].push_back(record);

  if (staticBody) StaticBodyShapeMap[&*staticBody].push_back(record);

  assert(record->btobject);
  assert(record->sshape);
  assert(record->shape);
  assert(record->btshape);
  assert(record->contactor);
  // assert(record->contactor->offset);
  // assert(record->contactor->offset->size() == 7);

  // Allow Bullet to report colliding DSs.  We need to access it from
  // the collision callback as the record base class so down-cast it.
  btobject->setUserPointer(
      reinterpret_cast<void*>(static_cast<siconos::collision::BodyShapeRecord*>(&*record)));

  // initial parameter update (change version to make something happen)
  record->shape_version -= 1;
  updateShape(record);

  return btobject;
}

template <typename DS, typename SHAPE>
void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::
    createCollisionObject(const std::shared_ptr<siconos::algebra::SiconosVector> base,
                          const std::shared_ptr<DS> ds, std::shared_ptr<SHAPE> plane,
                          std::shared_ptr<SiconosContactor> contactor,
                          const std::shared_ptr<StaticBody> staticBody) {
  THROW_EXCEPTION("Undefined creation process for this combination of DS and shape.");
}

void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::
    createCollisionObject(const std::shared_ptr<siconos::algebra::SiconosVector> base,
                          const std::shared_ptr<siconos::collision::RigidBodyDS> ds,
                          std::shared_ptr<siconos::collision::SiconosPlane> plane,
                          std::shared_ptr<siconos::collision::SiconosContactor> contactor,
                          const std::shared_ptr<siconos::collision::StaticBody> staticBody) {
  // create the initial plane with default parameters
#ifdef USE_BOX_FOR_PLANE
  btScalar h = (1000 + plane->outsideMargin()) * _options->worldScale;
  auto btplane = std::make_shared<btBoxShape>(btVector3(h, h, h));

  // Internal margin
  btplane->setMargin((plane->insideMargin() + plane->outsideMargin()) * _options->worldScale);
#else
#ifdef USE_CONVEXHULL_FOR_PLANE
  btScalar h = 1000 * _options->worldScale;
  const btScalar pts[] = {
      h, h, 0, h, -h, 0, -h, -h, 0, -h, h, 0,
  };
  auto btplane = std::make_shared<btConvexHullShape>(pts, 4, sizeof(pts[0]) * 3);

  // We ignore the desired internal margin for plane and just use a large one.
  plane->setInsideMargin(1000 * _options->worldScale);

  // External margin
  btplane->setMargin((plane->insideMargin() + plane->outsideMargin()) * _options->worldScale);
#else
  auto btplane = std::make_shared<btStaticPlaneShape>(btVector3(0, 0, 1), 0.0);
  btplane->setMargin(plane->outsideMargin() * _options->worldScale);
#endif
#endif

  // initialization
  createCollisionObjectHelper<siconos::collision::SiconosPlane, BTPLANESHAPE,
                              siconos::collision::RigidBodyDS, BodyPlaneRecord>(
      base, ds, plane, btplane, bodyShapeMap, contactor, staticBodyShapeMap, staticBody);
}

void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::updateShape(
    std::shared_ptr<BodyPlaneRecord> record) {
  auto plane = record->shape;
  auto btplane = record->btshape;

  auto o = record->contactor->offset;

  // Adjust the offset according to plane implementation
#ifdef USE_BOX_FOR_PLANE
  o(2) -= -plane->outsideMargin() + 1000;
#else
#ifdef USE_CONVEXHULL_FOR_PLANE
  o(2) -= plane->insideMargin();
#else  // USE_PLANE_FOR_PLANE
  o(2) -= plane->insideMargin();
#endif
#endif

  // Note, we do not use generic updateShapePosition for plane
  auto t = offsetTransform(*record->base, o);
  t.setOrigin(t.getOrigin() * _options->worldScale);
  record->btobject->setWorldTransform(t);
}

void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::
    createCollisionObject(
        const std::shared_ptr<siconos::algebra::SiconosVector> base,
        const std::shared_ptr<siconos::collision::RigidBodyDS> ds,
        std::shared_ptr<siconos::collision::SiconosSphere> sphere,
        const std::shared_ptr<siconos::collision::SiconosContactor> contactor,
        const std::shared_ptr<siconos::collision::StaticBody> staticBody) {
  // set radius to 1.0 and use scaling instead of setting radius
  // directly, makes it easier to change during update

#ifdef USE_CONVEXHULL_FOR_SPHERE
  // A sphere can be represented as a convex hull of a single point, with the
  // margin equal to the radius size
  auto btsphere = std::make_shared<btConvexHullShape>();
  {
    btsphere->addPoint(btVector3(0.0, 0.0, 0.0));
    btsphere->setMargin(0.0);
  }
#else
  auto btsphere = std::make_shared<btSphereShape>(1.0);
  btsphere->setMargin(0.0);
#endif

  // initialization
  createCollisionObjectHelper<siconos::collision::SiconosSphere, BTSPHERESHAPE,
                              siconos::collision::RigidBodyDS, BodySphereRecord>(
      base, ds, sphere, btsphere, bodyShapeMap, contactor, staticBodyShapeMap, staticBody);
}

void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::updateShape(
    std::shared_ptr<BodySphereRecord> record) {
  auto sphere = record->shape;
  auto btsphere = record->btshape;

  if (sphere->version() != record->shape_version) {
    double r = (sphere->radius() + sphere->outsideMargin()) * _options->worldScale;

    // Update shape parameters
#ifdef USE_CONVEXHULL_FOR_SPHERE
    btsphere->setMargin(r);
#else
    btsphere->setLocalScaling(btVector3(r, r, r));

    // btSphereShape has an internal margin
    btsphere->setMargin(sphere->insideMargin() * _options->worldScale);
#endif
    auto rbds = std::static_pointer_cast<RigidBodyDS>(record->ds);
    if (record->ds && rbds->useContactorInertia()) {
      updateContactorInertia(rbds, btsphere);
    }

    if (record->btobject->getBroadphaseHandle()) {
      _collisionWorld->updateSingleAabb(&*record->btobject);
      //      _collisionWorld->getBroadphase()->getOverlappingPairCache()->
      //        cleanProxyFromPairs(record->btobject->getBroadphaseHandle(),
      //        &*_dispatcher);
    }

    record->shape_version = sphere->version();
  }

  updateShapePosition(record);
}

void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::
    createCollisionObject(const std::shared_ptr<siconos::algebra::SiconosVector> base,
                          const std::shared_ptr<siconos::collision::RigidBodyDS> ds,
                          std::shared_ptr<siconos::collision::SiconosBox> box,
                          std::shared_ptr<siconos::collision::SiconosContactor> contactor,
                          const std::shared_ptr<siconos::collision::StaticBody> staticBody) {
  const btScalar half = 0.5;

#ifdef USE_CONVEXHULL_FOR_BOX
  const btScalar pts[] = {
      -half, half, -half, -half, -half, -half, -half, -half, half,  -half, half,  half,
      half,  half, half,  half,  half,  -half, half,  -half, -half, half,  -half, half,
  };
  auto btbox = std::make_shared<btConvexHullShape>(pts, 8, sizeof(pts[0]) * 3);

  // External margin
  btbox->setMargin(box->outsideMargin());
#else
  auto btbox = std::make_shared<btBoxShape>(btVector3(half, half, half));

  // btBoxShape has an internal margin
  btbox->setMargin(box->insideMargin() * _options->worldScale);
#endif

  // initialization
  createCollisionObjectHelper<siconos::collision::SiconosBox, BTBOXSHAPE,
                              siconos::collision::RigidBodyDS, BodyBoxRecord>(
      base, ds, box, btbox, bodyShapeMap, contactor, staticBodyShapeMap, staticBody);
}

void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::updateShape(
    std::shared_ptr<BodyBoxRecord> record) {
  auto box = record->shape;
  auto btbox = record->btshape;

  // Update shape parameters
  if (box->version() != record->shape_version) {
#ifdef USE_CONVEXHULL_FOR_BOX
    double m = -box->insideMargin();
#else
    double m = box->outsideMargin();
#endif

    double sx = ((*box->dimensions())(0) + m * 2) * _options->worldScale;
    double sy = ((*box->dimensions())(1) + m * 2) * _options->worldScale;
    double sz = ((*box->dimensions())(2) + m * 2) * _options->worldScale;

    assert(sx > 0 && sy > 0 && sz > 0);

    btbox->setLocalScaling(btVector3(sx, sy, sz));
    btbox->setMargin((box->insideMargin() + box->outsideMargin()) * _options->worldScale);

    auto rbds = std::static_pointer_cast<RigidBodyDS>(record->ds);
    if (record->ds && rbds->useContactorInertia()) updateContactorInertia(rbds, btbox);

    if (record->btobject->getBroadphaseHandle()) {
      _collisionWorld->updateSingleAabb(&*record->btobject);
      //      _collisionWorld->getBroadphase()->getOverlappingPairCache()->
      //        cleanProxyFromPairs(record->btobject->getBroadphaseHandle(),
      //        &*_dispatcher);
    }

    record->shape_version = box->version();
  }

  updateShapePosition(record);
}

void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::
    createCollisionObject(const std::shared_ptr<siconos::algebra::SiconosVector> base,
                          const std::shared_ptr<siconos::collision::RigidBodyDS> ds,
                          std::shared_ptr<siconos::collision::SiconosCylinder> cylinder,
                          std::shared_ptr<siconos::collision::SiconosContactor> contactor,
                          const std::shared_ptr<siconos::collision::StaticBody> staticBody) {
  auto btcylinder = std::make_shared<btCylinderShape>(btVector3(1.0, 1.0, 1.0));

  // initialization
  createCollisionObjectHelper<siconos::collision::SiconosCylinder, btCylinderShape,
                              siconos::collision::RigidBodyDS, BodyCylinderRecord>(
      base, ds, cylinder, btcylinder, bodyShapeMap, contactor, staticBodyShapeMap, staticBody);
}

void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::updateShape(
    std::shared_ptr<BodyCylinderRecord> record) {
  auto cyl = record->shape;
  auto btcyl = record->btshape;

  // Update shape parameters
  if (cyl->version() != record->shape_version) {
    // Bullet cylinder has an inside margin, so we add the outside
    // margin explicitly.
    double m = cyl->outsideMargin();

    double radius = (cyl->radius() + m) * _options->worldScale;
    double length = (cyl->length() / 2 + m) * _options->worldScale;

    assert(radius > 0 && length > 0);

    btcyl->setLocalScaling(btVector3(radius, length, radius));
    btcyl->setMargin((cyl->insideMargin() + cyl->outsideMargin()) * _options->worldScale);
    auto rbds = std::static_pointer_cast<RigidBodyDS>(record->ds);
    if (record->ds && rbds->useContactorInertia()) updateContactorInertia(rbds, btcyl);

    if (record->btobject->getBroadphaseHandle()) {
      _collisionWorld->updateSingleAabb(&*record->btobject);
      //      _collisionWorld->getBroadphase()->getOverlappingPairCache()->
      //        cleanProxyFromPairs(record->btobject->getBroadphaseHandle(),
      //        &*_dispatcher);
    }

    record->shape_version = cyl->version();
  }

  updateShapePosition(record);
}

void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::
    createCollisionObject(const std::shared_ptr<siconos::algebra::SiconosVector> base,
                          const std::shared_ptr<siconos::collision::RigidBodyDS> ds,
                          std::shared_ptr<siconos::collision::SiconosCone> cone,
                          std::shared_ptr<siconos::collision::SiconosContactor> contactor,
                          const std::shared_ptr<siconos::collision::StaticBody> staticBody) {
  auto btcone = std::make_shared<btConeShape>(1.0, 1.0);

  // initialization
  createCollisionObjectHelper<siconos::collision::SiconosCone, btConeShape,
                              siconos::collision::RigidBodyDS, BodyConeRecord>(
      base, ds, cone, btcone, bodyShapeMap, contactor, staticBodyShapeMap, staticBody);
}

void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::updateShape(
    std::shared_ptr<BodyConeRecord> record) {
  auto cone = record->shape;
  auto btcone = record->btshape;

  // Update shape parameters
  if (cone->version() != record->shape_version) {
    // Bullet cone has an inside margin, so we add the outside
    // margin explicitly.
    double m = cone->outsideMargin();

    double radius = (cone->radius() + m) * _options->worldScale;
    double length = (cone->length() / 2 + m) * _options->worldScale;

    assert(radius > 0 && length > 0);

    btcone->setLocalScaling(btVector3(radius, length, radius));
    btcone->setMargin((cone->insideMargin() + cone->outsideMargin()) * _options->worldScale);
    auto rbds = std::static_pointer_cast<RigidBodyDS>(record->ds);
    if (record->ds && rbds->useContactorInertia()) updateContactorInertia(rbds, btcone);

    if (record->btobject->getBroadphaseHandle()) {
      _collisionWorld->updateSingleAabb(&*record->btobject);
      //      _collisionWorld->getBroadphase()->getOverlappingPairCache()->
      //        cleanProxyFromPairs(record->btobject->getBroadphaseHandle(),
      //        &*_dispatcher);
    }

    record->shape_version = cone->version();
  }

  updateShapePosition(record);
}

void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::
    createCollisionObject(const std::shared_ptr<siconos::algebra::SiconosVector> base,
                          const std::shared_ptr<siconos::collision::RigidBodyDS> ds,
                          std::shared_ptr<siconos::collision::SiconosCapsule> capsule,
                          std::shared_ptr<siconos::collision::SiconosContactor> contactor,
                          const std::shared_ptr<siconos::collision::StaticBody> staticBody) {
  auto btcapsule = std::make_shared<btCapsuleShape>(1.0, 1.0);

  // initialization
  createCollisionObjectHelper<siconos::collision::SiconosCapsule, btCapsuleShape,
                              siconos::collision::RigidBodyDS, BodyCapsuleRecord>(
      base, ds, capsule, btcapsule, bodyShapeMap, contactor, staticBodyShapeMap, staticBody);
}

void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::updateShape(
    std::shared_ptr<BodyCapsuleRecord> record) {
  auto capsule = record->shape;
  auto btcapsule = record->btshape;

  // Update shape parameters
  if (capsule->version() != record->shape_version) {
    // Bullet capsule has an inside margin, so we add the outside
    // margin explicitly.
    double m = capsule->outsideMargin();

    double radius = (capsule->radius() + m) * _options->worldScale;
    double length = (capsule->length() / 2 + m) * _options->worldScale;

    assert(radius > 0 && length > 0);

    btcapsule->setLocalScaling(btVector3(radius, length, radius));
    btcapsule->setMargin((capsule->insideMargin() + capsule->outsideMargin()) *
                         _options->worldScale);
    auto rbds = std::static_pointer_cast<RigidBodyDS>(record->ds);
    if (record->ds && rbds->useContactorInertia()) updateContactorInertia(rbds, btcapsule);

    if (record->btobject->getBroadphaseHandle()) {
      _collisionWorld->updateSingleAabb(&*record->btobject);
      //      _collisionWorld->getBroadphase()->getOverlappingPairCache()->
      //        cleanProxyFromPairs(record->btobject->getBroadphaseHandle(),
      //        &*_dispatcher);
    }

    record->shape_version = capsule->version();
  }

  updateShapePosition(record);
}

void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::
    createCollisionObject(const std::shared_ptr<siconos::algebra::SiconosVector> base,
                          const std::shared_ptr<siconos::collision::RigidBodyDS> ds,
                          std::shared_ptr<siconos::collision::SiconosConvexHull> ch,
                          std::shared_ptr<siconos::collision::SiconosContactor> contactor,
                          const std::shared_ptr<siconos::collision::StaticBody> staticBody) {
  if (!ch->vertices()) THROW_EXCEPTION("No vertices matrix specified for convex hull.");

  if (ch->vertices()->cols() != 3)
    THROW_EXCEPTION("Convex hull vertices matrix must have 3 columns.");

  // Copy and scale the points
  int rows = ch->vertices()->rows();
  std::vector<btScalar> pts;
  pts.resize(rows * 3);
  for (int r = 0; r < rows; r++) {
    pts[r * 3 + 0] = (*ch->vertices())(r, 0) * _options->worldScale;
    pts[r * 3 + 1] = (*ch->vertices())(r, 1) * _options->worldScale;
    pts[r * 3 + 2] = (*ch->vertices())(r, 2) * _options->worldScale;
  }

  std::shared_ptr<btConvexHullShape> btch{nullptr};
  if (ch->insideMargin() == 0) {
    // Create a convex hull directly with no further processing.
    // TODO: In case of worldScale=1, maybe we could avoid the copy to pts.
    btch = std::make_shared<btConvexHullShape>(&pts[0], rows, sizeof(btScalar) * 3);
  } else {
    // Internal margin implemented by shrinking the hull
    // TODO: Do we need the shrink clamp? (last parameter)
    btConvexHullComputer shrinkCH;
    btScalar shrunkBy = shrinkCH.compute(&pts[0], sizeof(btScalar) * 3, rows,
                                         ch->insideMargin() * _options->worldScale, 0);
    if (shrunkBy < 0) {
      // TODO: Warning
      // "insideMargin is too large, convex hull would be too small.";
      btch = std::make_shared<btConvexHullShape>(&pts[0], rows, sizeof(btScalar) * 3);
      ch->setInsideMargin(0);
    } else {
      btch = std::make_shared<btConvexHullShape>();
      for (int i = 0; i < shrinkCH.vertices.size(); i++) {
        const btVector3& v(shrinkCH.vertices[i]);
#if defined(BT_BULLET_VERSION) && (BT_BULLET_VERSION <= 281)
        btch->addPoint(v);
#else
        btch->addPoint(v, false);
#endif
      }
      ch->setInsideMargin(shrunkBy / _options->worldScale);
    }
  }

  // Add external margin and recalc bounding box
  btch->setMargin((ch->insideMargin() + ch->outsideMargin()) * _options->worldScale);
  btch->recalcLocalAabb();

  // initialization
  createCollisionObjectHelper<siconos::collision::SiconosConvexHull, btConvexHullShape,
                              siconos::collision::RigidBodyDS, BodyCHRecord>(
      base, ds, ch, btch, bodyShapeMap, contactor, staticBodyShapeMap, staticBody);
}

void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::updateShape(
    std::shared_ptr<BodyCHRecord> record) {
  auto ch = record->shape;
  auto btch = record->btshape;

  // Update shape parameters
  if (ch->version() != record->shape_version) {
    // TODO
    // btbox->setLocalScaling(btVector3(sx, sy, sz));
    btch->setMargin((ch->insideMargin() + ch->outsideMargin()) * _options->worldScale);
    auto rbds = std::static_pointer_cast<RigidBodyDS>(record->ds);
    if (record->ds && rbds->useContactorInertia()) updateContactorInertia(rbds, btch);

    if (record->btobject->getBroadphaseHandle()) {
      _collisionWorld->updateSingleAabb(&*record->btobject);
      // _collisionWorld->getBroadphase()->getOverlappingPairCache()->
      // cleanProxyFromPairs(record->btobject->getBroadphaseHandle(),
      // &*_dispatcher);
    }

    record->shape_version = ch->version();
  }

  updateShapePosition(record);
}

void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::
    createCollisionObject(const std::shared_ptr<siconos::algebra::SiconosVector> base,
                          const std::shared_ptr<siconos::collision::RigidBodyDS> ds,
                          std::shared_ptr<siconos::collision::SiconosMesh> mesh,
                          std::shared_ptr<siconos::collision::SiconosContactor> contactor,
                          const std::shared_ptr<siconos::collision::StaticBody> staticBody) {
  if (!mesh->indexes()) THROW_EXCEPTION("No indexes matrix specified for mesh.");

  if ((mesh->indexes()->size() % 3) != 0)
    THROW_EXCEPTION("Mesh indexes size must be divisible by 3.");

  if (!mesh->vertices()) THROW_EXCEPTION("No vertices matrix specified for mesh.");

  if (mesh->vertices()->rows() != 3)
    THROW_EXCEPTION("Convex hull vertices matrix must have 3 columns.");

  // Create Bullet triangle list, either by copying on non-copying method
  // TODO: worldScale on vertices
  std::pair<std::shared_ptr<btTriangleIndexVertexArray>, btScalar*> datapair{
      make_bt_vertex_array(mesh, (btScalar)0, (*mesh->vertices())(0, 0))};
  std::shared_ptr<btTriangleIndexVertexArray> bttris(datapair.first);

  // Create Bullet mesh object
  auto btmesh =
      std::make_shared<siconos::collision::bullet::internal::SiconosMeshData>(&*bttris);

  // Hold on to the data since Bullet does not make a copy
  btmesh->btTriData = bttris;
  btmesh->btScalarVertices = datapair.second;

  // Initial bound update for btGImpaceMeshShape
  btmesh->updateBound();

  // initialization
  createCollisionObjectHelper<siconos::collision::SiconosMesh, internal::SiconosMeshData,
                              siconos::collision::RigidBodyDS, BodyMeshRecord>(
      base, ds, mesh, btmesh, bodyShapeMap, contactor, staticBodyShapeMap, staticBody);
}

void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::updateShape(
    std::shared_ptr<BodyMeshRecord> record) {
  auto mesh = record->shape;
  auto btmesh = record->btshape;

  // Update shape parameters
  if (mesh->version() != record->shape_version) {
    // btBvhTriangleMeshShape supports only outsideMargin.
    // TODO: support insideMargin, scale the points by their normals.
    btmesh->setMargin(mesh->outsideMargin() * _options->worldScale);
    btmesh->postUpdate();

    // TODO: Calculate inertia from a mesh.  For now we leave it at
    // identity, the user can provide an inertia if desired.

    // if (record->ds && record->ds->useContactorInertia())
    //   updateContactorInertia(record->ds, btmesh);

    if (record->btobject->getBroadphaseHandle()) {
      _collisionWorld->updateSingleAabb(&*record->btobject);
      // _collisionWorld->getBroadphase()->getOverlappingPairCache()->
      // cleanProxyFromPairs(record->btobject->getBroadphaseHandle(),
      // &*_dispatcher);
    }

    record->shape_version = mesh->version();
  }

  updateShapePosition(record);
}

void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::
    createCollisionObject(const std::shared_ptr<siconos::algebra::SiconosVector> base,
                          const std::shared_ptr<siconos::collision::RigidBodyDS> ds,
                          std::shared_ptr<siconos::collision::SiconosHeightMap> heightmap,
                          std::shared_ptr<siconos::collision::SiconosContactor> contactor,
                          const std::shared_ptr<siconos::collision::StaticBody> staticBody) {
  if (!heightmap->height_data()) THROW_EXCEPTION("No height matrix specified for heightmap.");

  auto data = heightmap->height_data();

  if (!data || data->rows() < 2 || data->cols() < 2)
    THROW_EXCEPTION(
        "Height matrix does not have sufficient dimensions "
        "to represent a plane.");

  // Create heightfield data for Bullet.  Make a copy in case data
  // type btScalar is different from SiconosMatrix.
  // Calculate min and max value at the same time.
  double vmin = std::numeric_limits<double>::infinity();
  double vmax = -vmin;
  std::shared_ptr<std::vector<btScalar>> heightfield(
      std::make_shared<std::vector<btScalar>>());

  heightfield->resize(data->rows() * data->cols());

  for (siconos::algebra::Index i = 0; i < data->rows(); i++) {
    for (siconos::algebra::Index j = 0; j < data->cols(); j++) {
      double v = (*data)(i, j);
      (*heightfield)[j * data->rows() + i] = v;
      if (v > vmax) vmax = v;
      if (v < vmin) vmin = v;
    }
  }

  // Create Bullet height object
  auto btheight = std::make_shared<siconos::collision::bullet::internal::SiconosHeightData>(
      data->rows(), heightfield, vmin, vmax);

  // initialization
  auto btobject =
      createCollisionObjectHelper<siconos::collision::SiconosHeightMap,
                                  siconos::collision::bullet::internal::SiconosHeightData,
                                  siconos::collision::RigidBodyDS, BodyHeightRecord>(
          base, ds, heightmap, btheight, bodyShapeMap, contactor, staticBodyShapeMap,
          staticBody);

  // this flag allows to call the call gContactAddedCallback when the callback
  // has just been in the manifold In the case of the heightmap, we use it to
  // tweak the normal to avoid internal edge contact.
  btobject->setCollisionFlags(btCollisionObject::CF_CUSTOM_MATERIAL_CALLBACK);
}

void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::updateShape(
    std::shared_ptr<BodyHeightRecord> record) {
  auto height = record->shape;
  auto btheight = record->btshape;

  // Update shape parameters
  if (height->version() != record->shape_version) {
    // btBvhTriangleHeightShape supports only outsideMargin.
    // TODO: support insideMargin, scale the points by their normals.
    btheight->setMargin((height->insideMargin() + height->outsideMargin()) *
                        _options->worldScale);

    // The local scaling determines the extents of the base of the heightmap
    btheight->setLocalScaling(
        btVector3(height->length_x() / (height->height_data()->rows() - 1),
                  height->length_y() / (height->height_data()->cols() - 1), 1));

    // TODO vertical position offset to compensate for Bullet's centering
    //  TODO: Calculate the local Aabb
    // btheight->recalcLocalAabb();

    // TODO: Calculate inertia from a height.  For now we leave it at
    // identity, the user can provide an inertia if desired.

    // if (record->ds && record->ds->useContactorInertia())
    //   updateContactorInertia(record->ds, btheight);

    if (record->btobject->getBroadphaseHandle()) {
      _collisionWorld->updateSingleAabb(&*record->btobject);
      // _collisionWorld->getBroadphase()->getOverlappingPairCache()->
      // cleanProxyFromPairs(record->btobject->getBroadphaseHandle(),
      // &*_dispatcher);
    }

    record->shape_version = height->version();
  }

  // Like updateShapePosition(record), but we have to pre-compensate
  // the Bullet vertical centering of btHeightfieldTerrainShape.  Must
  // be done before body rotation.
  // Bullet automatically moves the center of the object to the
  // vertical center of the heightfield, so combine it here with the
  // contactor offset.
  auto mnz = btheight->_min_height, mxz = btheight->_max_height;
  auto z_offset = (mxz - mnz) / 2 + mnz - height->insideMargin();
  siconos::algebra::SiconosVector o(7);
  o.setZero();
  o(2) = z_offset;
  o(3) = 1;

  auto t = offsetTransform(record->contactor->offset, o);

  o(0) = t.getOrigin().getX();
  o(1) = t.getOrigin().getY();
  o(2) = t.getOrigin().getZ();
  o(3) = t.getRotation().getW();
  o(4) = t.getRotation().getX();
  o(5) = t.getRotation().getY();
  o(6) = t.getRotation().getZ();

  // Now apply the combined height and contactor offset o to the base
  // transform q
  siconos::algebra::SiconosVector q(7);
  if (record->base)
    q = *record->base;
  else {
    q.setZero();
    q(3) = 1;
  }
  t = offsetTransform(q, o);

  // DEBUG_PRINTF("updating shape position: %p(%ld) - %f, %f, %f\n",
  //              &*box,box.use_count(), q(0), q(1), q(2));

  t.setOrigin(t.getOrigin() * _options->worldScale);
  record->btobject->setWorldTransform(t);
}

///////////////////////////////////////////////////////////////////////////
// 2D shapes
///////////////////////////////////////////////////////////////////////////

void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::
    createCollisionObject(
        const std::shared_ptr<siconos::algebra::SiconosVector> base,
        const std::shared_ptr<siconos::collision::RigidBody2dDS> ds,
        std::shared_ptr<siconos::collision::SiconosDisk> disk,
        const std::shared_ptr<siconos::collision::SiconosContactor> contactor,
        const std::shared_ptr<siconos::collision::StaticBody> staticBody) {
  DEBUG_BEGIN(
      "void "
      "siconos::collision::bullet::internal::SiconosBulletCollisionManager_"
      "impl::"
      "createCollisionObject(..., disk, ...)\n");
  // set radius to 1.0 and use scaling instead of setting radius
  // directly, makes it easier to change during update

  // This version is ok
  double SCALING = 1.0;

  auto* childShape2 = new btCylinderShapeZ(
      btVector3(btScalar(SCALING * 1), btScalar(SCALING * 1), btScalar(_options->Depth2D)));
  // btConvexShape* colShape3= new btConvex2dShape(childShape2);
  auto btconvex2d1 = std::make_shared<btConvex2dShape>(childShape2);

  // //This version not
  // auto btcylinder =
  // std::make_shared<btCylinderShapeZ(btVector3(1.0, 1.0, 1.0))); auto
  // btconvex2d=std::make_shared< btConvex2dShape(btcylinder.get()));

  // //This version not
  // btConvexShape* btcylinder = new btCylinderShapeZ(btVector3(1.0, 1.0,
  // 0.04)); auto btconvex2d= std::make_shared< btConvex2dShape(btcylinder));
  // btcylinder->setMargin(0.0);

  // initialization
  createCollisionObjectHelper<siconos::collision::SiconosDisk, btConvex2dShape,
                              siconos::collision::RigidBody2dDS, BodyDiskRecord>(
      base, ds, disk, btconvex2d1, bodyShapeMap, contactor, staticBodyShapeMap, staticBody);
  DEBUG_END(
      "void "
      "siconos::collision::bullet::internal::SiconosBulletCollisionManager_"
      "impl::"
      "createCollisionObject(..., disk, ..) \n");
}

void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::updateShape(
    std::shared_ptr<BodyDiskRecord> record) {
  DEBUG_BEGIN(
      "siconos::collision::bullet::internal::SiconosBulletCollisionManager_"
      "impl::"
      "updateShape("
      "BodyDiskRecord &record)\n");
  auto disk = record->shape;
  auto btconvex2d = record->btshape;

  // Update shape parameters
  if (disk->version() != record->shape_version) {
    // Bullet cylinder has an inside margin, so we add the outside
    // margin explicitly.
    double m = disk->outsideMargin();

    double radius = (disk->radius() + m) * _options->worldScale;

    assert(radius > 0);
    DEBUG_PRINTF("outside margin=%f \n", m);
    DEBUG_PRINTF("radius=%f \n", radius);
    DEBUG_PRINTF("_options->worldScale=%f \n", _options->worldScale);
    btconvex2d->setLocalScaling(btVector3(radius, radius, radius / 25.0));
    btconvex2d->setMargin((disk->insideMargin() + disk->outsideMargin()) *
                          _options->worldScale);

    auto rbds = std::static_pointer_cast<RigidBody2dDS>(record->ds);
    if (record->ds && rbds->useContactorInertia()) update2DContactorInertia(rbds, btconvex2d);

    if (record->btobject->getBroadphaseHandle()) {
      _collisionWorld->updateSingleAabb(&*record->btobject);
      // _collisionWorld->getBroadphase()->getOverlappingPairCache()->
      // cleanProxyFromPairs(record->btobject->getBroadphaseHandle(),
      // &*_dispatcher);
    }

    record->shape_version = disk->version();
  }

  updateShapePosition(record);
  DEBUG_END(
      "siconos::collision::bullet::internal::SiconosBulletCollisionManager_"
      "impl::"
      "updateShape("
      "BodyDiskRecord &record)\n");
}

void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::
    createCollisionObject(
        const std::shared_ptr<siconos::algebra::SiconosVector> base,
        const std::shared_ptr<siconos::collision::RigidBody2dDS> ds,
        std::shared_ptr<siconos::collision::SiconosBox2d> box2d,
        const std::shared_ptr<siconos::collision::SiconosContactor> contactor,
        const std::shared_ptr<siconos::collision::StaticBody> staticBody) {
  DEBUG_BEGIN(
      "void "
      "siconos::collision::bullet::internal::SiconosBulletCollisionManager_"
      "impl::"
      "createCollisionObject(..., box2d, ...)\n");
  // set radius to 1.0 and use scaling instead of setting radius
  // directly, makes it easier to change during update

  // This version is ok
  double SCALING = 1.0;
  btConvexShape* childShape0 = new btBoxShape(
      btVector3(btScalar(SCALING * 1), btScalar(SCALING * 1), btScalar(SCALING * 1)));
  // btConvexShape* colShape= new btConvex2dShape(childShape0);
  auto btconvex2d = std::make_shared<btConvex2dShape>(childShape0);

  // initialization
  createCollisionObjectHelper<siconos::collision::SiconosBox2d, btConvex2dShape,
                              siconos::collision::RigidBody2dDS, BodyBox2dRecord>(
      base, ds, box2d, btconvex2d, bodyShapeMap, contactor, staticBodyShapeMap, staticBody);
  DEBUG_END(
      "void  "
      "siconos::collision::bullet::internal::SiconosBulletCollisionManager_"
      "impl::"
      "createCollisionObject(..., box2d, ..) \n");
}

void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::updateShape(
    std::shared_ptr<BodyBox2dRecord> record) {
  DEBUG_BEGIN(
      "siconos::collision::bullet::internal::SiconosBulletCollisionManager_"
      "impl::"
      "updateShape("
      "BodyBox2dRecord &record)\n");
  auto box2d = record->shape;
  auto btconvex2d = record->btshape;

  // Update shape parameters
  if (box2d->version() != record->shape_version) {
    // Bullet cylinder has an inside margin, so we add the outside
    // margin explicitly.
    double m = box2d->outsideMargin();

    auto dimensions = box2d->dimensions();

    double width = ((*dimensions)(0) + m) * _options->worldScale;

    double height = ((*dimensions)(1) + m) * _options->worldScale;

    assert(width > 0);
    assert(height > 0);
    DEBUG_PRINTF("outside margin=%f \n", m);
    DEBUG_PRINTF("width=%f \n", width);
    DEBUG_PRINTF("height=%f \n", height);

    DEBUG_PRINTF("_options->worldScale=%f \n", _options->worldScale);
    btconvex2d->setLocalScaling(
        btVector3(width / 2.0, height / 2.0, _options->Depth2D * _options->worldScale / 2.0));
    btconvex2d->setMargin((box2d->insideMargin() + box2d->outsideMargin()) *
                          _options->worldScale);

    auto rbds = std::static_pointer_cast<RigidBody2dDS>(record->ds);
    if (record->ds && rbds->useContactorInertia()) update2DContactorInertia(rbds, btconvex2d);

    if (record->btobject->getBroadphaseHandle()) {
      _collisionWorld->updateSingleAabb(&*record->btobject);
      // _collisionWorld->getBroadphase()->getOverlappingPairCache()->
      // cleanProxyFromPairs(record->btobject->getBroadphaseHandle(),
      // &*_dispatcher);
    }

    record->shape_version = box2d->version();
  }

  updateShapePosition(record);
  DEBUG_END(
      "siconos::collision::bullet::internal::SiconosBulletCollisionManager_"
      "impl::"
      "updateShape("
      "BodyBox2dRecord &record)\n");
}

void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::
    createCollisionObject(
        const std::shared_ptr<siconos::algebra::SiconosVector> base,
        const std::shared_ptr<siconos::collision::RigidBody2dDS> ds,
        std::shared_ptr<siconos::collision::SiconosConvexHull2d> ch2d,
        const std::shared_ptr<siconos::collision::SiconosContactor> contactor,
        const std::shared_ptr<siconos::collision::StaticBody> staticBody) {
  DEBUG_BEGIN(
      "void  "
      "siconos::collision::bullet::internal::SiconosBulletCollisionManager_"
      "impl::"
      "createCollisionObject(..., ch2d, ...)\n");
  // set radius to 1.0 and use scaling instead of setting radius
  // directly, makes it easier to change during update
  if (!ch2d->vertices()) THROW_EXCEPTION("No vertices matrix specified for convex hull.");

  if (ch2d->vertices()->cols() != 2)
    THROW_EXCEPTION("2d Convex hull vertices matrix must have 2 columns.");

  // First way. We avoid to double the point
  // This works well if the  _options->worldScale is near to 1.
  // for a unknown reason
  // int rows = ch2d->vertices()->rows();
  // std::vector<btScalar> pts;
  // pts.resize(rows*3);
  // for(int r=0; r < rows; r++)
  // {
  //   pts[r*3+0] = (*ch2d->vertices())(r, 0) * _options->worldScale;
  //   pts[r*3+1] = (*ch2d->vertices())(r, 1) * _options->worldScale;
  //   pts[r*3+2] = 0.0;
  // }

  // Second way. We double the points
  // it seems to be more robust for the contact detection
  // it avoids to find contact on the edge of the convex hull "plate"
  // Copy and scale the points
  int rows2d = ch2d->vertices()->rows();
  int rows = rows2d * 2;

  std::vector<btScalar> pts;
  pts.resize(rows * 3);
  for (int r = 0; r < rows2d; r++) {
    pts[r * 3 + 0] = (*ch2d->vertices())(r, 0) * _options->worldScale;
    pts[r * 3 + 1] = (*ch2d->vertices())(r, 1) * _options->worldScale;
    pts[r * 3 + 2] = _options->Depth2D * _options->worldScale / 2.0;
  }
  for (int r = rows2d; r < rows; r++) {
    pts[r * 3 + 0] = (*ch2d->vertices())(r - rows2d, 0) * _options->worldScale;
    pts[r * 3 + 1] = (*ch2d->vertices())(r - rows2d, 1) * _options->worldScale;
    pts[r * 3 + 2] = -_options->Depth2D * _options->worldScale / 2.0;
  }

  DEBUG_EXPR_WE(for (int r = 0; r < rows; r++)
                    printf("pts[r*3+0] = %8.4e, pts[r*3+1] =%8.4e, pts[r*3+2] =%8.4e\n",
                           pts[r * 3 + 0], pts[r * 3 + 1], pts[r * 3 + 2]););

  // This version is ok
  // btConvexHullShape* childShape1 = new btConvexHullShape(&pts[0],rows,
  // sizeof(btScalar)*3);

  auto* btch =
      new btConvexHullShape(&pts[0], rows,
                            sizeof(btScalar) * 3);  // Warning: Possible loss of memory

  btScalar shrunkBy = 0.0;

  if (ch2d->insideMargin() == 0) {
    shrunkBy = -1.0;
  } else {
    // Internal margin implemented by shrinking the hull
    // TODO: Do we need the shrink clamp? (last parameter)
    // If "shrinkClamp" is positive, "shrink" is clamped to not exceed
    // "shrinkClamp * innerRadius", where "innerRadius"
    //  is the minimum distance of a face to the center of the convex hull.

    // btConvexHullComputer shrinkCH_0;
    // DEBUG_PRINTF("Internal margin implemented by shrinking the hull of %e\n",
    // 0.0); btScalar shrunkBy_0 = shrinkCH_0.compute(&pts[0],
    // sizeof(btScalar)*3, rows,
    //                                          0.0,
    //                                          0);
    // for(int i=0; i < shrinkCH_0.vertices.size(); i++)
    // {
    //   printf("shrinkCH_0.original_vertex_index[%i] %i \n", i,
    //   shrinkCH_0.original_vertex_index[i] );
    // }

    btConvexHullComputer shrinkCH;
    DEBUG_PRINTF("Internal margin implemented by shrinking the hull of %e\n",
                 ch2d->insideMargin() * _options->worldScale);
    shrunkBy = shrinkCH.compute(&pts[0], sizeof(btScalar) * 3, rows,
                                ch2d->insideMargin() * _options->worldScale, 0);
    if (shrunkBy < 0) {
      // TODO: Warning
      // "insideMargin is too large, convex hull would be too small.";
      DEBUG_PRINTF(
          "insideMargin is too large, convex hull would be too small. shrunkby "
          "%e\n ",
          shrunkBy);
      DEBUG_PRINT("come back to original convex hull\n");

      // btch = new btConvexHullShape(&pts[0], rows, sizeof(btScalar)*3);  //
      // Warning: Possible loss of memory since we cannot SP here
      ch2d->setInsideMargin(0.);
    } else {
      delete (btch);
      btch = new btConvexHullShape;
      for (int i = 0; i < shrinkCH.vertices.size(); i++) {
        // printf("shrinkCH.original_vertex_index[%i] %i \n", i,
        // shrinkCH_0.original_vertex_index[i] );
        // printf("shrinkCH.original_vertex_index[%i] %i \n", i,
        // shrinkCH.original_vertex_index[i] ); const btVector3
        // &v(shrinkCH.vertices[shrinkCH_0.original_vertex_index[i] ]);
        const btVector3& v(shrinkCH.vertices[i]);
#if defined(BT_BULLET_VERSION) && (BT_BULLET_VERSION <= 281)
        btch->addPoint(v);
#else
        btch->addPoint(v, false);
#endif
      }
      DEBUG_PRINTF("shrinking by : %e\n", shrunkBy / _options->worldScale);
      ch2d->setInsideMargin(shrunkBy / _options->worldScale);
    }
  }

  DEBUG_EXPR(display_info_btConvexHullShape(*btch););

  // recalc bounding box
  btch->recalcLocalAabb();

  auto btconvex2d = std::make_shared<btConvex2dShape>(btch);

  // Add external margin
  //  set the margin of the  btConvex2dShape. This will set the margin of the
  //  child shape btConvexHullShape
  DEBUG_PRINTF("ch2d->insideMargin() = %8.4e\t, ch2d->outsideMargin() = %8.4e\n",
               ch2d->insideMargin(), ch2d->outsideMargin());
  btconvex2d->setMargin((ch2d->insideMargin() + ch2d->outsideMargin()) * _options->worldScale);

  // initialization
  auto btobject =
      createCollisionObjectHelper<siconos::collision::SiconosConvexHull2d, btConvex2dShape,
                                  siconos::collision::RigidBody2dDS, BodyCH2dRecord>(
          base, ds, ch2d, btconvex2d, bodyShapeMap, contactor, staticBodyShapeMap, staticBody);

  if (ch2d->avoidInternalEdgeContact()) {
    // this flag allows to call the call gContactAddedCallback when the callback
    // has just been in the manifold In the case of the heightmap, we use it to
    // tweak the normal to avoid internal edge contact.
    btobject->setCollisionFlags(btCollisionObject::CF_CUSTOM_MATERIAL_CALLBACK);

    // keep track of the points of the selected edge
    // UGLY. the only way we found to keep track of the edge since the shinkring
    // method re-sorts the vertices.
    if (shrunkBy > 0) {
      // find the closest point to original point of index _normal_edge_pointA
      int r = ch2d->_normal_edge_pointA;
      btVector3 pointA{pts[r * 3 + 0], pts[r * 3 + 1], pts[r * 3 + 2]};
      ch2d->_normal_edge_pointA = find_index_closest_point_btConvexHullShape(pointA, *btch);
      // find the closest point to original point of index _normal_edge_pointB
      r = ch2d->_normal_edge_pointB;
      btVector3 pointB{pts[r * 3 + 0], pts[r * 3 + 1], pts[r * 3 + 2]};
      ch2d->_normal_edge_pointB = find_index_closest_point_btConvexHullShape(pointB, *btch);
    }
  }

  DEBUG_END(
      "void  "
      "siconos::collision::bullet::internal::SiconosBulletCollisionManager_"
      "impl::"
      "createCollisionObject(..., ch2d, ..) \n");
}

void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::updateShape(
    std::shared_ptr<BodyCH2dRecord> record) {
  DEBUG_BEGIN(
      "siconos::collision::bullet::internal::SiconosBulletCollisionManager_"
      "impl::"
      "updateShape("
      "BodyCH2dRecord &record)\n");
  auto ch2d = record->shape;
  auto btconvex2d = record->btshape;

  // Update shape parameters
  if (ch2d->version() != record->shape_version) {
    DEBUG_PRINT("We update the shape parameters");
    // TODO
    // btbox->setLocalScaling(btVector3(sx, sy, sz));
    btconvex2d->setMargin((ch2d->insideMargin() + ch2d->outsideMargin()) *
                          _options->worldScale);
    // btconvex2d->setMargin((ch2d->outsideMargin()) * _options->worldScale);

    auto rbds = std::static_pointer_cast<RigidBody2dDS>(record->ds);
    if (record->ds && rbds->useContactorInertia()) {
      DEBUG_PRINT("We update the inertia using the contactor inertia");
      update2DContactorInertia(rbds, btconvex2d);
    }

    if (record->btobject->getBroadphaseHandle()) {
      _collisionWorld->updateSingleAabb(&*record->btobject);
      // _collisionWorld->getBroadphase()->getOverlappingPairCache()->
      // cleanProxyFromPairs(record->btobject->getBroadphaseHandle(),
      // &*_dispatcher);
    }

    record->shape_version = ch2d->version();
  }

  updateShapePosition(record);
  DEBUG_END(
      "siconos::collision::bullet::internal::SiconosBulletCollisionManager_"
      "impl::"
      "updateShape("
      "BodyCH2dRecord &record)\n");
}

void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::
    updateShapePosition(std::shared_ptr<BodyBulletShapeRecord> record) {
  DEBUG_BEGIN(
      "siconos::collision::bullet::internal::SiconosBulletCollisionManager_"
      "impl::"
      "updateShapePosition(...)\n");
  siconos::algebra::SiconosVector q(7);
  if (record->base) {
    DEBUG_EXPR(record->base->display(););
    if (record->base->size() == 7) {
      DEBUG_PRINT("3D DS\n");
      if (_options->extrapolationCoefficient) {
        if (record->ds) {
          // record.ds->display();
          auto rbds = std::dynamic_pointer_cast<RigidBodyDS>(record->ds);
          rbds->compute_extrapolated_position(_options->extrapolationCoefficient);
        }

        q = *record->base;
      } else
        q = *record->base;

    } else if (record->base->size() == 3) {
      DEBUG_PRINT("2D DS\n");
      q(0) = (*record->base)(0);
      q(1) = (*record->base)(1);
      q(2) = 0.0;
      double angle = (*record->base)(2);
      q(3) = cos(angle / 2.);
      q(4) = 0.0;
      q(5) = 0.0;
      q(6) = sin(angle / 2.);
    }
  } else {
    q.setZero();
    q(3) = 1;
  }
  DEBUG_PRINT("Position of the shape given to bullet:")
  DEBUG_EXPR_WE(siconos::algebra::print(q););

  auto t = offsetTransform(q, record->contactor->offset);

  t.setOrigin(t.getOrigin() * _options->worldScale);
  DEBUG_PRINTF("transformation = %f,%f,%f\n", float(t.getOrigin().getX()),
               float(t.getOrigin().getY()), float(t.getOrigin().getZ()));
  DEBUG_PRINTF("Rotation axis = %f,%f,%f\n", float(t.getRotation().getAxis().getX()),
               float(t.getRotation().getAxis().getY()),
               float(t.getRotation().getAxis().getZ()));
  record->btobject->setWorldTransform(t);
  DEBUG_END(
      "siconos::collision::bullet::internal::SiconosBulletCollisionManager_"
      "impl::"
      "updateShapePosition(...)\n");
}

void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::
    createCollisionObjectsForStaticBodyContactorSet(
        const std::shared_ptr<StaticBody> staticBody,
        std::shared_ptr<SiconosContactorSet> contactors) {
  DEBUG_BEGIN(
      "siconos::collision::bullet::internal::SiconosBulletCollisionManager_"
      "impl::"
      "createCollisionObjectsForStaticBodyContactorSet(...)\n");
  assert(contactors);

  std::shared_ptr<SiconosContactorSet> con{contactors};
  if (!con) {
    DEBUG_PRINT("No contactors");
    DEBUG_BEGIN(
        "siconos::collision::bullet::internal::SiconosBulletCollisionManager_"
        "impl::"
        "createCollisionObjectsForStaticBodyContactorSet(...)"
        "\n");
    return;
  }
  auto ccosv = std::make_shared<CreateCollisionObjectShapeVisitor>(
      this->shared_from_this(), nullptr, staticBody->base, staticBody);

  /* Call createCollisionObject for each shape type using the visitor
   * defined above */
  for (auto& it : *con->vector()) {
    // special collision group -1 = do not collide, thus we can skip
    // creation of associated collision objects
    if (it->collision_group == -1) continue;

    // otherwise visit the object with createCollisionObject
    ccosv->contactor = it;
    ccosv->contactor->shape->accept(ccosv);
  }
  DEBUG_END(
      "siconos::collision::bullet::internal::SiconosBulletCollisionManager_"
      "impl::"
      "createCollisionObjectsForStaticBodyContactorSet(...)\n");
}

void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::
    createCollisionObjectsForBodyContactorSetFromDS(
        const std::shared_ptr<siconos::modeling::SecondOrderDS> ds) {
  DEBUG_BEGIN(
      "siconos::collision::bullet::internal::SiconosBulletCollisionManager_"
      "impl::"
      "createCollisionObjectsForBodyContactorSetFromDS(...)\n");
  assert(ds);

  auto rbds = std::dynamic_pointer_cast<RigidBodyDS>(ds);
  auto rb2dds = std::dynamic_pointer_cast<RigidBody2dDS>(ds);
  std::shared_ptr<SiconosContactorSet> con{nullptr};

  std::shared_ptr<siconos::algebra::SiconosVector> base;

  if (rbds) {
    DEBUG_PRINT("RigidBodyDS case");
    con = rbds->contactors();
    // printf("\n\n extrapolation %e\n", _options.extrapolationCoefficient);
    // getchar();
    if (_options->extrapolationCoefficient > 0.) {
      rbds->compute_extrapolated_position(_options->extrapolationCoefficient);
      base = rbds->base_extrapolated_position();
    } else {
      base = rbds->base_position();
      // siconos::algebra::print(*base);
    }
  }
  if (rb2dds) {
    DEBUG_PRINT("RigidBody2dDS case");
    con = rb2dds->contactors();
    base = rb2dds->q();
  }
  if (!con) {
    DEBUG_PRINT("No contactors");
    DEBUG_BEGIN(
        "siconos::collision::bullet::internal::SiconosBulletCollisionManager_"
        "impl::"
        "createCollisionObjectsForBodyContactorSetFromDS(...)"
        "\n");
    return;
  }

  auto ccosv = std::make_shared<CreateCollisionObjectShapeVisitor>(this->shared_from_this(),
                                                                   ds, base, nullptr);

  /* Call createCollisionObject for each shape type using the visitor
   * defined above */
  for (auto it : *con->vector()) {
    // special collision group -1 = do not collide, thus we can skip
    // creation of associated collision objects
    if (it->collision_group == -1) continue;

    // otherwise visit the object with createCollisionObject
    ccosv->contactor = it;
    ccosv->contactor->shape->accept(ccosv);

    DEBUG_END(
        "siconos::collision::bullet::internal::SiconosBulletCollisionManager_"
        "impl::"
        "createCollisionObjectsForBodyContactorSetFromDS(...)\n");
  }
}

btTransform
siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::offsetTransform(
    const siconos::algebra::SiconosVector& position,
    const siconos::algebra::SiconosVector& offset) {
  /* Adjust offset position according to current rotation */
  btQuaternion rbase(position(4), position(5), position(6), position(3));
  btVector3 rboffset = quatRotate(rbase, btVector3(offset(0), offset(1), offset(2)));

  /* Calculate total orientation */
  btQuaternion roffset(offset(4), offset(5), offset(6), offset(3));

  /* Set the absolute shape position */
  return btTransform(rbase * roffset,
                     btVector3(position(0), position(1), position(2)) + rboffset);
}

btTransform
siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::offsetTransform(
    const siconos::algebra::SiconosVector& position) {
  /* Adjust offset position according to current rotation */
  btQuaternion rbase(position(4), position(5), position(6), position(3));

  /* Set the absolute shape position */
  return btTransform(rbase, btVector3(position(0), position(1), position(2)));
}

void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::
    updateContactorInertia(std::shared_ptr<siconos::modeling::NewtonEulerDS> ds,
                           std::shared_ptr<btCollisionShape> btshape) {
  btVector3 localinertia;
  double scale_factor;
  scale_factor = 1. / (_options->worldScale * _options->worldScale);
  localinertia[0] = std::numeric_limits<double>::signaling_NaN();
  localinertia[1] = std::numeric_limits<double>::signaling_NaN();
  localinertia[2] = std::numeric_limits<double>::signaling_NaN();
  btshape->calculateLocalInertia(ds->scalarMass(), localinertia);

  localinertia[0] *= scale_factor;
  localinertia[1] *= scale_factor;
  localinertia[2] *= scale_factor;
  assert(!((localinertia.x() == 0.0 && localinertia.y() == 0.0 && localinertia.z() == 0.0) ||
           std::isinf(localinertia.x()) || std::isinf(localinertia.y()) ||
           std::isinf(localinertia.z())) &&
         "calculateLocalInertia() returned garbage");
  ds->setInertia(localinertia[0], localinertia[1], localinertia[2]);
}

void siconos::collision::bullet::internal::SiconosBulletCollisionManager_impl::
    update2DContactorInertia(std::shared_ptr<siconos::modeling::LagrangianDS> ds,
                             std::shared_ptr<btCollisionShape> btshape) {
  // TBD Warningx
}
