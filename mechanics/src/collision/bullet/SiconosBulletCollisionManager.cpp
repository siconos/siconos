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

/*! \file SiconosBulletCollisionManager.cpp

    \brief Definition of a Bullet-based interaction handler for contact
    detection.
*/

// Note, in general the "outside margin" is not implemented.  What is
// needed is a way to project the point detected on the external shell
// back to the shape surface.  This could be for example the closest
// point on the convex hull.  (For convex shapes.)

#include "SiconosBulletCollisionManager.hpp"

#include <BulletCollision/BroadphaseCollision/btAxisSweep3.h>
#include <BulletCollision/BroadphaseCollision/btDbvtBroadphase.h>
#include <BulletCollision/CollisionDispatch/btBox2dBox2dCollisionAlgorithm.h>
#include <BulletCollision/CollisionDispatch/btConvex2dConvex2dAlgorithm.h>
#include <BulletCollision/CollisionShapes/btHeightfieldTerrainShape.h>
#include <BulletCollision/CollisionShapes/btTriangleInfoMap.h>
#include <BulletCollision/Gimpact/btGImpactCollisionAlgorithm.h>
#include <BulletCollision/NarrowPhaseCollision/btMinkowskiPenetrationDepthSolver.h>

#include <algorithm>
#include <map>

#include "Bullet1DR.hpp"
#include "Bullet2d3DR.hpp"
#include "Bullet2dR.hpp"
#include "Bullet5DR.hpp"
#include "BulletR.hpp"
#include "Interaction.hpp"
#include "IterateContactPoint.hpp"
#include "NewtonEulerJointR.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include <FremondImpactFrictionNSL.hpp>
#include "NewtonImpactRollingFrictionNSL.hpp"
#include "RigidBody2dDS.hpp"
#include "RigidBodyDS.hpp"
#include "SiconosBulletCollisionManager_impl.hpp"
#include "SiconosCollisionQueryResult.hpp"
#include "SiconosContactor.hpp"
#include "SiconosPointers.hpp"
#include "SiconosVector.hpp"
#include "Simulation.hpp"
#include "SimulationGraphs.hpp"
#include "StaticBody.hpp"
#include "Topology.hpp"
#include "siconos_debug.h"
// #include "SiconosBulletDefines.h"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES

siconos::collision::bullet::SiconosBulletCollisionManager::SiconosBulletCollisionManager(
    std::shared_ptr<SiconosBulletOptions> options)
    : SiconosCollisionManager{}, _options{options} {
  assert(options);
  // if(not _options)
  //   _options = std::make_shared<SiconosBulletOptions>(); // Default options set
  initialize_impl();
}

siconos::collision::bullet::SiconosBulletCollisionManager::SiconosBulletCollisionManager()
    : SiconosBulletCollisionManager{std::make_shared<SiconosBulletOptions>()} {}

void siconos::collision::bullet::SiconosBulletCollisionManager::initialize_impl() {
  _impl = std::make_shared<internal::SiconosBulletCollisionManager_impl>(_options);

  // collision configuration contains default setup for memory, collision setup
  _impl->_collisionConfiguration.reset(new btDefaultCollisionConfiguration());

  if (_options->perturbationIterations > 0 ||
      _options->minimumPointsPerturbationThreshold > 0) {
    _impl->_collisionConfiguration->setConvexConvexMultipointIterations(
        _options->perturbationIterations, _options->minimumPointsPerturbationThreshold);
    _impl->_collisionConfiguration->setPlaneConvexMultipointIterations(
        _options->perturbationIterations, _options->minimumPointsPerturbationThreshold);
  }

  // use the default collision dispatcher. For parallel processing you can use a diffent
  // dispatcher (see Extras/BulletMultiThreaded)
  _impl->_dispatcher.reset(new btCollisionDispatcher(&*_impl->_collisionConfiguration));

  if (_options->useAxisSweep3)
    _impl->_broadphase.reset(new btAxisSweep3(btVector3(), btVector3()));
  else
    _impl->_broadphase.reset(new btDbvtBroadphase());

  _impl->_collisionWorld.reset(new btCollisionWorld(&*_impl->_dispatcher, &*_impl->_broadphase,
                                                    &*_impl->_collisionConfiguration));

  btOverlapFilterCallback *filterCallback = new internal::SiconosBulletFilterCallback();
  reinterpret_cast<internal::SiconosBulletFilterCallback *>(filterCallback)
      ->interactionManager = this;
  _impl->_collisionWorld->getPairCache()->setOverlapFilterCallback(filterCallback);

  DEBUG_PRINTF("_options->dimension = %i", _options->dimension);

  // 2D specific
  if (_options->dimension == SiconosBulletDimension::TwoD) {
    btVoronoiSimplexSolver *m_simplexSolver = new btVoronoiSimplexSolver();
    btMinkowskiPenetrationDepthSolver *m_pdSolver = new btMinkowskiPenetrationDepthSolver();

    btConvex2dConvex2dAlgorithm::CreateFunc *m_convexAlgo2d =
        new btConvex2dConvex2dAlgorithm::CreateFunc(m_simplexSolver, m_pdSolver);
    btBox2dBox2dCollisionAlgorithm::CreateFunc *m_box2dbox2dAlgo =
        new btBox2dBox2dCollisionAlgorithm::CreateFunc();

    _impl->_dispatcher->registerCollisionCreateFunc(CONVEX_2D_SHAPE_PROXYTYPE,
                                                    CONVEX_2D_SHAPE_PROXYTYPE, m_convexAlgo2d);
    _impl->_dispatcher->registerCollisionCreateFunc(BOX_2D_SHAPE_PROXYTYPE,
                                                    CONVEX_2D_SHAPE_PROXYTYPE, m_convexAlgo2d);
    _impl->_dispatcher->registerCollisionCreateFunc(CONVEX_2D_SHAPE_PROXYTYPE,
                                                    BOX_2D_SHAPE_PROXYTYPE, m_convexAlgo2d);
    _impl->_dispatcher->registerCollisionCreateFunc(BOX_2D_SHAPE_PROXYTYPE,
                                                    BOX_2D_SHAPE_PROXYTYPE, m_box2dbox2dAlgo);
  } else
    btGImpactCollisionAlgorithm::registerAlgorithm(&*_impl->_dispatcher);

  _impl->_collisionWorld->getDispatchInfo().m_useContinuous = false;
  _impl->_collisionWorld->getDispatchInfo().m_enableSatConvex = _options->enableSatConvex;
}

siconos::collision::bullet::SiconosBulletCollisionManager::
    ~SiconosBulletCollisionManager() noexcept {
  // unlink() will be called on all remaining
  // contact points when world is destroyed

  // must be the first de-allocated, otherwise segfault
  _impl->_collisionWorld.reset();
}

std::shared_ptr<siconos::collision::StaticBody>
siconos::collision::bullet::SiconosBulletCollisionManager::addStaticBody(
    std::shared_ptr<siconos::collision::SiconosContactorSet> cs,
    const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>> &position,
    int number) {
  auto rec = std::make_shared<siconos::collision::StaticBody>();
  rec->contactorSet = cs;
  if (!position) {
    rec->base = std::make_shared<siconos::algebra::SiconosVector>(7);
    rec->base->setZero();
    (*(rec->base))(3) = 1.0;  // we give a unit identity quarternion
  } else {
    rec->base = std::make_shared<siconos::algebra::SiconosVector>(*position);
  }
  rec->number = number;
  // std::cout << "siconos::collision::bullet::SiconosBulletCollisionManager::addStaticBody
  // number : " << number << std::endl;
  _impl->createCollisionObjectsForStaticBodyContactorSet(rec, cs);
  return rec;
}

void siconos::collision::bullet::SiconosBulletCollisionManager::removeStaticBody(
    const std::shared_ptr<siconos::collision::StaticBody> &body) {
  auto it = _impl->staticBodyShapeMap.find(&*body);
  if (it == _impl->staticBodyShapeMap.end()) return;

  for (auto &it2 : it->second)
    _impl->_collisionWorld->removeCollisionObject(&*(it2)->btobject);

  _impl->staticBodyShapeMap.erase(it);
}

void siconos::collision::bullet::SiconosBulletCollisionManager::removeBody(
    const std::shared_ptr<siconos::modeling::SecondOrderDS> &body) {
  auto it = _impl->bodyShapeMap.find(&*body);
  if (it == _impl->bodyShapeMap.end()) return;

  // For the given body, loop through records in the bodyShapeMap
  for (auto &records : it->second) {
    std::visit(
        [this](auto rec) { _impl->_collisionWorld->removeCollisionObject(&*(rec)->btobject); },
        records);
  }
  _impl->bodyShapeMap.erase(it);
}

// called once for each contact point as it is destroyed
siconos::simulation::Simulation
    *siconos::collision::bullet::SiconosBulletCollisionManager::gSimulation;

bool siconos::collision::bullet::SiconosBulletCollisionManager::bulletContactClear(
    void *userPersistentData) {
  /* note: stored pointer to shared_ptr! */
  std::shared_ptr<siconos::modeling::Interaction> *p_inter =
      (std::shared_ptr<siconos::modeling::Interaction> *)userPersistentData;
  assert(p_inter && "Contact point's stored (Interaction*) is null!");
  DEBUG_PRINTF("unlinking interaction %p, number %zu \n", &**p_inter, (*p_inter)->number());

  // std::shared_ptr<BulletR>
  // rel_bulletR(std::dynamic_pointer_cast<BulletR>((*p_inter)->relation())); Bullet5DR
  // rel_bullet5DR(std::dynamic_pointer_cast<Bullet5DR>((*p_inter)->relation())); Bullet2dR
  // rel_bullet2dR(std::dynamic_pointer_cast<Bullet2dR>((*p_inter)->relation()));
  // auto
  // rel_bullet2d3DR(std::dynamic_pointer_cast<Bullet2d3DR>((*p_inter)->relation())); if
  // (rel_bulletR)
  //   rel_bulletR->preDelete();
  // else if (rel_bullet5DR)
  //   rel_bullet5DR->preDelete();
  // else if (rel_bullet2dR)
  //   rel_bullet2dR->preDelete();
  // else if (rel_bullet2d3DR)
  //   rel_bullet2d3DR->preDelete();
  // std::static_pointer_cast<BulletR>((*p_inter)->relation())->preDelete();
  //_stats.interaction_destroyed++;
  gSimulation->unlink(*p_inter);
  delete p_inter;
  return false;
}

namespace {
void siconosBulletAdjustInternalEdgeContacts(btManifoldPoint &cp,
                                             const btCollisionObjectWrapper *colObj0Wrap,
                                             const btCollisionObjectWrapper *colObj1Wrap,
                                             int partId0, int index0) {
  DEBUG_BEGIN("siconosBulletAdjustInternalEdgeContacts \n");

  DEBUG_EXPR(display_info_contact_point(cp);
             display_info_collision_object(colObj0Wrap->getCollisionObject());
             display_info_collision_object(colObj1Wrap->getCollisionObject()););

  // btAssert(colObj0->getCollisionShape()->getShapeType() == TRIANGLE_SHAPE_PROXYTYPE);
  //  if (colObj0Wrap->getCollisionShape()->getShapeType() != TRIANGLE_SHAPE_PROXYTYPE)
  //  	return;
  if (colObj0Wrap->getCollisionObject()->getCollisionShape()->getShapeType() ==
      CONVEX_2D_SHAPE_PROXYTYPE) {
    // printf("CONVEX_2D_SHAPE_PROXYTYPE  : %i\n",  CONVEX_2D_SHAPE_PROXYTYPE );
    const btCollisionObject *objectA = colObj0Wrap->getCollisionObject();
    const auto *pairA =
        reinterpret_cast<const siconos::collision::bullet::BodyBulletShapeRecord *>(
            objectA->getUserPointer());
    auto sshape = pairA->sshape;

    // pairA->display();
    // const btCollisionObject*  objectB= colObj1Wrap->getCollisionObject();
    // const BodyBulletShapeRecord *pairB = reinterpret_cast<const
    // BodyBulletShapeRecord*>(objectB->getUserPointer()); pairB->display();

    auto ch2d = std::dynamic_pointer_cast<siconos::collision::SiconosConvexHull2d>(sshape);

    if (ch2d && ch2d->avoidInternalEdgeContact()) {
      DEBUG_PRINTF("a Siconos ch2d shape and ch2d->avoidInternalEdgeContact() true \n");

      // Retrieve the first two points assuming that it corresponds to the edge of interest

      btConvex2dShape *btConvex2d =
          (btConvex2dShape *)colObj0Wrap->getCollisionObject()->getCollisionShape();
      btConvexShape *btconvex = (btConvexShape *)(btConvex2d->getChildShape());
      btConvexHullShape *btch = (btConvexHullShape *)btconvex;

      // printf("number of points in convex hull : %i\n ", numPoints);
      const btVector3 *points = btch->getPoints();
      // if (numPoints  > 4)
      // {
      //   printf("Warning: number of points in convex hull is more than 2\n     We consider
      //   the two first point as the contact edge.\n"); int p=ch2d->_normal_edge_pointA;
      //   printf("   point # %i x , y, x : %e\t, %e\t, %e\t \n", p,  points[p].x(),
      //   points[p].y(), points[p].z()); p=ch2d->_normal_edge_pointB; printf("   point # %i
      //   x , y, x : %e\t, %e\t, %e\t \n", p,  points[p].x(), points[p].y(), points[p].z());
      // }

      // for (int p = 0 ; p < numPoints; p++)
      // {
      //   printf("   point # %i x , y, x : %e\t, %e\t, %e\t \n", p,  points[p].x(),
      //   points[p].y(), points[p].z());
      // }

      // Compute the normal to the selected edge
      int idx_A = ch2d->_normal_edge_pointA;
      int idx_B = ch2d->_normal_edge_pointB;

      btScalar AB_x = points[idx_B].x() - points[idx_A].x();
      btScalar AB_y = points[idx_B].y() - points[idx_A].y();
      btVector3 normal = btVector3(-AB_y, AB_x, btScalar(0.f));
      normal.safeNormalize();

      DEBUG_PRINTF(" new  normal x , y, z : %e\t, %e\t, %e\t \n", normal.x(), normal.y(),
                   normal.z());

      cp.m_normalWorldOnB = normal;

      // Reproject collision point along normal. (what about cp.m_distance1?)
      cp.m_positionWorldOnB = cp.m_positionWorldOnA - cp.m_normalWorldOnB * cp.m_distance1;
      cp.m_localPointB = colObj0Wrap->getWorldTransform().invXform(cp.m_positionWorldOnB);
    }

    DEBUG_END("siconosBulletAdjustInternalEdgeContacts \n");
    return;
    // getchar();
  }

  if (colObj0Wrap->getCollisionObject()->getCollisionShape()->getShapeType() ==
      TERRAIN_SHAPE_PROXYTYPE) {
    // btHeightfieldTerrainShape *heightfield =
    //    (btHeightfieldTerrainShape *)colObj0Wrap->getCollisionObject()->getCollisionShape();
    // auto triangleInfoMapPtr = heightfield->getTriangleInfoMap();

    btVector3 newNormal = btVector3(0, 0, 1);

    const btTriangleShape *tri_shape =
        static_cast<const btTriangleShape *>(colObj0Wrap->getCollisionShape());
    btVector3 tri_normal;
    tri_shape->calcNormal(tri_normal);

    newNormal = tri_normal;
    //					cp.m_distance1 = cp.m_distance1 *
    // newNormal.dot(cp.m_normalWorldOnB);
    // printf("old normal %e\t%e\t%e\n", oldNormal.x(),  oldNormal.y(), oldNormal.z());
    // printf("new normal %e\t%e\t%e\n", newNormal.x(),  newNormal.y(), newNormal.z());
    // printf("cp.m_distance1 = %e\n", cp.m_distance1 );

    // Option 1 - we test if the normal are similar or not before changing the normal

    // btScalar cosine =  oldNormal.dot(newNormal);
    // //printf("cosine %e\n", cosine);
    // if (cosine < 0.0)
    // {
    //   newNormal = -1.0*tri_normal;
    //   cosine =  oldNormal.dot(newNormal);
    // }
    // //btScalar diff  = oldNormal.distance(newNormal);
    // //printf("diff %e\n", diff);
    // if ((1.0 - cosine) > 3e-03 ) // around 5 degrees
    // {
    //   //printf("--------------------------------------> change edge  normal to triangle
    //   normal\n"); cp.m_normalWorldOnB = newNormal;
    // }
    // else return;

    // Option 2 - we take in any cases the normal to the triangle face

    // btScalar cosine =  oldNormal.dot(newNormal);
    // if (cosine < 0.0)
    // We assume that the normal to the triangle face must be upward.
    if (tri_normal.z() < 0.0) {
      newNormal = -1.0 * tri_normal;
    }
    cp.m_normalWorldOnB = newNormal;

    // Reproject collision point along normal. (what about cp.m_distance1?)
    cp.m_positionWorldOnB = cp.m_positionWorldOnA - cp.m_normalWorldOnB * cp.m_distance1;
    cp.m_localPointB = colObj0Wrap->getWorldTransform().invXform(cp.m_positionWorldOnB);
    DEBUG_END("siconosBulletAdjustInternalEdgeContacts \n");
    return;
  }
}
}  // namespace
bool siconos::collision::bullet::SiconosBulletCollisionManager::bulletContactAddedCallback(
    btManifoldPoint &cp, const btCollisionObjectWrapper *colObj0Wrap, int partId0, int index0,
    const btCollisionObjectWrapper *colObj1Wrap, int partId1, int index1) {
  // printf("--------- bulletContactAddedCallback start\n");
  // btAdjustInternalEdgeContacts(cp, colObj1Wrap, colObj0Wrap, partId1, index1);
  siconosBulletAdjustInternalEdgeContacts(cp, colObj1Wrap, colObj0Wrap, partId1, index1);
  // printf("--------- bulletContactAddedCallback end\n");

  return true;
}

std::shared_ptr<siconos::collision::bullet::BulletR>
siconos::collision::bullet::SiconosBulletCollisionManager::makeBulletR(
    std::shared_ptr<siconos::collision::RigidBodyDS> ds1,
    std::shared_ptr<siconos::collision::SiconosShape> shape1,
    std::shared_ptr<siconos::collision::RigidBodyDS> ds2,
    std::shared_ptr<siconos::collision::SiconosShape> shape2, const btManifoldPoint &p) {
  return std::make_shared<BulletR>();
}

std::shared_ptr<siconos::collision::bullet::Bullet5DR>
siconos::collision::bullet::SiconosBulletCollisionManager::makeBullet5DR(
    std::shared_ptr<siconos::collision::RigidBodyDS> ds1,
    std::shared_ptr<siconos::collision::SiconosShape> shape1,
    std::shared_ptr<siconos::collision::RigidBodyDS> ds2,
    std::shared_ptr<siconos::collision::SiconosShape> shape2, const btManifoldPoint &p) {
  return std::make_shared<Bullet5DR>();
}

std::shared_ptr<siconos::collision::bullet::Bullet2dR>
siconos::collision::bullet::SiconosBulletCollisionManager::makeBullet2dR(
    std::shared_ptr<siconos::collision::RigidBody2dDS> ds1,
    std::shared_ptr<siconos::collision::SiconosShape> shape1,
    std::shared_ptr<siconos::collision::RigidBody2dDS> ds2,
    std::shared_ptr<siconos::collision::SiconosShape> shape2, const btManifoldPoint &p) {
  return std::make_shared<Bullet2dR>();
}

std::shared_ptr<siconos::collision::bullet::Bullet2d3DR>
siconos::collision::bullet::SiconosBulletCollisionManager::makeBullet2d3DR(
    std::shared_ptr<siconos::collision::RigidBody2dDS> ds1,
    std::shared_ptr<siconos::collision::SiconosShape> shape1,
    std::shared_ptr<siconos::collision::RigidBody2dDS> ds2,
    std::shared_ptr<siconos::collision::SiconosShape> shape2, const btManifoldPoint &p) {
  return std::make_shared<Bullet2d3DR>();
}

void siconos::collision::bullet::SiconosBulletCollisionManager::updateInteractions(
    std::shared_ptr<siconos::simulation::Simulation> simulation) {
  DEBUG_BEGIN(
      "siconos::collision::bullet::SiconosBulletCollisionManager::updateInteractions(std::"
      "shared_ptr<siconos::simulation::"
      "Simulation> simulation)\n");

#ifdef BULLET_TIMER
  //  CProfileManager::Start_Profile("bullet_profile.txt");
  //  CProfileManager::Reset();
  auto start = std::chrono::system_clock::now();
#endif

  //  update collision objects from all RigidBodyDS dynamical systems
  auto dsg = simulation->nonSmoothDynamicalSystem()->dynamicalSystems();
  siconos::graphs::DynamicalSystemsGraph::VIterator dsi, dsiend;
  std::tie(dsi, dsiend) = dsg->vertices();

  // for (; dsi != dsiend; ++dsi) {
  //   BodiesVariant bv = dsg->bundle(*dsi);
  //   std::visit(
  //       [this](auto ds) {
  //         _impl->updateShapes(ds);
  //       },
  //       bv);
  // }

  for (; dsi != dsiend; ++dsi) {
    auto ds = dsg->bundle(*dsi);
    if (auto bds = std::dynamic_pointer_cast<RigidBodyDS>(ds)) {
      _impl->updateShapes(bds);
    } else if (auto bds = std::dynamic_pointer_cast<RigidBody2dDS>(ds)) {
      _impl->updateShapes(bds);
    } else {
      THROW_EXCEPTION("SiconosBulletManager, works only for RigidBodyDS or RigidBody2dDS.")
    }
  }

#ifdef BULLET_TIMER
  auto end = std::chrono::system_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  std::cout << "\n[mechanics] -2 : visit " << elapsed << " ms \n";
#endif
  // Clear cache automatically before collision detection if requested
  if (_options->clearOverlappingPairCache) clearOverlappingPairCache();

  if (!_impl->_queuedCollisionObjects.empty()) {
    int collisionFilterMask = 1;

    for (auto &it : _impl->_queuedCollisionObjects) {
      // std::pair<std::shared_ptr<btCollisionObject>, int> p = *it;
      int collisionFilterGroup = it.second;
      auto collisionObject = it.first;
      _impl->_collisionWorld->addCollisionObject(&*collisionObject, collisionFilterGroup,
                                                 collisionFilterMask);
    }
    _impl->_queuedCollisionObjects.clear();
  }

  // -1. reset statistical counters
  resetStatistics();

#ifdef BULLET_TIMER
  auto end_old = end;
  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - end_old).count();
  std::cout << "[mechanics] -1 : addCollisionObject " << elapsed << " ms" << std::endl;
#endif
  // 0. set up bullet callbacks
  gSimulation = &*simulation;
  gContactDestroyedCallback = this->bulletContactClear;
  gContactAddedCallback = this->bulletContactAddedCallback;

  // Important parameter controlling contact point making and breaking
  gContactBreakingThreshold = _options->contactBreakingThreshold;

  // 1. perform bullet collision detection
  _impl->_collisionWorld->performDiscreteCollisionDetection();

#ifdef BULLET_TIMER
  end_old = end;
  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - end_old).count();
  std::cout << "[mechanics]  1 : collisionDectection " << elapsed << " ms" << std::endl;
//  CProfileManager::dumpAll();
//  CProfileManager::Stop_Profile();
#endif

  DEBUG_PRINT("SiconosBulletCollisionManager :: iterating contact points:\n");
  // getchar();
  //  2. deleted contact points have been removed from the graph during the
  //     bullet collision detection callbacks

  // 3. for each contact point, if there is no interaction, create one
  internal::IterateContactPoints t{_impl->_collisionWorld};
  auto itend = t.end();
  DEBUG_EXPR_WE(

      int num_contact_points = 0;

      for (auto it : t) { num_contact_points++; }

      std::cout
      << "Number of contacts points detected by bullet: " << num_contact_points << "\n";);

  for (auto it = t.begin(); it != itend; ++it) {
    DEBUG_PRINTF("\n\n\nSiconosBulletCollisionManager ::   -- %p, %p, %p\n", it->objectA,
                 it->objectB, it->point);

    // Get the RigidBodyDS and SiconosShape pointers

    auto const *pairA =
        reinterpret_cast<const BodyBulletShapeRecord *>(it->objectA->getUserPointer());
    auto const *pairB =
        reinterpret_cast<const BodyBulletShapeRecord *>(it->objectB->getUserPointer());
    assert(pairA && pairB && "btCollisionObject had a null user pointer!");

    // The first pair will always be the non-static object
    // As a consequence, if there is a static body, it is always associated with second pair
    // pairB
    bool flip = false;
    if (pairB->ds && !pairA->ds) {
      pairA = reinterpret_cast<const BodyBulletShapeRecord *>(it->objectB->getUserPointer());
      pairB = reinterpret_cast<const BodyBulletShapeRecord *>(it->objectA->getUserPointer());
      flip = true;
    }
    DEBUG_PRINTF("SiconosBulletCollisionManager :: flip = %i \n", flip);
    // If both collision objects belong to the same body (or no body),
    // no interaction is created.
    if (pairA->ds == pairB->ds) continue;

    // If the two bodies are already connected by another type of
    // relation (e.g. EqualityCondition == they have a joint between
    // them), then don't create contact constraints, because it leads
    // to an ill-conditioned problem.

    DEBUG_EXPR_WE(
        if (pairA->ds && pairB->ds) {
          DEBUG_PRINTF("SiconosBulletCollisionManager ::   -- ds1 :  %zu,  ds2: %zu\n",
                       pairA->ds->number(), pairB->ds->number());
        } if (pairA->ds && pairB->staticBody) {
          DEBUG_PRINTF("SiconosBulletCollisionManager ::   -- ds1 :  %zu  staticbody: %i\n",
                       pairA->ds->number(), pairB->staticBody->number);
        });

    DEBUG_PRINTF("SiconosBulletCollisionManager :: _with_equality_constraints  -- %i\n",
                 _with_equality_constraints);

    if (_with_equality_constraints && pairA->ds && pairB->ds) {
      siconos::graphs::InteractionsGraph::VIterator ui, uiend;
      auto indexSet0 = simulation->nonSmoothDynamicalSystem()->topology()->indexSet0();
      bool match = false;
      for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui) {
        auto inter = indexSet0->bundle(*ui);
        auto ds1 = std::dynamic_pointer_cast<siconos::modeling::SecondOrderDS>(
            indexSet0->properties(*ui).source);
        auto ds2 = std::dynamic_pointer_cast<siconos::modeling::SecondOrderDS>(
            indexSet0->properties(*ui).target);
        if (ds1 && ds2 &&
            (((&*ds1 == &*pairA->ds) && (&*ds2 == &*pairB->ds)) ||
             ((&*ds1 == &*pairB->ds) && (&*ds2 == &*pairA->ds)))) {
          auto br = std::dynamic_pointer_cast<BulletR>(inter->relation());
          DEBUG_EXPR(std::cout << "br" << br << std::endl;);
          if (!br) {
            DEBUG_PRINT(
                "Only match on non-BulletR interactions, i.e. non-contact relations\n");
            auto jr = std::dynamic_pointer_cast<siconos::joints::NewtonEulerJointR>(
                inter->relation());

            /* If it is a joint, check the joint self-collide property */
            if (jr && !jr->allowSelfCollide()) match = true;

            /* If any non-contact relation is found, both bodies must
             * allow self-collide */
            // We need to check for other type of dynamical systems.
            auto rbdsA = std::static_pointer_cast<RigidBodyDS>(pairA->ds);
            auto rbdsB = std::static_pointer_cast<RigidBodyDS>(pairB->ds);
            if (!rbdsA->allowSelfCollide() || !rbdsB->allowSelfCollide()) match = true;
          }
          if (match) break;
        }
      }
      if (match) continue;
    }
    DEBUG_PRINTF("SiconosBulletCollisionManager :: it->point->m_userPersistentData  %p \n",
                 it->point->m_userPersistentData);
    if (it->point->m_userPersistentData) {
      /* interaction already exists */
      DEBUG_PRINT("SiconosBulletCollisionManager :: interaction already exists \n");
      std::shared_ptr<siconos::modeling::Interaction> *p_inter =
          (std::shared_ptr<siconos::modeling::Interaction> *)it->point->m_userPersistentData;

      auto rel_bulletR = std::dynamic_pointer_cast<BulletR>((*p_inter)->relation());
      auto rel_bullet5DR = std::dynamic_pointer_cast<Bullet5DR>((*p_inter)->relation());
      auto rel_bullet2dR = std::dynamic_pointer_cast<Bullet2dR>((*p_inter)->relation());
      auto rel_bullet2d3DR = std::dynamic_pointer_cast<Bullet2d3DR>((*p_inter)->relation());

      if (rel_bulletR || rel_bullet5DR) {
        DEBUG_PRINT("SiconosBulletCollisionManager :: BulletR case || rel_bullet5DR\n");
        // We need to check for other type of dynamical systems.
        auto rbdsA = std::static_pointer_cast<RigidBodyDS>(pairA->ds);
        auto rbdsB = std::static_pointer_cast<RigidBodyDS>(pairB->ds);

        /* update the relation */
        auto rel = std::static_pointer_cast<BulletR>((*p_inter)->relation());
        rel->updateContactPointsFromManifoldPoint(
            *it->manifold, *it->point, flip, _options->worldScale, rbdsA,
            rbdsB ? rbdsB : std::shared_ptr<siconos::modeling::NewtonEulerDS>());
      } else if (rel_bullet2dR) {
        DEBUG_PRINT("SiconosBulletCollisionManager :: Bullet2dR case");
        // We need to check for other type of dynamical systems.
        auto rbdsA = std::static_pointer_cast<RigidBody2dDS>(pairA->ds);
        auto rbdsB = std::static_pointer_cast<RigidBody2dDS>(pairB->ds);

        /* update the relation */
        rel_bullet2dR->updateContactPointsFromManifoldPoint(
            *it->manifold, *it->point, flip, _options->worldScale, rbdsA,
            rbdsB ? rbdsB : std::shared_ptr<siconos::collision::RigidBody2dDS>());
      } else if (rel_bullet2d3DR) {
        DEBUG_PRINT("SiconosBulletCollisionManager :: Bullet2d3DR case");
        // We need to check for other type of dynamical systems.
        auto rbdsA = std::static_pointer_cast<RigidBody2dDS>(pairA->ds);
        auto rbdsB = std::static_pointer_cast<RigidBody2dDS>(pairB->ds);

        /* update the relation */
        rel_bullet2d3DR->updateContactPointsFromManifoldPoint(
            *it->manifold, *it->point, flip, _options->worldScale, rbdsA,
            rbdsB ? rbdsB : std::shared_ptr<siconos::collision::RigidBody2dDS>());
      }

      else {
        THROW_EXCEPTION("Unknown relation type");
      }

      _stats.existing_interactions_processed++;
    } else {
      /* new interaction */
      DEBUG_PRINT("SiconosBulletCollisionManager :: New interaction\n");
      std::shared_ptr<siconos::modeling::Interaction> inter;

      auto g1 = pairA->contactor->collision_group;
      auto g2 = pairB->contactor->collision_group;
      auto nslaw = nonSmoothLaw(g1, g2);

      /* test nslaw type and then deduce the type of relation to be created */
      auto nslaw_NewtonImpactFrictionNSL =
          std::dynamic_pointer_cast<siconos::modeling::NewtonImpactFrictionNSL>(nslaw);
      auto nslaw_FremondImpactFrictionNSL =
          std::dynamic_pointer_cast<siconos::modeling::FremondImpactFrictionNSL>(nslaw);
      auto nslaw_NewtonImpactRollingFrictionNSL =
          std::dynamic_pointer_cast<siconos::modeling::NewtonImpactRollingFrictionNSL>(nslaw);

      // DEBUG_EXPR(std::cout << nslaw_NewtonImpactFrictionNSL << std::endl;);
      // DEBUG_EXPR(std::cout << nslaw_NewtonImpactRollingFrictionNSL << std::endl;);

      // we assume that this test checks if  we deal with 3D problem with RigidBodies
      // Clearly, it will not be sufficient with meshed FE bodies.
      if (nslaw && nslaw_NewtonImpactFrictionNSL  || nslaw_FremondImpactFrictionNSL) {
        if (nslaw->size() == 3) {
          DEBUG_PRINT("Creation of a relation for 3D frictional contact\n");
          auto rbdsA = std::static_pointer_cast<RigidBodyDS>(pairA->ds);
          auto rbdsB = std::static_pointer_cast<RigidBodyDS>(pairB->ds);

          auto rel = makeBulletR(rbdsA, pairA->sshape, rbdsB, pairB->sshape, *it->point);

          if (!rel) continue;

          // Fill in extra contact information
          rel->bodyShapeRecordA =
              siconos::pointers::createSPtr(*const_cast<BodyBulletShapeRecord *>(pairA));
          rel->bodyShapeRecordB =
              siconos::pointers::createSPtr(*const_cast<BodyBulletShapeRecord *>(pairB));
          rel->btObject[0] = pairA->btobject;
          rel->btObject[1] = pairB->btobject;

          // TODO cast down btshape from BodyShapeRecord-derived classes
          // rel->btShape[0] = pairA->btshape;
          // rel->btShape[1] = pairB->btshape;

          rel->updateContactPointsFromManifoldPoint(
              *it->manifold, *it->point, flip, _options->worldScale,
              rbdsA ? rbdsA : std::shared_ptr<siconos::modeling::NewtonEulerDS>(),
              rbdsB ? rbdsB : std::shared_ptr<siconos::modeling::NewtonEulerDS>());

          // We wish to be sure that no Interactions are created without
          // sufficient warning before contact.  TODO: Replace with exception or
          // flag.
          if (rel->distance() < -WARNING_TOLERANCE_AT_CREATION_INTERACTION) {
            DEBUG_PRINTF(
                "SiconosBulletCollisionManager :: Interactions must be created with positive "
                "distance (%f).\n",
                rel->distance());
            _stats.interaction_warnings++;
          }

          inter = std::make_shared<siconos::modeling::Interaction>(nslaw, rel);
          _stats.new_interactions_created++;
        } else if (nslaw && nslaw->size() == 2) {
          DEBUG_PRINT("Creation of a relation for 2D frictional contact\n");
          auto rbdsA = std::static_pointer_cast<RigidBody2dDS>(pairA->ds);
          auto rbdsB = std::static_pointer_cast<RigidBody2dDS>(pairB->ds);

          auto rel(makeBullet2dR(rbdsA, pairA->sshape, rbdsB, pairB->sshape, *it->point));

          if (!rel) continue;

          // Fill in extra contact information
          rel->bodyShapeRecordA =
              siconos::pointers::createSPtr(*const_cast<BodyBulletShapeRecord *>(pairA));
          rel->bodyShapeRecordB =
              siconos::pointers::createSPtr(*const_cast<BodyBulletShapeRecord *>(pairB));
          rel->btObject[0] = pairA->btobject;
          rel->btObject[1] = pairB->btobject;

          // TODO cast down btshape from BodyShapeRecord-derived classes
          // rel->btShape[0] = pairA->btshape;
          // rel->btShape[1] = pairB->btshape;

          rel->updateContactPointsFromManifoldPoint(
              *it->manifold, *it->point, flip, _options->worldScale,
              rbdsA ? rbdsA : std::shared_ptr<siconos::collision::RigidBody2dDS>(),
              rbdsB ? rbdsB : std::shared_ptr<siconos::collision::RigidBody2dDS>());

          // We wish to be sure that no Interactions are created without
          // sufficient warning before contact.  TODO: Replace with exception or
          // flag.
          if (rel->distance() < -WARNING_TOLERANCE_AT_CREATION_INTERACTION) {
            DEBUG_PRINTF(
                "SiconosBulletCollisionManager :: Interactions must be created with positive "
                "distance (%f).\n",
                rel->distance());
            _stats.interaction_warnings++;
          }
          DEBUG_PRINT("SiconosBulletCollisionManager :: create 2d interaction\n");
          inter = std::make_shared<siconos::modeling::Interaction>(nslaw, rel);
          _stats.new_interactions_created++;
        }

      } else if (nslaw && nslaw_NewtonImpactRollingFrictionNSL) {
        if (nslaw && nslaw->size() == 5) {
          DEBUG_PRINT("Creation of a relation for 3D Rolling frictional contact\n");
          auto rbdsA = std::static_pointer_cast<RigidBodyDS>(pairA->ds);
          auto rbdsB = std::static_pointer_cast<RigidBodyDS>(pairB->ds);

          auto rel = makeBullet5DR(rbdsA, pairA->sshape, rbdsB, pairB->sshape, *it->point);

          if (!rel) continue;

          // Fill in extra contact information
          rel->bodyShapeRecordA =
              siconos::pointers::createSPtr(*const_cast<BodyBulletShapeRecord *>(pairA));
          rel->bodyShapeRecordB =
              siconos::pointers::createSPtr(*const_cast<BodyBulletShapeRecord *>(pairB));
          rel->btObject[0] = pairA->btobject;
          rel->btObject[1] = pairB->btobject;

          // TODO cast down btshape from BodyShapeRecord-derived classes
          // rel->btShape[0] = pairA->btshape;
          // rel->btShape[1] = pairB->btshape;

          rel->updateContactPointsFromManifoldPoint(
              *it->manifold, *it->point, flip, _options->worldScale,
              rbdsA ? rbdsA : std::shared_ptr<siconos::modeling::NewtonEulerDS>(),
              rbdsB ? rbdsB : std::shared_ptr<siconos::modeling::NewtonEulerDS>());

          // We wish to be sure that no Interactions are created without
          // sufficient warning before contact.  TODO: Replace with exception or
          // flag.
          if (rel->distance() < -WARNING_TOLERANCE_AT_CREATION_INTERACTION) {
            DEBUG_PRINTF(
                "Interactions must be created with positive "
                "distance (%f).\n",
                rel->distance());
            _stats.interaction_warnings++;
          }

          inter = std::make_shared<siconos::modeling::Interaction>(nslaw, rel);
          _stats.new_interactions_created++;
        } else if (nslaw && nslaw->size() == 3) {
          DEBUG_PRINT("Creation of a relation for 2D rolling frictional contact\n");
          auto rbdsA = std::static_pointer_cast<RigidBody2dDS>(pairA->ds);
          auto rbdsB = std::static_pointer_cast<RigidBody2dDS>(pairB->ds);

          auto rel(makeBullet2d3DR(rbdsA, pairA->sshape, rbdsB, pairB->sshape, *it->point));

          if (!rel) continue;

          // Fill in extra contact information
          rel->bodyShapeRecordA =
              siconos::pointers::createSPtr(*const_cast<BodyBulletShapeRecord *>(pairA));
          rel->bodyShapeRecordB =
              siconos::pointers::createSPtr(*const_cast<BodyBulletShapeRecord *>(pairB));
          rel->btObject[0] = pairA->btobject;
          rel->btObject[1] = pairB->btobject;

          // TODO cast down btshape from BodyShapeRecord-derived classes
          // rel->btShape[0] = pairA->btshape;
          // rel->btShape[1] = pairB->btshape;

          // TODO cast down btshape from BodyShapeRecord-derived classes
          // rel->btShape[0] = pairA->btshape;
          // rel->btShape[1] = pairB->btshape;

          rel->updateContactPointsFromManifoldPoint(
              *it->manifold, *it->point, flip, _options->worldScale,
              rbdsA ? rbdsA : std::shared_ptr<siconos::collision::RigidBody2dDS>(),
              rbdsB ? rbdsB : std::shared_ptr<siconos::collision::RigidBody2dDS>());

          // We wish to be sure that no Interactions are created without
          // sufficient warning before contact.  TODO: Replace with exception or
          // flag.
          if (rel->distance() < -WARNING_TOLERANCE_AT_CREATION_INTERACTION) {
            DEBUG_PRINTF(
                "SiconosBulletCollisionManager :: Interactions must be created with positive "
                "distance (%f).\n",
                rel->distance());
            _stats.interaction_warnings++;
          }
          DEBUG_PRINT("SiconosBulletCollisionManager :: create 2d interaction\n");
          inter = std::make_shared<siconos::modeling::Interaction>(nslaw, rel);
          _stats.new_interactions_created++;
        }
      } else {
        if (nslaw && nslaw->size() == 1) {
          auto rel = std::make_shared<Bullet1DR>(siconos::pointers::createSPtr(*it->point));
          inter = std::make_shared<siconos::modeling::Interaction>(nslaw, rel);
        }
      }

      if (inter) {
        /* store interaction in the contact point data, it will be freed by the
         * Bullet callback gContactDestroyedCallback */
        /* note: storing pointer to shared_ptr! */
        it->point->m_userPersistentData =
            (void *)(new std::shared_ptr<siconos::modeling::Interaction>(inter));

        DEBUG_PRINT("SiconosBulletCollisionManager :: link the interaction\n");
        /* link bodies by the new interaction */
        simulation->link(inter, pairA->ds, pairB->ds);
      }
    }
    // getchar();
  }
// getchar();
#ifdef BULLET_TIMER
  end_old = end;
  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - end_old).count();
  std::cout << "[mechanics]  2 : creation of interaction " << elapsed << " ms" << std::endl;
#endif
  DEBUG_END(
      "siconos::collision::bullet::SiconosBulletCollisionManager::updateInteractions(std::"
      "shared_ptr<siconos::simulation::"
      "Simulation> simulation)\n");
}

void siconos::collision::bullet::SiconosBulletCollisionManager::clearOverlappingPairCache() {
  if (!_impl->_collisionWorld) return;

  btOverlappingPairCache *pairCache =
      _impl->_collisionWorld->getBroadphase()->getOverlappingPairCache();

  for (auto it = _impl->bodyShapeMap.begin(); it != _impl->bodyShapeMap.end(); it++) {
    for (auto &recV : it->second) {
      std::visit(
          [this, pairCache](auto rec) {
            if (rec->btobject)
              pairCache->cleanProxyFromPairs(rec->btobject->getBroadphaseHandle(),
                                             &*_impl->_dispatcher);
          },
          recV);
    }
  }

  // for (it = _impl->bodyShapeMap.begin(); it != _impl->bodyShapeMap.end(); it++) {
  //   std::vector<std::shared_ptr<BodyBulletShapeRecord>>::iterator rec;
  //   for (rec = it->second.begin(); rec != it->second.end(); rec++) {
  //     if ((*rec)->btobject) {
  //       pairCache->cleanProxyFromPairs((*rec)->btobject->getBroadphaseHandle(),
  //                                      &*_impl->_dispatcher);
  //     }
  //   }
  // }
}

std::vector<std::shared_ptr<siconos::collision::SiconosCollisionQueryResult>>
siconos::collision::bullet::SiconosBulletCollisionManager::lineIntersectionQuery(
    const siconos::algebra::SiconosVector &start, const siconos::algebra::SiconosVector &end,
    bool closestOnly, bool sorted) {
  std::vector<std::shared_ptr<siconos::collision::SiconosCollisionQueryResult>> result_list;

  btVector3 btstart(start(0), start(1), start(2));
  btVector3 btend(end(0), end(1), end(2));

  // Return at most one object
  if (closestOnly) {
    btCollisionWorld::ClosestRayResultCallback rayResult(btstart, btend);
    _impl->_collisionWorld->rayTest(btstart, btend, rayResult);

    if (rayResult.hasHit()) {
      const BodyShapeRecord *rec = reinterpret_cast<const BodyShapeRecord *>(
          rayResult.m_collisionObject->getUserPointer());

      if (rec) {
        auto result = std::make_shared<siconos::collision::SiconosCollisionQueryResult>();
        result->point = std::make_shared<siconos::algebra::SiconosVector>(3);
        (*result->point)(0) = rayResult.m_hitPointWorld.getX();
        (*result->point)(1) = rayResult.m_hitPointWorld.getY();
        (*result->point)(2) = rayResult.m_hitPointWorld.getZ();
        result->distance = (*result->point - start).norm();
        result->body = rec->ds;  // note: may be null for static contactors
        result->shape = rec->sshape;
        result->contactor = rec->contactor;
        result_list.push_back(result);
      } else {
        DEBUG_PRINT(
            "SiconosBulletCollisionManager :: BodyShapeRecord found by intersection was "
            "null.");
      }
    }
  }

  // Return more than one object, Bullet provides a different
  // interface
  else {
    btCollisionWorld::AllHitsRayResultCallback rayResult(btstart, btend);
    _impl->_collisionWorld->rayTest(btstart, btend, rayResult);

    if (rayResult.hasHit()) {
      for (int i = 0; i < rayResult.m_collisionObjects.size(); i++) {
        const BodyShapeRecord *rec = reinterpret_cast<const BodyShapeRecord *>(
            rayResult.m_collisionObjects[i]->getUserPointer());

        if (rec) {
          auto result(std::make_shared<SiconosCollisionQueryResult>());
          result->point = std::make_shared<siconos::algebra::SiconosVector>(3);
          result->point->resize(3);
          (*result->point)(0) = rayResult.m_hitPointWorld[i].getX();
          (*result->point)(1) = rayResult.m_hitPointWorld[i].getY();
          (*result->point)(2) = rayResult.m_hitPointWorld[i].getZ();
          result->distance = (*result->point - start).norm();
          result->body = rec->ds;  // note: null for static contactors
          result->shape = rec->sshape;
          result->contactor = rec->contactor;
          result_list.push_back(result);
        } else {
          DEBUG_PRINT("Siconos warning: BodyShapeRecord found by intersection was null.");
        }
      }
    }
  }

  if (sorted && result_list.size() > 1)
    std::sort(result_list.begin(), result_list.end(),
              [](const std::shared_ptr<siconos::collision::SiconosCollisionQueryResult> &a,
                 const std::shared_ptr<siconos::collision::SiconosCollisionQueryResult> &b) {
                return a->distance < b->distance;
              });

  return result_list;
}
