/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

#include <algorithm>
#include <MechanicsFwd.hpp>
#include <BulletSiconosFwd.hpp>

#undef SICONOS_VISITABLES
#define SICONOS_VISITABLES()                    \
  KERNEL_CLASSES()                              \
  MECHANICS_CLASSES()                           \
  REGISTER(BodyBulletShapeRecord)               \
  REGISTER(BodyBoxRecord)                       \
  REGISTER(BodySphereRecord)                    \
  REGISTER(BodyCHRecord)                        \
  REGISTER(BodyPlaneRecord)                     \
  REGISTER(BodyCylinderRecord)                  \
  REGISTER(BodyConeRecord)                      \
  REGISTER(BodyCapsuleRecord)                   \
  REGISTER(BodyMeshRecord)                      \
  REGISTER(BodyHeightRecord)                    \
  REGISTER(BodyDiskRecord)                      \
  REGISTER(BodyBox2dRecord)                     \
  REGISTER(BodyCH2dRecord)                      \


DEFINE_SPTR(UpdateShapeVisitor)

#include "SiconosBulletCollisionManager.hpp"
#include "RigidBodyDS.hpp"
#include "RigidBody2dDS.hpp"
#include "BulletR.hpp"
#include "Bullet5DR.hpp"
#include "Bullet1DR.hpp"
#include "Bullet2dR.hpp"
#include "Bullet2d3DR.hpp"
#include "StaticBody.hpp"

#include "BodyShapeRecord.hpp"


#include <map>
#include <limits>
#include <boost/format.hpp>

#include <Relation.hpp>
#include <Simulation.hpp>
#include <NonSmoothDynamicalSystem.hpp>
#include <SimulationTypeDef.hpp>
#include <NonSmoothLaw.hpp>
#include <OneStepIntegrator.hpp>
#include <NewtonImpactFrictionNSL.hpp>
#include <NewtonImpactRollingFrictionNSL.hpp>

#include <SiconosMatrix.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <NewtonEulerJointR.hpp>

#include <Question.hpp>

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunreachable-code"
#pragma clang diagnostic ignored "-Woverloaded-virtual"
#elif !(__INTEL_COMPILER || __APPLE__ )
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#endif

#include <BulletCollision/CollisionDispatch/btCollisionWorld.h>
#include <BulletCollision/Gimpact/btGImpactCollisionAlgorithm.h>
#include <BulletCollision/CollisionDispatch/btDefaultCollisionConfiguration.h>
#include <BulletCollision/BroadphaseCollision/btDbvtBroadphase.h>
#include <BulletCollision/BroadphaseCollision/btAxisSweep3.h>

#include <BulletCollision/CollisionShapes/btStaticPlaneShape.h>
#include <BulletCollision/CollisionShapes/btSphereShape.h>
#include <BulletCollision/CollisionShapes/btBoxShape.h>
#include <BulletCollision/CollisionShapes/btCylinderShape.h>
#include <BulletCollision/CollisionShapes/btConeShape.h>
#include <BulletCollision/CollisionShapes/btCapsuleShape.h>
#include <BulletCollision/CollisionShapes/btConvexHullShape.h>
#include <BulletCollision/CollisionShapes/btTriangleIndexVertexArray.h>
#include <BulletCollision/CollisionShapes/btHeightfieldTerrainShape.h>
#include <BulletCollision/CollisionDispatch/btInternalEdgeUtility.h>

// 2D shapes
#include "BulletCollision/CollisionShapes/btConvexShape.h"
#include "BulletCollision/CollisionShapes/btBox2dShape.h"
#include "BulletCollision/CollisionShapes/btConvex2dShape.h"

// 2D specific contact detection algorithm (Takane from bullet Planar2D.cpp example)
#include "BulletCollision/CollisionDispatch/btBox2dBox2dCollisionAlgorithm.h"
#include "BulletCollision/CollisionDispatch/btConvex2dConvex2dAlgorithm.h"
#include "BulletCollision/NarrowPhaseCollision/btMinkowskiPenetrationDepthSolver.h"



#include <LinearMath/btConvexHullComputer.h>
#include <BulletCollision/Gimpact/btGImpactShape.h>

#include <LinearMath/btQuaternion.h>
#include <LinearMath/btVector3.h>
//#define BULLET_TIMER
#ifdef BULLET_TIMER
#define BT_ENABLE_PROFILE 1
#include <LinearMath/btQuickprof.h>
#endif
#if defined(__clang__)
#pragma clang diagnostic pop
#elif !(__INTEL_COMPILER || __APPLE__ )
#pragma GCC diagnostic pop
#endif


// Comment this to try un-queued static contactor behaviour
#define QUEUE_STATIC_CONTACTORS 1

// We can replace the primitives by alternative implementations.  To date,
// everything works (under test conditions, very tentative) except
// btStaticPlaneShape, so we replace it with a large box.

#define USE_CONVEXHULL_FOR_BOX 1
// #define USE_CONVEXHULL_FOR_SPHERE 1
// #define USE_BOX_FOR_PLANE 1
#define USE_CONVEXHULL_FOR_PLANE 1

// Bullet types (from USE flags above)
#ifdef USE_CONVEXHULL_FOR_BOX
#define BTBOXSHAPE btConvexHullShape
#else
#define BTBOXSHAPE btBoxShape
#endif

#define BTCYLSHAPE btCylinderShape
#define BTCONSHAPE btConeShape
#define BTCAPSHAPE btCapsuleShape
#define BTCHSHAPE btConvexHullShape

#ifdef USE_CONVEXHULL_FOR_SPHERE
#define BTSPHERESHAPE btConvexHullShape
#else
#define BTSPHERESHAPE btSphereShape
#endif

#ifdef USE_BOX_FOR_PLANE
#define BTPLANESHAPE btBoxShape
#else
#ifdef USE_CONVEXHULL_FOR_PLANE
#define BTPLANESHAPE btConvexHullShape
#else
#define BTPLANESHAPE btStaticPlaneShape
#endif
#endif



// We need a bit more space to hold mesh data
class btSiconosMeshData : public btGImpactMeshShape
{
public:
  btSiconosMeshData(btStridingMeshInterface*i)
    : btGImpactMeshShape(i), btScalarVertices(nullptr) {}
  ~btSiconosMeshData()
  {
    if(btScalarVertices) delete[] btScalarVertices;
  }
  btScalar* btScalarVertices;
  SP::btTriangleIndexVertexArray btTriData;
};
#define BTMESHSHAPE btSiconosMeshData

// Similarly, place to store height matrix
class btSiconosHeightData : public btHeightfieldTerrainShape
{
public:
  btSiconosHeightData(int width, std::shared_ptr< std::vector<btScalar> > data,
                      btScalar min_height, btScalar max_height)
    : btHeightfieldTerrainShape(width, data->size()/width, data->data(),
                                0, // scale ignored for PHY_FLOAT
                                min_height, max_height,
                                2, PHY_FLOAT, false), // up = z, flip = false
      _data(data), _min_height(min_height), _max_height(max_height) {}
  std::shared_ptr< std::vector<btScalar> > _data;
  btScalar _min_height, _max_height;
};
#define BTHEIGHTSHAPE btSiconosHeightData

// Default Bullet options
SiconosBulletOptions::SiconosBulletOptions()
  : dimension(SICONOS_BULLET_3D)
  , contactBreakingThreshold(0.02)
  , contactProcessingThreshold(0.03)
  , worldScale(1.0)
  , useAxisSweep3(false)
  , clearOverlappingPairCache(false)
  , perturbationIterations(3)
  , minimumPointsPerturbationThreshold(3)
  , enableSatConvex(false)
  , enablePolyhedralContactClipping(false)
  , Depth2D(0.04)
{
}



// // We need to maintain a record associating each body with a shape,
// // contactor, and collision object for each shape type.  We also need
// // to access generic shape stuff (group, margin) by a pointer from the
// // collision callback, so we need a record base class.
// class BodyShapeRecord
// {
// public:
//   BodyShapeRecord(SP::SiconosVector b, SP::SecondOrderDS d, SP::SiconosShape sh,
//                   SP::SiconosContactor con, SP::StaticBody staticCSR)
//     : base(b), ds(d), sshape(sh), contactor(con), staticBody(staticCSR),
//       shape_version(sh->version()) {}
//   virtual ~BodyShapeRecord() {}

//   SP::SiconosVector base;
//   SP::SecondOrderDS ds;
//   SP::SiconosShape sshape;
//   SP::SiconosContactor contactor;
//   unsigned int shape_version;
//   SP::StaticBody staticBody;

//   VIRTUAL_ACCEPT_VISITORS();
// };


class BodyBulletShapeRecord : public  BodyShapeRecord
{
public:
  BodyBulletShapeRecord(SP::SiconosVector b, SP::SecondOrderDS d, SP::SiconosShape sh,
                        SP::btCollisionObject btobj,SP::SiconosContactor con, SP::StaticBody staticCSR):
    BodyShapeRecord(b, d, sh, con, staticCSR), btobject(btobj) {}
  SP::btCollisionObject btobject;
};

typedef std::map<const StaticBody*, std::vector<std::shared_ptr<BodyBulletShapeRecord> > >
StaticBodyShapeMap;

// template <typename SICONOSSHAPE, typename BULLETSHAPE>
// class BodyShapeRecordT : BodyShapeRecord
// {
// public:
//   BodyShapeRecordT(SP::SiconosVector base, SP::RigidBodyDS ds,
//                    SICONOSSHAPE sh, BULLETSHAPE btsh,
//                    SP::btCollisionObject btobj, SP::SiconosContactor con)
//     : BodyShapeRecord(base, ds, sh, btobj, con), shape(sh), btshape(btsh) {}
//   SICONOSSHAPE shape;
//   BULLETSHAPE btshape;
// };

#define SHAPE_RECORD(X, BODYDS, SICONOSSHAPE,BULLETSHAPE)     \
  class X : public BodyBulletShapeRecord,                           \
            public std::enable_shared_from_this<X> {          \
  public:                                                     \
    X(SP::SiconosVector base, BODYDS ds,                      \
      SICONOSSHAPE sh, BULLETSHAPE btsh,                      \
      SP::btCollisionObject btobj, SP::SiconosContactor con,  \
      SP::StaticBody staticCSR)                               \
      : BodyBulletShapeRecord(base, ds, sh, btobj, con, staticCSR), \
        shape(sh), btshape(btsh) {}                           \
    SICONOSSHAPE shape;                                       \
    BULLETSHAPE btshape;                                      \
    ACCEPT_VISITORS();                                        \
    };

// Body-Shape record types
SHAPE_RECORD(BodyBoxRecord, SP::RigidBodyDS, SP::SiconosBox, SP::BTBOXSHAPE);
SHAPE_RECORD(BodySphereRecord, SP::RigidBodyDS, SP::SiconosSphere, SP::BTSPHERESHAPE);
SHAPE_RECORD(BodyCHRecord, SP::RigidBodyDS, SP::SiconosConvexHull, SP::BTCHSHAPE);
SHAPE_RECORD(BodyPlaneRecord, SP::RigidBodyDS, SP::SiconosPlane, SP::BTPLANESHAPE);
SHAPE_RECORD(BodyCylinderRecord, SP::RigidBodyDS, SP::SiconosCylinder, SP::BTCYLSHAPE);
SHAPE_RECORD(BodyConeRecord, SP::RigidBodyDS, SP::SiconosCone, SP::BTCONSHAPE);
SHAPE_RECORD(BodyCapsuleRecord, SP::RigidBodyDS, SP::SiconosCapsule, SP::BTCAPSHAPE);
SHAPE_RECORD(BodyMeshRecord, SP::RigidBodyDS,  SP::SiconosMesh, std::shared_ptr<btSiconosMeshData>);
SHAPE_RECORD(BodyHeightRecord, SP::RigidBodyDS, SP::SiconosHeightMap, std::shared_ptr<btSiconosHeightData>);

SHAPE_RECORD(BodyDiskRecord, SP::RigidBody2dDS, SP::SiconosDisk,  SP::btConvex2dShape);
SHAPE_RECORD(BodyBox2dRecord, SP::RigidBody2dDS, SP::SiconosBox2d,  SP::btConvex2dShape);
SHAPE_RECORD(BodyCH2dRecord, SP::RigidBody2dDS, SP::SiconosConvexHull2d,  SP::btConvex2dShape);


typedef std::map<const SecondOrderDS*, std::vector<std::shared_ptr<BodyBulletShapeRecord> > >
BodyShapeMap;

class CollisionUpdater;

class SiconosBulletCollisionManager_impl
{
protected:
  SP::btCollisionWorld _collisionWorld;
  SP::btDefaultCollisionConfiguration _collisionConfiguration;
  SP::btCollisionDispatcher _dispatcher;
  SP::btBroadphaseInterface _broadphase;


  /* During iteration over DSs for position updates we need to access
   * btCollisionObject, so need a map DS->btXShape. */
  BodyShapeMap bodyShapeMap;

  StaticBodyShapeMap  staticBodyShapeMap;

  SP::Simulation _simulation;

  /* Create collision objects for each shape type */
  void createCollisionObject(const SP::SiconosVector base,
                             const SP::RigidBodyDS ds,
                             const SP::SiconosPlane plane,
                             const SP::SiconosContactor contactor,
                             const SP::StaticBody staticBody);
  void createCollisionObject(const SP::SiconosVector base,
                             const SP::RigidBodyDS ds,
                             const SP::SiconosSphere sphere,
                             const SP::SiconosContactor contactor,
                             const SP::StaticBody staticBody);
  void createCollisionObject(const SP::SiconosVector base,
                             const SP::RigidBodyDS ds,
                             const SP::SiconosBox box,
                             const SP::SiconosContactor contactor,
                             const SP::StaticBody staticBody);
  void createCollisionObject(const SP::SiconosVector base,
                             const SP::RigidBodyDS ds,
                             const SP::SiconosCylinder cyl,
                             const SP::SiconosContactor contactor,
                             const SP::StaticBody staticBody);
  void createCollisionObject(const SP::SiconosVector base,
                             const SP::RigidBodyDS ds,
                             const SP::SiconosCone cone,
                             const SP::SiconosContactor contactor,
                             const SP::StaticBody staticBody);
  void createCollisionObject(const SP::SiconosVector base,
                             const SP::RigidBodyDS ds,
                             const SP::SiconosCapsule capsule,
                             const SP::SiconosContactor contactor,
                             const SP::StaticBody staticBody);
  void createCollisionObject(const SP::SiconosVector base,
                             const SP::RigidBodyDS ds,
                             const SP::SiconosConvexHull ch,
                             const SP::SiconosContactor contactor,
                             const SP::StaticBody staticBody);
  void createCollisionObject(const SP::SiconosVector base,
                             const SP::RigidBodyDS ds,
                             const SP::SiconosMesh mesh,
                             const SP::SiconosContactor contactor,
                             const SP::StaticBody staticBody);
  void createCollisionObject(const SP::SiconosVector base,
                             const SP::RigidBodyDS ds,
                             const SP::SiconosHeightMap height,
                             const SP::SiconosContactor contactor,
                             const SP::StaticBody staticBody);

  void createCollisionObject(const SP::SiconosVector base,
                             const SP::RigidBody2dDS ds,
                             const SP::SiconosDisk disk,
                             const SP::SiconosContactor contactor,
                             const SP::StaticBody staticBody);
  void createCollisionObject(const SP::SiconosVector base,
                             const SP::RigidBody2dDS ds,
                             const SP::SiconosBox2d box2d,
                             const SP::SiconosContactor contactor,
                             const SP::StaticBody staticBody);
  void createCollisionObject(const SP::SiconosVector base,
                             const SP::RigidBody2dDS ds,
                             const SP::SiconosConvexHull2d ch2d,
                             const SP::SiconosContactor contactor,
                             const SP::StaticBody staticBody);


  /* Call the above functions for each shape associated with a body or contactor. */
  void createCollisionObjectsForBodyContactorSet(
    const SP::SecondOrderDS ds,
    const SP::StaticBody staticBody = SP::StaticBody(),
    const SP::SiconosVector base = SP::SiconosVector(),
    const SP::SiconosContactorSet contactor = SP::SiconosContactorSet());

  /* A helper function used to initialise new shapes, generic to the
   * shape type */
  template<typename ST, typename BT, typename DST, typename BR>
  SP::btCollisionObject  createCollisionObjectHelper(SP::SiconosVector base, const DST& ds,
                                   ST& shape, BT& btshape, BodyShapeMap& bodyShapeMap,
                                   SP::SiconosContactor contactor,
                                   StaticBodyShapeMap &StaticBodyShapeMap, SP::StaticBody staticBody);


  void updateShape(BodySphereRecord &record);
  void updateShape(BodyPlaneRecord &record);
  void updateShape(BodyBoxRecord &record);
  void updateShape(BodyCylinderRecord &record);
  void updateShape(BodyConeRecord &record);
  void updateShape(BodyCapsuleRecord &record);
  void updateShape(BodyCHRecord &record);
  void updateShape(BodyMeshRecord &record);
  void updateShape(BodyHeightRecord &record);

  void updateShape(BodyDiskRecord &record);
  void updateShape(BodyBox2dRecord &record);
  void updateShape(BodyCH2dRecord &record);

  void updateAllShapesForDS(const SecondOrderDS &bds);
  void updateShapePosition(const BodyBulletShapeRecord &record);

  /* Helper to apply an offset transform to a position and return as a
   * btTransform */
  btTransform offsetTransform(const SiconosVector& position,
                              const SiconosVector& offset);

  /** Helper to set the inertia of a NewtonEulerDS based on a
   * btCollisionShape */
  void updateContactorInertia(SP::NewtonEulerDS ds, SP::btCollisionShape btshape);

  /** Helper to set the inertia of a LagrangianDS based on a
   * btCollisionShape */
  void update2DContactorInertia(SP::LagrangianDS ds, SP::btCollisionShape btshape);

  SiconosBulletOptions &_options;

  std::vector<std::pair<SP::btCollisionObject,int>> _queuedCollisionObjects;

public:
  SiconosBulletCollisionManager_impl(SiconosBulletOptions &op) : _options(op) {}
  ~SiconosBulletCollisionManager_impl() {}

  friend class SiconosBulletCollisionManager;
  friend class CollisionUpdateVisitor;
  friend class CreateCollisionObjectShapeVisitor;
  friend class UpdateShapeVisitor;
};

SP::StaticBody SiconosBulletCollisionManager::addStaticBody(
  SP::SiconosContactorSet cs,
  SP::SiconosVector position,
  int number)
{
  SP::StaticBody rec(std::make_shared<StaticBody>());
  rec->contactorSet = cs;
  if(!position)
  {
    // Default at center
    position = std::make_shared<SiconosVector>(7);
    position->zero();
    (*position)(3) = 1.0; // we give a unit identity quarternion
  }

  rec->base = position;
  rec->number=number;
  //std::cout << "SiconosBulletCollisionManager::addStaticBody number : " << number <<  std::endl;
  _impl->createCollisionObjectsForBodyContactorSet(SP::SecondOrderDS(), rec, rec->base, cs);

  return rec;
}

void SiconosBulletCollisionManager::removeStaticBody(const SP::StaticBody& body)
{

  StaticBodyShapeMap::iterator it(_impl->staticBodyShapeMap.find(&*body));
  if(it == _impl->staticBodyShapeMap.end())
    return;

  std::vector<std::shared_ptr<BodyBulletShapeRecord> >::iterator it2;
  for(it2 = it->second.begin(); it2 != it->second.end(); it2++)
  {
    _impl->_collisionWorld->removeCollisionObject(&* (*it2)->btobject);
  }

  _impl->staticBodyShapeMap.erase(it);
}

/* We derive a specific callback for filtering the broadphase of Bullet
 * based on collision group */
struct SiconosBulletFilterCallback : public btOverlapFilterCallback
{
  InteractionManager * _interactionManager;

// return true when pairs need collision
  virtual bool needBroadphaseCollision(btBroadphaseProxy* proxy0,btBroadphaseProxy* proxy1) const
  {
    DEBUG_BEGIN("SiconosBulletFilterCallback :: needBroadphaseCollision\n");

    /* standard filter in Bullet */
    // bool collides = (proxy0->m_collisionFilterGroup & proxy1->m_collisionFilterMask) != 0;
    // collides = collides && (proxy1->m_collisionFilterGroup & proxy0->m_collisionFilterMask);

    //add some additional logic here that modified 'collides'
    SP::NonSmoothLaw nslaw = _interactionManager->nonSmoothLaw(proxy0->m_collisionFilterGroup,proxy1->m_collisionFilterGroup);
    bool collides = (bool)nslaw;

    DEBUG_END("SiconosBulletFilterCallback :: needBroadphaseCollision\n");
    return collides;
  }
};

void SiconosBulletCollisionManager::initialize_impl()
{
  _impl.reset(new SiconosBulletCollisionManager_impl(_options));

  //collision configuration contains default setup for memory, collision setup
  _impl->_collisionConfiguration.reset(
    new btDefaultCollisionConfiguration());

  if(_options.perturbationIterations > 0
      || _options.minimumPointsPerturbationThreshold > 0)
  {
    _impl->_collisionConfiguration->setConvexConvexMultipointIterations(
      _options.perturbationIterations,
      _options.minimumPointsPerturbationThreshold);
    _impl->_collisionConfiguration->setPlaneConvexMultipointIterations(
      _options.perturbationIterations,
      _options.minimumPointsPerturbationThreshold);
  }

  //use the default collision dispatcher. For parallel processing you can use a diffent dispatcher (see Extras/BulletMultiThreaded)
  _impl->_dispatcher.reset(
    new btCollisionDispatcher(&*_impl->_collisionConfiguration));


  if(_options.useAxisSweep3)
    _impl->_broadphase.reset(new btAxisSweep3(btVector3(), btVector3()));
  else
    _impl->_broadphase.reset(new btDbvtBroadphase());

  _impl->_collisionWorld.reset(
    new btCollisionWorld(&*_impl->_dispatcher, &*_impl->_broadphase,
                         &*_impl->_collisionConfiguration));


  btOverlapFilterCallback * filterCallback = new SiconosBulletFilterCallback();
  reinterpret_cast<SiconosBulletFilterCallback*>(filterCallback)->_interactionManager = this;
  _impl->_collisionWorld->getPairCache()->setOverlapFilterCallback(filterCallback);

  DEBUG_PRINTF("_options.dimension = %i", _options.dimension);

  //2D specific
  if(_options.dimension == SICONOS_BULLET_2D)
  {
    btVoronoiSimplexSolver* m_simplexSolver = new btVoronoiSimplexSolver();
    btMinkowskiPenetrationDepthSolver* m_pdSolver = new btMinkowskiPenetrationDepthSolver();

    btConvex2dConvex2dAlgorithm::CreateFunc* m_convexAlgo2d = new btConvex2dConvex2dAlgorithm::CreateFunc(m_simplexSolver,m_pdSolver);
    btBox2dBox2dCollisionAlgorithm::CreateFunc* m_box2dbox2dAlgo = new btBox2dBox2dCollisionAlgorithm::CreateFunc();

    _impl->_dispatcher->registerCollisionCreateFunc(CONVEX_2D_SHAPE_PROXYTYPE,CONVEX_2D_SHAPE_PROXYTYPE,m_convexAlgo2d);
    _impl->_dispatcher->registerCollisionCreateFunc(BOX_2D_SHAPE_PROXYTYPE,CONVEX_2D_SHAPE_PROXYTYPE,m_convexAlgo2d);
    _impl->_dispatcher->registerCollisionCreateFunc(CONVEX_2D_SHAPE_PROXYTYPE,BOX_2D_SHAPE_PROXYTYPE,m_convexAlgo2d);
    _impl->_dispatcher->registerCollisionCreateFunc(BOX_2D_SHAPE_PROXYTYPE,BOX_2D_SHAPE_PROXYTYPE,m_box2dbox2dAlgo);
  }
  else
    btGImpactCollisionAlgorithm::registerAlgorithm(&*_impl->_dispatcher);

  _impl->_collisionWorld->getDispatchInfo().m_useContinuous = false;
  _impl->_collisionWorld->getDispatchInfo().m_enableSatConvex = _options.enableSatConvex;
}

SiconosBulletCollisionManager::SiconosBulletCollisionManager()
  : SiconosCollisionManager(),_with_equality_constraints(false)
{
  initialize_impl();
}



SiconosBulletCollisionManager::SiconosBulletCollisionManager(const SiconosBulletOptions &options)
  : _with_equality_constraints(false),_options(options)
{
  initialize_impl();
}

SiconosBulletCollisionManager::~SiconosBulletCollisionManager()
{
  // unlink() will be called on all remaining
  // contact points when world is destroyed

  // must be the first de-allocated, otherwise segfault
  _impl->_collisionWorld.reset();
}

class UpdateShapeVisitor : public SiconosVisitor
{
public:
  using SiconosVisitor::visit;
  SiconosBulletCollisionManager_impl &impl;

  UpdateShapeVisitor(SiconosBulletCollisionManager_impl &_impl)
    : impl(_impl) {}

  void visit(std::shared_ptr<BodyPlaneRecord> record)
  {
    impl.updateShape(*record);
  }
  void visit(std::shared_ptr<BodySphereRecord> record)
  {
    impl.updateShape(*record);
  }
  void visit(std::shared_ptr<BodyBoxRecord> record)
  {
    impl.updateShape(*record);
  }
  void visit(std::shared_ptr<BodyCylinderRecord> record)
  {
    impl.updateShape(*record);
  }
  void visit(std::shared_ptr<BodyCapsuleRecord> record)
  {
    impl.updateShape(*record);
  }
  void visit(std::shared_ptr<BodyConeRecord> record)
  {
    impl.updateShape(*record);
  }
  void visit(std::shared_ptr<BodyCHRecord> record)
  {
    impl.updateShape(*record);
  }
  void visit(std::shared_ptr<BodyMeshRecord> record)
  {
    impl.updateShape(*record);
  }
  void visit(std::shared_ptr<BodyHeightRecord> record)
  {
    impl.updateShape(*record);
  }
  void visit(std::shared_ptr<BodyDiskRecord> record)
  {
    impl.updateShape(*record);
  }
  void visit(std::shared_ptr<BodyBox2dRecord> record)
  {
    impl.updateShape(*record);
  }
  void visit(std::shared_ptr<BodyCH2dRecord> record)
  {
    impl.updateShape(*record);
  }

};

void SiconosBulletCollisionManager_impl::updateAllShapesForDS(const SecondOrderDS &bds)
{
  SP::UpdateShapeVisitor updateShapeVisitor(new UpdateShapeVisitor(*this));
  std::vector<std::shared_ptr<BodyBulletShapeRecord> >::iterator it;
  for(it = bodyShapeMap[&bds].begin(); it != bodyShapeMap[&bds].end(); it++)
    (*it)->acceptSP(updateShapeVisitor);
}

// helper for enabling polyhedral contact clipping for shape types
// derived from btPolyhedralConvexShape
static void initPolyhedralFeatures(btPolyhedralConvexShape& btshape)
{
  btshape.initializePolyhedralFeatures();
}
static void initPolyhedralFeatures(btCollisionShape& btshape) {}

template<typename ST, typename BT, typename DST, typename BR>
SP::btCollisionObject SiconosBulletCollisionManager_impl::createCollisionObjectHelper(
  SP::SiconosVector base, const DST& ds, ST& shape, BT& btshape,
  BodyShapeMap& bodyShapeMap, SP::SiconosContactor contactor,
  StaticBodyShapeMap &StaticBodyShapeMap, SP::StaticBody staticBody)
{
  assert(base && "Collision objects must have a base position.");

  // create corresponding Bullet object and shape
  SP::btCollisionObject btobject(new btCollisionObject());

  // default parameters
  btobject->setContactProcessingThreshold(_options.contactProcessingThreshold);

  // associate the shape with the object
  btobject->setCollisionShape(&*btshape);

  // enable contact clipping for SAT
  if(_options.enablePolyhedralContactClipping)
    initPolyhedralFeatures(*btshape);

  if(!ds)
    btobject->setCollisionFlags(btCollisionObject::CF_STATIC_OBJECT);
  else
    btobject->setCollisionFlags(btCollisionObject::CF_KINEMATIC_OBJECT);

  // put it in the world
  int collisionFilterGroup = contactor->collision_group;
  int collisionFilterMask  = 1;

#ifdef QUEUE_STATIC_CONTACTORS
  if(!ds)
  {
    _queuedCollisionObjects.push_back(std::make_pair(btobject,collisionFilterGroup));
  }
  else
    _collisionWorld->addCollisionObject(&*btobject,collisionFilterGroup,collisionFilterMask);
#else
  _collisionWorld->addCollisionObject(&*btobject,collisionFilterGroup,collisionFilterMask);
#endif

  // create a record to keep track of things
  // (for static contactor, ds=nil)
  std::shared_ptr<BR> record(
    std::make_shared<BR>(base, ds, shape, btshape, btobject, contactor, staticBody));

  bodyShapeMap[ds ? &*ds : nullptr].push_back(record);


  if(staticBody)
    StaticBodyShapeMap[&*staticBody].push_back(record);


  assert(record->btobject);
  assert(record->sshape);
  assert(record->shape);
  assert(record->btshape);
  assert(record->contactor);
  assert(record->contactor->offset);
  assert(record->contactor->offset->size() == 7);

  // Allow Bullet to report colliding DSs.  We need to access it from
  // the collision callback as the record base class so down-cast it.
  btobject->setUserPointer(
    reinterpret_cast<void*>(
      static_cast<BodyShapeRecord*>(&*record)));

  // initial parameter update (change version to make something happen)
  record->shape_version -= 1;
  updateShape(*record);


  return btobject;

}

btTransform SiconosBulletCollisionManager_impl::offsetTransform(const SiconosVector& position,
    const SiconosVector& offset)
{
  /* Adjust offset position according to current rotation */
  btQuaternion rbase(position(4), position(5),
                     position(6), position(3));
  btVector3 rboffset = quatRotate(rbase, btVector3(offset(0),
                                  offset(1),
                                  offset(2)));

  /* Calculate total orientation */
  btQuaternion roffset(offset(4), offset(5), offset(6), offset(3));

  /* Set the absolute shape position */
  return btTransform(rbase * roffset,
                     btVector3(position(0), position(1), position(2)) + rboffset);
}

void SiconosBulletCollisionManager_impl::updateShapePosition(const BodyBulletShapeRecord &record)
{
  DEBUG_BEGIN("SiconosBulletCollisionManager_impl::updateShapePosition(...)\n");
  SiconosVector q(7);
  if(record.base)
  {
    DEBUG_EXPR(record.base->display(););
    if(record.base->size() ==7)
    {
      DEBUG_PRINT("3D DS\n");
      q = *record.base;
    }
    else if(record.base->size() ==3)
    {
      DEBUG_PRINT("2D DS\n");
      q(0) = record.base->getValue(0);
      q(1) = record.base->getValue(1);
      q(2) = 0.0;
      double angle = record.base->getValue(2);
      q(3) = cos(angle/2.);
      q(4) = 0.0;
      q(5) = 0.0;
      q(6) = sin(angle/2.);
    }
  }
  else
  {
    q.zero();
    q(3) = 1;
  }
  DEBUG_PRINT("Position of the shape given to bullet:")
  DEBUG_EXPR_WE(q.display(););

  btTransform t = offsetTransform(q, *record.contactor->offset);
  t.setOrigin(t.getOrigin() * _options.worldScale);
  DEBUG_PRINTF("transformation = %f,%f,%f\n", float(t.getOrigin().getX()), float(t.getOrigin().getY()), float(t.getOrigin().getZ()));
  DEBUG_PRINTF("Rotation axis = %f,%f,%f\n", float(t.getRotation().getAxis().getX()), float(t.getRotation().getAxis().getY()), float(t.getRotation().getAxis().getZ()));
  record.btobject->setWorldTransform(t);
  DEBUG_END("SiconosBulletCollisionManager_impl::updateShapePosition(...)\n");
}

void SiconosBulletCollisionManager_impl::createCollisionObject(
  const SP::SiconosVector base,
  const SP::RigidBodyDS ds,
  SP::SiconosSphere sphere,
  const SP::SiconosContactor contactor,
  const SP::StaticBody staticBody)
{
  // set radius to 1.0 and use scaling instead of setting radius
  // directly, makes it easier to change during update

#ifdef USE_CONVEXHULL_FOR_SPHERE
  // A sphere can be represented as a convex hull of a single point, with the
  // margin equal to the radius size
  SP::btConvexHullShape btsphere(new btConvexHullShape());
  {
    btsphere->addPoint(btVector3(0.0, 0.0, 0.0));
    btsphere->setMargin(0.0);
  }
#else
  SP::btSphereShape btsphere(new btSphereShape(1.0));
  btsphere->setMargin(0.0);
#endif

  // initialization
  createCollisionObjectHelper<SP::SiconosSphere, SP::BTSPHERESHAPE,
                              SP::RigidBodyDS, BodySphereRecord>
    (base, ds, sphere, btsphere, bodyShapeMap, contactor,
     staticBodyShapeMap, staticBody );
}

void SiconosBulletCollisionManager_impl::updateContactorInertia(
  SP::NewtonEulerDS ds, SP::btCollisionShape btshape)
{
  btVector3 localinertia;
  double scale_factor;
  scale_factor = 1./(_options.worldScale*_options.worldScale);
  localinertia[0] = std::numeric_limits<double>::signaling_NaN();
  localinertia[1] = std::numeric_limits<double>::signaling_NaN();
  localinertia[2] = std::numeric_limits<double>::signaling_NaN();
  btshape->calculateLocalInertia(ds->scalarMass(), localinertia);

  localinertia[0]*=scale_factor;
  localinertia[1]*=scale_factor;
  localinertia[2]*=scale_factor;
  assert(!((localinertia.x() == 0.0
            && localinertia.y() == 0.0
            && localinertia.z() == 0.0)
           || std::isinf(localinertia.x())
           || std::isinf(localinertia.y())
           || std::isinf(localinertia.z()))
         && "calculateLocalInertia() returned garbage");
  ds->setInertia(localinertia[0],
                 localinertia[1],
                 localinertia[2]);

}

void SiconosBulletCollisionManager_impl::update2DContactorInertia(
  SP::LagrangianDS ds, SP::btCollisionShape btshape)
{
  //TBD Warningx
}
void SiconosBulletCollisionManager_impl::updateShape(BodySphereRecord &record)
{
  SP::SiconosSphere sphere(record.shape);
  SP::BTSPHERESHAPE btsphere(record.btshape);

  if(sphere->version() != record.shape_version)
  {
    double r = (sphere->radius() + sphere->outsideMargin()) * _options.worldScale;

    // Update shape parameters
#ifdef USE_CONVEXHULL_FOR_SPHERE
    btsphere->setMargin(r);
#else
    btsphere->setLocalScaling(btVector3(r, r, r));

    // btSphereShape has an internal margin
    btsphere->setMargin(sphere->insideMargin() * _options.worldScale);
#endif
    SP::RigidBodyDS rbds=std::static_pointer_cast<RigidBodyDS>(record.ds);
    if(record.ds && rbds->useContactorInertia())
    {
      updateContactorInertia(rbds, btsphere);
    }

    if(record.btobject->getBroadphaseHandle())
    {
      _collisionWorld->updateSingleAabb(&*record.btobject);
//      _collisionWorld->getBroadphase()->getOverlappingPairCache()->
//        cleanProxyFromPairs(record.btobject->getBroadphaseHandle(), &*_dispatcher);
    }

    record.shape_version = sphere->version();
  }

  updateShapePosition(record);
}

void SiconosBulletCollisionManager_impl::createCollisionObject(
  const SP::SiconosVector base,
  const SP::RigidBodyDS ds,
  SP::SiconosPlane plane,
  SP::SiconosContactor contactor,
  const SP::StaticBody staticBody)
{
  // create the initial plane with default parameters
#ifdef USE_BOX_FOR_PLANE
  btScalar h = (1000 + plane->outsideMargin()) * _options.worldScale;
  SP::btBoxShape btplane(new btBoxShape(btVector3(h, h, h)));
  // Internal margin
  btplane->setMargin((plane->insideMargin() + plane->outsideMargin())
                     * _options.worldScale);
#else
#ifdef USE_CONVEXHULL_FOR_PLANE
  btScalar h = 1000 * _options.worldScale;
  const btScalar pts[] =
  {
    h, h, 0,
    h, -h, 0,
    -h, -h, 0,
    -h, h, 0,
  };
  SP::btConvexHullShape btplane(
    new btConvexHullShape(pts, 4, sizeof(pts[0])*3));

  // We ignore the desired internal margin for plane and just use a large one.
  plane->setInsideMargin(1000 * _options.worldScale);

  // External margin
  btplane->setMargin((plane->insideMargin() + plane->outsideMargin())
                     * _options.worldScale);
#else
  SP::btStaticPlaneShape btplane(
    new btStaticPlaneShape(btVector3(0, 0, 1), 0.0));
  btplane->setMargin(plane->outsideMargin() * _options.worldScale);
#endif
#endif

  // initialization
  createCollisionObjectHelper<SP::SiconosPlane, SP::BTPLANESHAPE, SP::RigidBodyDS, BodyPlaneRecord>
  (base, ds, plane, btplane, bodyShapeMap, contactor,
   staticBodyShapeMap, staticBody);
}

void SiconosBulletCollisionManager_impl::updateShape(BodyPlaneRecord& record)
{
  SP::SiconosPlane plane(record.shape);
  SP::BTPLANESHAPE btplane(record.btshape);

  SiconosVector o(7);
  o = *record.contactor->offset;

  // Adjust the offset according to plane implementation
#ifdef USE_BOX_FOR_PLANE
  o(2) -= -plane->outsideMargin() + 1000;
#else
#ifdef USE_CONVEXHULL_FOR_PLANE
  o(2) -= plane->insideMargin();
#else // USE_PLANE_FOR_PLANE
  o(2) -= plane->insideMargin();
#endif
#endif

  // Note, we do not use generic updateShapePosition for plane
  btTransform t = offsetTransform(*record.base, o);
  t.setOrigin(t.getOrigin() * _options.worldScale);
  record.btobject->setWorldTransform(t);
}

void SiconosBulletCollisionManager_impl::createCollisionObject(
  const SP::SiconosVector base,
  const SP::RigidBodyDS ds,
  SP::SiconosBox box,
  SP::SiconosContactor contactor,
  const SP::StaticBody staticBody)
{
  const btScalar half = 0.5;

#ifdef USE_CONVEXHULL_FOR_BOX
  const btScalar pts[] =
  {
    -half, half, -half,
      -half, -half, -half,
      -half, -half, half,
      -half, half, half,
      half, half, half,
      half, half, -half,
      half, -half, -half,
      half, -half, half,
    };
  SP::btConvexHullShape btbox(
    new btConvexHullShape(pts, 8, sizeof(pts[0])*3));

  // External margin
  btbox->setMargin(box->outsideMargin());
#else
  SP::btBoxShape btbox(new btBoxShape(btVector3(half, half, half)));

  // btBoxShape has an internal margin
  btbox->setMargin(box->insideMargin() * _options.worldScale);
#endif

  // initialization
  createCollisionObjectHelper<SP::SiconosBox, SP::BTBOXSHAPE, SP::RigidBodyDS, BodyBoxRecord>
  (base, ds, box, btbox, bodyShapeMap, contactor,
   staticBodyShapeMap, staticBody);
}

void SiconosBulletCollisionManager_impl::updateShape(BodyBoxRecord &record)
{
  SP::SiconosBox box(record.shape);
  SP::BTBOXSHAPE btbox(record.btshape);

  // Update shape parameters
  if(box->version() != record.shape_version)
  {
#ifdef USE_CONVEXHULL_FOR_BOX
    double m = -box->insideMargin();
#else
    double m = box->outsideMargin();
#endif

    double sx = ((*box->dimensions())(0) + m*2) * _options.worldScale;
    double sy = ((*box->dimensions())(1) + m*2) * _options.worldScale;
    double sz = ((*box->dimensions())(2) + m*2) * _options.worldScale;

    assert(sx > 0 && sy > 0 && sz > 0);

    btbox->setLocalScaling(btVector3(sx, sy, sz));
    btbox->setMargin((box->insideMargin() + box->outsideMargin()) * _options.worldScale);

    SP::RigidBodyDS rbds=std::static_pointer_cast<RigidBodyDS>(record.ds);
    if(record.ds && rbds->useContactorInertia())
      updateContactorInertia(rbds, btbox);

    if(record.btobject->getBroadphaseHandle())
    {
      _collisionWorld->updateSingleAabb(&*record.btobject);
//      _collisionWorld->getBroadphase()->getOverlappingPairCache()->
//        cleanProxyFromPairs(record.btobject->getBroadphaseHandle(), &*_dispatcher);
    }

    record.shape_version = box->version();
  }

  updateShapePosition(record);
}
void SiconosBulletCollisionManager_impl::createCollisionObject(
  const SP::SiconosVector base,
  const SP::RigidBodyDS ds,
  SP::SiconosCylinder cylinder,
  SP::SiconosContactor contactor,
  const SP::StaticBody staticBody)
{
  SP::BTCYLSHAPE btcylinder(new BTCYLSHAPE(btVector3(1.0, 1.0, 1.0)));

  // initialization
  createCollisionObjectHelper<SP::SiconosCylinder, SP::BTCYLSHAPE, SP::RigidBodyDS, BodyCylinderRecord>
  (base, ds, cylinder, btcylinder, bodyShapeMap, contactor,
   staticBodyShapeMap, staticBody);
}

void SiconosBulletCollisionManager_impl::updateShape(BodyCylinderRecord &record)
{
  SP::SiconosCylinder cyl(record.shape);
  SP::BTCYLSHAPE btcyl(record.btshape);

  // Update shape parameters
  if(cyl->version() != record.shape_version)
  {
    // Bullet cylinder has an inside margin, so we add the outside
    // margin explicitly.
    double m = cyl->outsideMargin();

    double radius = (cyl->radius() + m) * _options.worldScale;
    double length = (cyl->length()/2 + m) * _options.worldScale;

    assert(radius > 0 && length > 0);

    btcyl->setLocalScaling(btVector3(radius, length, radius));
    btcyl->setMargin((cyl->insideMargin() + cyl->outsideMargin()) * _options.worldScale);
    SP::RigidBodyDS rbds=std::static_pointer_cast<RigidBodyDS>(record.ds);
    if(record.ds && rbds->useContactorInertia())
      updateContactorInertia(rbds, btcyl);

    if(record.btobject->getBroadphaseHandle())
    {
      _collisionWorld->updateSingleAabb(&*record.btobject);
//      _collisionWorld->getBroadphase()->getOverlappingPairCache()->
//        cleanProxyFromPairs(record.btobject->getBroadphaseHandle(), &*_dispatcher);
    }

    record.shape_version = cyl->version();
  }

  updateShapePosition(record);
}


void SiconosBulletCollisionManager_impl::createCollisionObject(
  const SP::SiconosVector base,
  const SP::RigidBodyDS ds,
  SP::SiconosCone cone,
  SP::SiconosContactor contactor,
  const SP::StaticBody staticBody)
{
  SP::BTCONSHAPE btcone(new BTCONSHAPE(1.0, 1.0));

  // initialization
  createCollisionObjectHelper<SP::SiconosCone, SP::BTCONSHAPE, SP::RigidBodyDS, BodyConeRecord>
  (base, ds, cone, btcone, bodyShapeMap, contactor,
   staticBodyShapeMap, staticBody);
}

void SiconosBulletCollisionManager_impl::updateShape(BodyConeRecord &record)
{
  SP::SiconosCone cone(record.shape);
  SP::BTCONSHAPE btcone(record.btshape);

  // Update shape parameters
  if(cone->version() != record.shape_version)
  {
    // Bullet cone has an inside margin, so we add the outside
    // margin explicitly.
    double m = cone->outsideMargin();

    double radius = (cone->radius() + m) * _options.worldScale;
    double length = (cone->length()/2 + m) * _options.worldScale;

    assert(radius > 0 && length > 0);

    btcone->setLocalScaling(btVector3(radius, length, radius));
    btcone->setMargin((cone->insideMargin() + cone->outsideMargin()) * _options.worldScale);
    SP::RigidBodyDS rbds=std::static_pointer_cast<RigidBodyDS>(record.ds);
    if(record.ds && rbds->useContactorInertia())
      updateContactorInertia(rbds, btcone);

    if(record.btobject->getBroadphaseHandle())
    {
      _collisionWorld->updateSingleAabb(&*record.btobject);
//      _collisionWorld->getBroadphase()->getOverlappingPairCache()->
//        cleanProxyFromPairs(record.btobject->getBroadphaseHandle(), &*_dispatcher);
    }

    record.shape_version = cone->version();
  }

  updateShapePosition(record);
}

void SiconosBulletCollisionManager_impl::createCollisionObject(
  const SP::SiconosVector base,
  const SP::RigidBodyDS ds,
  SP::SiconosCapsule capsule,
  SP::SiconosContactor contactor,
  const SP::StaticBody staticBody)
{
  SP::BTCAPSHAPE btcapsule(new BTCAPSHAPE(1.0, 1.0));

  // initialization
  createCollisionObjectHelper<SP::SiconosCapsule, SP::BTCAPSHAPE, SP::RigidBodyDS, BodyCapsuleRecord>
  (base, ds, capsule, btcapsule, bodyShapeMap, contactor,
   staticBodyShapeMap, staticBody);
}

void SiconosBulletCollisionManager_impl::updateShape(BodyCapsuleRecord &record)
{
  SP::SiconosCapsule capsule(record.shape);
  SP::BTCAPSHAPE btcapsule(record.btshape);

  // Update shape parameters
  if(capsule->version() != record.shape_version)
  {
    // Bullet capsule has an inside margin, so we add the outside
    // margin explicitly.
    double m = capsule->outsideMargin();

    double radius = (capsule->radius() + m) * _options.worldScale;
    double length = (capsule->length()/2 + m) * _options.worldScale;

    assert(radius > 0 && length > 0);

    btcapsule->setLocalScaling(btVector3(radius, length, radius));
    btcapsule->setMargin((capsule->insideMargin() + capsule->outsideMargin()) * _options.worldScale);
    SP::RigidBodyDS rbds=std::static_pointer_cast<RigidBodyDS>(record.ds);
    if(record.ds && rbds->useContactorInertia())
      updateContactorInertia(rbds, btcapsule);

    if(record.btobject->getBroadphaseHandle())
    {
      _collisionWorld->updateSingleAabb(&*record.btobject);
//      _collisionWorld->getBroadphase()->getOverlappingPairCache()->
//        cleanProxyFromPairs(record.btobject->getBroadphaseHandle(), &*_dispatcher);
    }

    record.shape_version = capsule->version();
  }

  updateShapePosition(record);
}







void SiconosBulletCollisionManager_impl::createCollisionObject(
  const SP::SiconosVector base,
  const SP::RigidBodyDS ds,
  SP::SiconosConvexHull ch,
  SP::SiconosContactor contactor,
  const SP::StaticBody staticBody)
{
  if(!ch->vertices())
    THROW_EXCEPTION("No vertices matrix specified for convex hull.");

  if(ch->vertices()->size(1) != 3)
    THROW_EXCEPTION("Convex hull vertices matrix must have 3 columns.");

  // Copy and scale the points
  int rows = ch->vertices()->size(0);
  std::vector<btScalar> pts;
  pts.resize(rows*3);
  for(int r=0; r < rows; r++)
  {
    pts[r*3+0] = (*ch->vertices())(r, 0) * _options.worldScale;
    pts[r*3+1] = (*ch->vertices())(r, 1) * _options.worldScale;
    pts[r*3+2] = (*ch->vertices())(r, 2) * _options.worldScale;
  }

  SP::btConvexHullShape btch;
  if(ch->insideMargin() == 0)
  {
    // Create a convex hull directly with no further processing.
    // TODO: In case of worldScale=1, maybe we could avoid the copy to pts.
    btch.reset(new btConvexHullShape(&pts[0], rows, sizeof(btScalar)*3));
  }
  else
  {
    // Internal margin implemented by shrinking the hull
    // TODO: Do we need the shrink clamp? (last parameter)
    btConvexHullComputer shrinkCH;
    btScalar shrunkBy = shrinkCH.compute(&pts[0], sizeof(btScalar)*3, rows,
                                         ch->insideMargin() * _options.worldScale,
                                         0);
    if(shrunkBy < 0)
    {
      // TODO: Warning
      // "insideMargin is too large, convex hull would be too small.";
      btch.reset(new btConvexHullShape(&pts[0], rows, sizeof(btScalar)*3));
      ch->setInsideMargin(0);
    }
    else
    {
      btch.reset(new btConvexHullShape);
      for(int i=0; i < shrinkCH.vertices.size(); i++)
      {
        const btVector3 &v(shrinkCH.vertices[i]);
#if defined(BT_BULLET_VERSION) && (BT_BULLET_VERSION <= 281)
        btch->addPoint(v);
#else
        btch->addPoint(v, false);
#endif
      }
      ch->setInsideMargin(shrunkBy / _options.worldScale);
    }
  }

  // Add external margin and recalc bounding box
  btch->setMargin((ch->insideMargin() + ch->outsideMargin()) * _options.worldScale);
  btch->recalcLocalAabb();

  // initialization
  createCollisionObjectHelper<SP::SiconosConvexHull, SP::BTCHSHAPE, SP::RigidBodyDS,BodyCHRecord>
  (base, ds, ch, btch, bodyShapeMap, contactor,
   staticBodyShapeMap, staticBody);
}

void SiconosBulletCollisionManager_impl::updateShape(BodyCHRecord &record)
{
  SP::SiconosConvexHull ch(record.shape);
  SP::BTCHSHAPE btch(record.btshape);

  // Update shape parameters
  if(ch->version() != record.shape_version)
  {
    // TODO
    //btbox->setLocalScaling(btVector3(sx, sy, sz));
    btch->setMargin((ch->insideMargin() + ch->outsideMargin()) * _options.worldScale);
    SP::RigidBodyDS rbds=std::static_pointer_cast<RigidBodyDS>(record.ds);
    if(record.ds && rbds->useContactorInertia())
      updateContactorInertia(rbds, btch);

    if(record.btobject->getBroadphaseHandle())
    {
      _collisionWorld->updateSingleAabb(&*record.btobject);
      // _collisionWorld->getBroadphase()->getOverlappingPairCache()->
      // cleanProxyFromPairs(record.btobject->getBroadphaseHandle(), &*_dispatcher);
    }

    record.shape_version = ch->version();
  }

  updateShapePosition(record);
}

// If type of SiconosMatrix is the same as btScalar, we can avoid a copy
template<typename SCALAR>
std::pair<SP::btTriangleIndexVertexArray, SCALAR*>
make_bt_vertex_array(SP::SiconosMesh mesh,
                     SCALAR _s1, SCALAR _s2)
{
  assert(mesh->vertices()->size(0) == 3);
  SP::btTriangleIndexVertexArray bttris(
    std::make_shared<btTriangleIndexVertexArray>(
      mesh->indexes()->size()/3,
      (int*)mesh->indexes()->data(),
      sizeof(int)*3,
      mesh->vertices()->size(1),
      mesh->vertices()->getArray(),
      sizeof(btScalar)*3));

  return std::make_pair(bttris, (btScalar*)nullptr);
}

// If type of SiconosMatrix is not the same as btScalar, we must copy
template<typename SCALAR1, typename SCALAR2>
std::pair<SP::btTriangleIndexVertexArray, btScalar*>
make_bt_vertex_array(SP::SiconosMesh mesh, SCALAR1 _s1, SCALAR2 _s2)
{
  assert(mesh->vertices()->size(0) == 3);
  unsigned int numIndices = mesh->indexes()->size();
  unsigned int numVertices = mesh->vertices()->size(1);
  btScalar *vertices = new btScalar[numVertices*3];
  for(unsigned int i=0; i < numVertices; i++)
  {
    vertices[i*3+0] = (*mesh->vertices())(0,i);
    vertices[i*3+1] = (*mesh->vertices())(1,i);
    vertices[i*3+2] = (*mesh->vertices())(2,i);
  }
  SP::btTriangleIndexVertexArray bttris(
    std::make_shared<btTriangleIndexVertexArray>(
      numIndices/3,
      (int*)mesh->indexes()->data(),
      sizeof(int)*3,
      numVertices,
      vertices,
      sizeof(btScalar)*3));
  return std::make_pair(bttris, vertices);
}

void SiconosBulletCollisionManager_impl::createCollisionObject(
  const SP::SiconosVector base,
  const SP::RigidBodyDS ds,
  SP::SiconosMesh mesh,
  SP::SiconosContactor contactor,
  const SP::StaticBody staticBody)
{
  if(!mesh->indexes())
    THROW_EXCEPTION("No indexes matrix specified for mesh.");

  if((mesh->indexes()->size() % 3) != 0)
    THROW_EXCEPTION("Mesh indexes size must be divisible by 3.");

  if(!mesh->vertices())
    THROW_EXCEPTION("No vertices matrix specified for mesh.");

  if(mesh->vertices()->size(0) != 3)
    THROW_EXCEPTION("Convex hull vertices matrix must have 3 columns.");

  // Create Bullet triangle list, either by copying on non-copying method
  // TODO: worldScale on vertices
  std::pair<SP::btTriangleIndexVertexArray, btScalar*> datapair(
    make_bt_vertex_array(mesh, (btScalar)0, (*mesh->vertices())(0,0)));
  SP::btTriangleIndexVertexArray bttris(datapair.first);

  // Create Bullet mesh object
  SP::BTMESHSHAPE btmesh(std::make_shared<BTMESHSHAPE>(&*bttris));

  // Hold on to the data since Bullet does not make a copy
  btmesh->btTriData = bttris;
  btmesh->btScalarVertices = datapair.second;

  // Initial bound update for btGImpaceMeshShape
  btmesh->updateBound();

  // initialization
  createCollisionObjectHelper<SP::SiconosMesh, SP::BTMESHSHAPE, SP::RigidBodyDS, BodyMeshRecord>
  (base, ds, mesh, btmesh, bodyShapeMap, contactor,
   staticBodyShapeMap, staticBody);
}

void SiconosBulletCollisionManager_impl::updateShape(BodyMeshRecord &record)
{
  SP::SiconosMesh mesh(record.shape);
  SP::BTMESHSHAPE btmesh(record.btshape);

  // Update shape parameters
  if(mesh->version() != record.shape_version)
  {
    // btBvhTriangleMeshShape supports only outsideMargin.
    // TODO: support insideMargin, scale the points by their normals.
    btmesh->setMargin(mesh->outsideMargin() * _options.worldScale);
    btmesh->postUpdate();

    // TODO: Calculate inertia from a mesh.  For now we leave it at
    // identity, the user can provide an inertia if desired.

    // if (record.ds && record.ds->useContactorInertia())
    //   updateContactorInertia(record.ds, btmesh);

    if(record.btobject->getBroadphaseHandle())
    {
      _collisionWorld->updateSingleAabb(&*record.btobject);
      // _collisionWorld->getBroadphase()->getOverlappingPairCache()->
      // cleanProxyFromPairs(record.btobject->getBroadphaseHandle(), &*_dispatcher);
    }

    record.shape_version = mesh->version();
  }

  updateShapePosition(record);
}

void SiconosBulletCollisionManager_impl::createCollisionObject(
  const SP::SiconosVector base,
  const SP::RigidBodyDS ds,
  SP::SiconosHeightMap heightmap,
  SP::SiconosContactor contactor,
  const SP::StaticBody staticBody)
{
  if(!heightmap->height_data())
    THROW_EXCEPTION("No height matrix specified for heightmap.");

  SP::SiconosMatrix data = heightmap->height_data();

  if(!data || data->size(0) < 2 || data->size(1) < 2)
    THROW_EXCEPTION("Height matrix does not have sufficient dimensions "
                           "to represent a plane.");

  // Create heightfield data for Bullet.  Make a copy in case data
  // type btScalar is different from SiconosMatrix.
  // Calculate min and max value at the same time.
  double vmin = std::numeric_limits<double>::infinity();
  double vmax = -vmin;
  std::shared_ptr< std::vector<btScalar> > heightfield(
    std::make_shared< std::vector<btScalar> >());

  heightfield->resize(data->size(0) * data->size(1));

  for(unsigned int i=0; i < data->size(0); i++)
  {
    for(unsigned int j=0; j < data->size(1); j++)
    {
      double v = data->getValue(i,j);
      (*heightfield)[j*data->size(0)+i] = v;
      if(v > vmax) vmax = v;
      if(v < vmin) vmin = v;
    }
  }

  // Create Bullet height object
  SP::BTHEIGHTSHAPE btheight(std::make_shared<BTHEIGHTSHAPE>(
                               data->size(0), heightfield, vmin, vmax));

  // initialization
  SP::btCollisionObject btobject = createCollisionObjectHelper<SP::SiconosHeightMap,
                                                               SP::BTHEIGHTSHAPE,
                                                               SP::RigidBodyDS,
                                                               BodyHeightRecord>
    (base, ds, heightmap, btheight, bodyShapeMap, contactor,
     staticBodyShapeMap, staticBody);


  // this flag allows to call the call gContactAddedCallback when the callback has just been in the manifold
  // In the case of the heightmap, we use it to tweak the normal to avoid internal edge contact.
  btobject->setCollisionFlags(btCollisionObject::CF_CUSTOM_MATERIAL_CALLBACK);
}

void SiconosBulletCollisionManager_impl::updateShape(BodyHeightRecord &record)
{
  SP::SiconosHeightMap height(record.shape);
  SP::BTHEIGHTSHAPE btheight(record.btshape);

  // Update shape parameters
  if(height->version() != record.shape_version)
  {
    // btBvhTriangleHeightShape supports only outsideMargin.
    // TODO: support insideMargin, scale the points by their normals.
    btheight->setMargin((height->insideMargin() + height->outsideMargin())
                        * _options.worldScale);

    // The local scaling determines the extents of the base of the heightmap
    btheight->setLocalScaling(btVector3(
                                height->length_x() / (height->height_data()->size(0)-1),
                                height->length_y() / (height->height_data()->size(1)-1), 1));

    //TODO vertical position offset to compensate for Bullet's centering
    // TODO: Calculate the local Aabb
    //btheight->recalcLocalAabb();

    // TODO: Calculate inertia from a height.  For now we leave it at
    // identity, the user can provide an inertia if desired.

    // if (record.ds && record.ds->useContactorInertia())
    //   updateContactorInertia(record.ds, btheight);

    if(record.btobject->getBroadphaseHandle())
    {
      _collisionWorld->updateSingleAabb(&*record.btobject);
      // _collisionWorld->getBroadphase()->getOverlappingPairCache()->
      // cleanProxyFromPairs(record.btobject->getBroadphaseHandle(), &*_dispatcher);
    }

    record.shape_version = height->version();
  }

  // Like updateShapePosition(record), but we have to pre-compensate
  // the Bullet vertical centering of btHeightfieldTerrainShape.  Must
  // be done before body rotation.
  // Bullet automatically moves the center of the object to the
  // vertical center of the heightfield, so combine it here with the
  // contactor offset.
  btScalar mnz = btheight->_min_height, mxz = btheight->_max_height;
  btScalar z_offset = (mxz-mnz)/2 + mnz - height->insideMargin();
  SiconosVector o(7);
  o.zero();
  o(2) = z_offset;
  o(3) = 1;

  btTransform t = offsetTransform(*record.contactor->offset, o);
  o(0) = t.getOrigin().getX();
  o(1) = t.getOrigin().getY();
  o(2) = t.getOrigin().getZ();
  o(3) = t.getRotation().getW();
  o(4) = t.getRotation().getX();
  o(5) = t.getRotation().getY();
  o(6) = t.getRotation().getZ();

  // Now apply the combined height and contactor offset o to the base
  // transform q
  SiconosVector q(7);
  if(record.base)
    q = *record.base;
  else
  {
    q.zero();
    q(3) = 1;
  }
  t = offsetTransform(q, o);

  // DEBUG_PRINTF("updating shape position: %p(%ld) - %f, %f, %f\n",
  //              &*box,box.use_count(), q(0), q(1), q(2));

  t.setOrigin(t.getOrigin() * _options.worldScale);
  record.btobject->setWorldTransform(t);
}


///////////////////////////////////////////////////////////////////////////
// 2D shapes
///////////////////////////////////////////////////////////////////////////

void SiconosBulletCollisionManager_impl::createCollisionObject(
  const SP::SiconosVector base,
  const SP::RigidBody2dDS ds,
  SP::SiconosDisk disk,
  const SP::SiconosContactor contactor,
  const SP::StaticBody staticBody)
{
  DEBUG_BEGIN("void SiconosBulletCollisionManager_impl::createCollisionObject(..., disk, ...)\n");
  // set radius to 1.0 and use scaling instead of setting radius
  // directly, makes it easier to change during update


  //This version is ok
  double SCALING =1.0;

  btConvexShape* childShape2 = new btCylinderShapeZ(btVector3(btScalar(SCALING*1),btScalar(SCALING*1),btScalar(_options.Depth2D)));
  //btConvexShape* colShape3= new btConvex2dShape(childShape2);
  SP::btConvex2dShape btconvex2d1(new btConvex2dShape(childShape2));

  // //This version not
  // SP::btConvexShape btcylinder(new btCylinderShapeZ(btVector3(1.0, 1.0, 1.0)));
  // SP::btConvex2dShape btconvex2d(new btConvex2dShape(btcylinder.get()));

  // //This version not
  // btConvexShape* btcylinder = new btCylinderShapeZ(btVector3(1.0, 1.0, 0.04));
  // SP::btConvex2dShape btconvex2d(new btConvex2dShape(btcylinder));
  // btcylinder->setMargin(0.0);



  // initialization
  createCollisionObjectHelper<SP::SiconosDisk,
                              SP::btConvex2dShape,
                              SP::RigidBody2dDS,
                              BodyDiskRecord>
                              (base, ds, disk, btconvex2d1, bodyShapeMap, contactor,
                               staticBodyShapeMap, staticBody);
  DEBUG_END("void SiconosBulletCollisionManager_impl::createCollisionObject(..., disk, ..) \n");
}

void SiconosBulletCollisionManager_impl::updateShape(BodyDiskRecord &record)
{
  DEBUG_BEGIN("SiconosBulletCollisionManager_impl::updateShape(BodyDiskRecord &record)\n");
  SP::SiconosDisk disk(record.shape);
  SP::btConvex2dShape btconvex2d(record.btshape);

  // Update shape parameters
  if(disk->version() != record.shape_version)
  {
    // Bullet cylinder has an inside margin, so we add the outside
    // margin explicitly.
    double m = disk->outsideMargin();

    double radius = (disk->radius() + m) * _options.worldScale;

    assert(radius > 0);
    DEBUG_PRINTF("outside margin=%f \n", m);
    DEBUG_PRINTF("radius=%f \n", radius);
    DEBUG_PRINTF("_options.worldScale=%f \n", _options.worldScale);
    btconvex2d->setLocalScaling(btVector3(radius, radius, radius/25.0));
    btconvex2d->setMargin((disk->insideMargin() + disk->outsideMargin()) * _options.worldScale);


    SP::RigidBody2dDS rbds=std::static_pointer_cast<RigidBody2dDS>(record.ds);
    if(record.ds && rbds->useContactorInertia())
      update2DContactorInertia(rbds, btconvex2d);

    if(record.btobject->getBroadphaseHandle())
    {
      _collisionWorld->updateSingleAabb(&*record.btobject);
      // _collisionWorld->getBroadphase()->getOverlappingPairCache()->
      // cleanProxyFromPairs(record.btobject->getBroadphaseHandle(), &*_dispatcher);
    }

    record.shape_version = disk->version();
  }

  updateShapePosition(record);
  DEBUG_END("SiconosBulletCollisionManager_impl::updateShape(BodyDiskRecord &record)\n");
}



void SiconosBulletCollisionManager_impl::createCollisionObject(
  const SP::SiconosVector base,
  const SP::RigidBody2dDS ds,
  SP::SiconosBox2d box2d,
  const SP::SiconosContactor contactor,
  const SP::StaticBody staticBody)
{
  DEBUG_BEGIN("void SiconosBulletCollisionManager_impl::createCollisionObject(..., box2d, ...)\n");
  // set radius to 1.0 and use scaling instead of setting radius
  // directly, makes it easier to change during update


  //This version is ok
  double SCALING =1.0;
  btConvexShape* childShape0 = new btBoxShape(btVector3(btScalar(SCALING*1),btScalar(SCALING*1),btScalar(SCALING*1)));
  //btConvexShape* colShape= new btConvex2dShape(childShape0);
  SP::btConvex2dShape btconvex2d(new btConvex2dShape(childShape0));

  // initialization
  createCollisionObjectHelper<SP::SiconosBox2d,
                              SP::btConvex2dShape,
                              SP::RigidBody2dDS,
                              BodyBox2dRecord>
                              (base, ds, box2d, btconvex2d, bodyShapeMap, contactor,
                               staticBodyShapeMap, staticBody);
  DEBUG_END("void SiconosBulletCollisionManager_impl::createCollisionObject(..., box2d, ..) \n");
}

void SiconosBulletCollisionManager_impl::updateShape(BodyBox2dRecord &record)
{
  DEBUG_BEGIN("SiconosBulletCollisionManager_impl::updateShape(BodyBox2dRecord &record)\n");
  SP::SiconosBox2d box2d(record.shape);
  SP::btConvex2dShape btconvex2d(record.btshape);

  // Update shape parameters
  if(box2d->version() != record.shape_version)
  {
    // Bullet cylinder has an inside margin, so we add the outside
    // margin explicitly.
    double m = box2d->outsideMargin();

    SP::SiconosVector dimensions = box2d->dimensions();

    double width = ((*dimensions)(0) + m) * _options.worldScale;

    double height = ((*dimensions)(1) + m) * _options.worldScale;

    assert(width > 0);
    assert(height > 0);
    DEBUG_PRINTF("outside margin=%f \n", m);
    DEBUG_PRINTF("width=%f \n", width);
    DEBUG_PRINTF("height=%f \n", height);

    DEBUG_PRINTF("_options.worldScale=%f \n", _options.worldScale);
    btconvex2d->setLocalScaling(btVector3(width/2.0, height/2.0, _options.Depth2D * _options.worldScale/2.0));
    btconvex2d->setMargin((box2d->insideMargin() + box2d->outsideMargin()) * _options.worldScale);


    SP::RigidBody2dDS rbds=std::static_pointer_cast<RigidBody2dDS>(record.ds);
    if(record.ds && rbds->useContactorInertia())
      update2DContactorInertia(rbds, btconvex2d);

    if(record.btobject->getBroadphaseHandle())
    {
      _collisionWorld->updateSingleAabb(&*record.btobject);
      // _collisionWorld->getBroadphase()->getOverlappingPairCache()->
      // cleanProxyFromPairs(record.btobject->getBroadphaseHandle(), &*_dispatcher);
    }

    record.shape_version = box2d->version();
  }

  updateShapePosition(record);
  DEBUG_END("SiconosBulletCollisionManager_impl::updateShape(BodyBox2dRecord &record)\n");
}



void SiconosBulletCollisionManager_impl::createCollisionObject(
  const SP::SiconosVector base,
  const SP::RigidBody2dDS ds,
  SP::SiconosConvexHull2d ch2d,
  const SP::SiconosContactor contactor,
  const SP::StaticBody staticBody)
{
  DEBUG_BEGIN("void SiconosBulletCollisionManager_impl::createCollisionObject(..., ch2d, ...)\n");
  // set radius to 1.0 and use scaling instead of setting radius
  // directly, makes it easier to change during update
  if(!ch2d->vertices())
    THROW_EXCEPTION("No vertices matrix specified for convex hull.");

  if(ch2d->vertices()->size(1) != 2)
    THROW_EXCEPTION("2d Convex hull vertices matrix must have 2 columns.");

  // First way. We avoid to double the point
  // This works well if the  _options.worldScale is near to 1.
  // for a unknown reason
  // int rows = ch2d->vertices()->size(0);
  // std::vector<btScalar> pts;
  // pts.resize(rows*3);
  // for(int r=0; r < rows; r++)
  // {
  //   pts[r*3+0] = (*ch2d->vertices())(r, 0) * _options.worldScale;
  //   pts[r*3+1] = (*ch2d->vertices())(r, 1) * _options.worldScale;
  //   pts[r*3+2] = _options.Depth2D * _options.worldScale;
  // }

  // Second way. We double the points
  // it seems to be more robust for the contact detection
  // it avoids to find contact on the edge of the convex hull "plate"
  // Copy and scale the points
  int rows2d = ch2d->vertices()->size(0);
  int rows = rows2d *2;

  std::vector<btScalar> pts;
  pts.resize(rows*3);
  for(int r=0; r < rows2d; r++)
  {
    pts[r*3+0] = (*ch2d->vertices())(r, 0) * _options.worldScale;
    pts[r*3+1] = (*ch2d->vertices())(r, 1) * _options.worldScale;
    pts[r*3+2] = _options.Depth2D * _options.worldScale/2.0;
  }
  for(int r=rows2d; r < rows; r++)
  {
    pts[r*3+0] = (*ch2d->vertices())(r-rows2d, 0) * _options.worldScale;
    pts[r*3+1] = (*ch2d->vertices())(r-rows2d, 1) * _options.worldScale;
    pts[r*3+2] =  - _options.Depth2D * _options.worldScale/2.0;
  }

  DEBUG_EXPR_WE(
    for(int r=0; r < rows; r++)
      printf("pts[r*3+0] = %8.4e, pts[r*3+1] =%8.4e, pts[r*3+2] =%8.4e\n",pts[r*3+0], pts[r*3+1], pts[r*3+2]);
    );


  // This version is ok
  // btConvexHullShape* childShape1 = new btConvexHullShape(&pts[0],rows, sizeof(btScalar)*3);

  btConvexHullShape * btch;
  btch = new btConvexHullShape(&pts[0], rows, sizeof(btScalar)*3);  // Warning: Possible loss of memory since we cannot SP here

// Warning inside margin is not taken into account as in 3D
  if(ch2d->insideMargin() == 0)
  {
    // Create a convex hull directly with no further processing.
    // TODO: In case of worldScale=1, maybe we could avoid the copy to pts.
    btch = new btConvexHullShape(&pts[0], rows, sizeof(btScalar)*3);  // Warning: Possible loss of memory since we cannot SP here
  }
  else
  {
    // Internal margin implemented by shrinking the hull
    // TODO: Do we need the shrink clamp? (last parameter)
    btConvexHullComputer shrinkCH;
    btScalar shrunkBy = shrinkCH.compute(&pts[0], sizeof(btScalar)*3, rows,
                                         ch2d->insideMargin() * _options.worldScale,
                                         0);
    if(shrunkBy < 0)
    {
      // TODO: Warning
      // "insideMargin is too large, convex hull would be too small.";
      btch = new btConvexHullShape(&pts[0], rows, sizeof(btScalar)*3);  // Warning: Possible loss of memory since we cannot SP here
      ch2d->setInsideMargin(0);
    }
    else
    {
      btch = new btConvexHullShape;
      for(int i=0; i < shrinkCH.vertices.size(); i++)
      {
        const btVector3 &v(shrinkCH.vertices[i]);
#if defined(BT_BULLET_VERSION) && (BT_BULLET_VERSION <= 281)
        btch->addPoint(v);
#else
        btch->addPoint(v, false);
#endif
      }
      ch2d->setInsideMargin(shrunkBy / _options.worldScale);
    }
  }

  // Add external margin and recalc bounding box
  DEBUG_PRINTF("ch2d->insideMargin() = %8.4e\t, ch2d->outsideMargin() = %8.4e\n", ch2d->insideMargin(), ch2d->outsideMargin());
  btch->setMargin((ch2d->insideMargin() + ch2d->outsideMargin()) * _options.worldScale);
  btch->recalcLocalAabb();

  SP::btConvex2dShape btconvex2d(new btConvex2dShape(btch));
  btconvex2d->setMargin((ch2d->insideMargin() + ch2d->outsideMargin()) * _options.worldScale);

  // initialization
  createCollisionObjectHelper<SP::SiconosConvexHull2d,
                              SP::btConvex2dShape,
                              SP::RigidBody2dDS,
                              BodyCH2dRecord>
                              (base, ds, ch2d, btconvex2d, bodyShapeMap, contactor,
                               staticBodyShapeMap, staticBody);
  DEBUG_END("void SiconosBulletCollisionManager_impl::createCollisionObject(..., ch2d, ..) \n");
}

void SiconosBulletCollisionManager_impl::updateShape(BodyCH2dRecord &record)
{
  DEBUG_BEGIN("SiconosBulletCollisionManager_impl::updateShape(BodyCH2dRecord &record)\n");
  SP::SiconosConvexHull2d ch2d(record.shape);
  SP::btConvex2dShape btconvex2d(record.btshape);

  // Update shape parameters
  if(ch2d->version() != record.shape_version)
  {
    DEBUG_PRINT("We update the shape parameters");
    // TODO
    //btbox->setLocalScaling(btVector3(sx, sy, sz));
    btconvex2d->setMargin((ch2d->insideMargin() + ch2d->outsideMargin()) * _options.worldScale);

    SP::RigidBody2dDS rbds=std::static_pointer_cast<RigidBody2dDS>(record.ds);
    if(record.ds && rbds->useContactorInertia())
    {
      DEBUG_PRINT("We update the inertia using the contactor inertia");
      update2DContactorInertia(rbds, btconvex2d);
    }

    if(record.btobject->getBroadphaseHandle())
    {
      _collisionWorld->updateSingleAabb(&*record.btobject);
      // _collisionWorld->getBroadphase()->getOverlappingPairCache()->
      // cleanProxyFromPairs(record.btobject->getBroadphaseHandle(), &*_dispatcher);
    }

    record.shape_version = ch2d->version();
  }

  updateShapePosition(record);
  DEBUG_END("SiconosBulletCollisionManager_impl::updateShape(BodyCH2dRecord &record)\n");
}




class CreateCollisionObjectShapeVisitor : public SiconosVisitor
{
public:
  using SiconosVisitor::visit;
  SiconosBulletCollisionManager_impl &impl;
  const SP::SecondOrderDS ds;
  const SP::SiconosVector base;
  SP::SiconosContactor contactor;
  SP::StaticBody staticBody;

  CreateCollisionObjectShapeVisitor(SiconosBulletCollisionManager_impl &_impl,
                                    const SP::SecondOrderDS _ds,
                                    const SP::SiconosVector _base,
                                    const SP::StaticBody _staticBody)
    : impl(_impl), ds(_ds), base(_base), staticBody(_staticBody) {}

  void visit(SP::SiconosPlane shape)
  {
    SP::RigidBodyDS rbds =  std::static_pointer_cast<RigidBodyDS>(ds);
    impl.createCollisionObject(base, rbds, shape, contactor, staticBody);
  }
  void visit(SP::SiconosSphere shape)
  {
    SP::RigidBodyDS rbds =  std::static_pointer_cast<RigidBodyDS>(ds);
    impl.createCollisionObject(base, rbds, shape, contactor, staticBody);
  }
  void visit(SP::SiconosBox shape)
  {
    SP::RigidBodyDS rbds =  std::static_pointer_cast<RigidBodyDS>(ds);
    impl.createCollisionObject(base, rbds, shape, contactor, staticBody);
  }
  void visit(SP::SiconosCylinder shape)
  {
    SP::RigidBodyDS rbds =  std::static_pointer_cast<RigidBodyDS>(ds);
    impl.createCollisionObject(base, rbds, shape, contactor, staticBody);
  }
  void visit(SP::SiconosCone shape)
  {
    SP::RigidBodyDS rbds =  std::static_pointer_cast<RigidBodyDS>(ds);
    impl.createCollisionObject(base, rbds, shape, contactor, staticBody);
  }
  void visit(SP::SiconosCapsule shape)
  {
    SP::RigidBodyDS rbds =  std::static_pointer_cast<RigidBodyDS>(ds);
    impl.createCollisionObject(base, rbds, shape, contactor, staticBody);
  }
  void visit(SP::SiconosConvexHull shape)
  {
    SP::RigidBodyDS rbds =  std::static_pointer_cast<RigidBodyDS>(ds);
    impl.createCollisionObject(base, rbds, shape, contactor, staticBody);
  }
  void visit(SP::SiconosMesh shape)
  {
    SP::RigidBodyDS rbds =  std::static_pointer_cast<RigidBodyDS>(ds);
    impl.createCollisionObject(base, rbds, shape, contactor, staticBody);
  }
  void visit(SP::SiconosHeightMap shape)
  {
    SP::RigidBodyDS rbds =  std::static_pointer_cast<RigidBodyDS>(ds);
    impl.createCollisionObject(base, rbds, shape, contactor, staticBody);
  }
  void visit(SP::SiconosDisk shape)
  {
    SP::RigidBody2dDS rb2dds =  std::static_pointer_cast<RigidBody2dDS>(ds);
    impl.createCollisionObject(base, rb2dds, shape, contactor, staticBody);
  }
  void visit(SP::SiconosBox2d shape)
  {
    SP::RigidBody2dDS rb2dds =  std::static_pointer_cast<RigidBody2dDS>(ds);
    impl.createCollisionObject(base, rb2dds, shape, contactor, staticBody);
  }
  void visit(SP::SiconosConvexHull2d shape)
  {
    SP::RigidBody2dDS rb2dds =  std::static_pointer_cast<RigidBody2dDS>(ds);
    impl.createCollisionObject(base, rb2dds, shape, contactor, staticBody);
  }
};

void SiconosBulletCollisionManager_impl::createCollisionObjectsForBodyContactorSet(
  const SP::SecondOrderDS ds,
  const SP::StaticBody staticBody,
  SP::SiconosVector base,
  SP::SiconosContactorSet contactors)
{
  DEBUG_BEGIN("SiconosBulletCollisionManager_impl::createCollisionObjectsForBodyContactorSet(...)\n");
  // ensure consistency between ds and base and contactor -- if they
  // can be taken from the DS, ensure nothing else is provided.
  assert(!(ds && base) && "Provide either DS or base, but not both!");
  assert(!(ds && contactors) && "Provide either DS or contactor set, but not both!");

  // of course, we need at least one of the two combinations
  assert(ds || (base && contactors));

  SP::RigidBodyDS rbds(std::dynamic_pointer_cast<RigidBodyDS>(ds));
  SP::RigidBody2dDS rb2dds(std::dynamic_pointer_cast<RigidBody2dDS>(ds));


  SP::SiconosContactorSet con(contactors);


  if(rbds)
  {
    DEBUG_PRINT("RigidBodyDS case");
    con = rbds->contactors();
    base = rbds->q();
  }
  if(rb2dds)
  {
    DEBUG_PRINT("RigidBody2dDS case");
    con = rb2dds->contactors();
    base = rb2dds->q();
  }
  // if ((!rbds) and (!rb2dds))
  // {
  //   std::cout << "createCollisionObjectsForBodyContactorSet for static objects" << std::endl;
  //   if(!staticBody)
  //   {
  //     std::cout << "createCollisionObjectsForBodyContactorSet for static objects. a staticBody is required" << std::endl;
  //   }
  // }

  if(!con)
  {
    DEBUG_PRINT("No contactors");
    DEBUG_BEGIN("SiconosBulletCollisionManager_impl::createCollisionObjectsForBodyContactorSet(...)\n");
    return;
  }
  std::shared_ptr<CreateCollisionObjectShapeVisitor>
    ccosv(new CreateCollisionObjectShapeVisitor(*this, ds, base, staticBody));

  /* Call createCollisionObject for each shape type using the visitor
   * defined above */
  std::vector< SP::SiconosContactor >::const_iterator it;
  for(it=con->begin(); it!=con->end(); it++)
  {
    // special collision group -1 = do not collide, thus we can skip
    // creation of associated collision objects
    if((*it)->collision_group == -1) continue;

    // otherwise visit the object with createCollisionObject
    ccosv->contactor = *it;
    ccosv->contactor->shape->acceptSP(ccosv);
  }
  DEBUG_END("SiconosBulletCollisionManager_impl::createCollisionObjectsForBodyContactorSet(...)\n");
}

void SiconosBulletCollisionManager::removeBody(const SP::SecondOrderDS& body)
{

  BodyShapeMap::iterator it(_impl->bodyShapeMap.find(&*body));
  if(it == _impl->bodyShapeMap.end())
    return;

  std::vector<std::shared_ptr<BodyBulletShapeRecord> >::iterator it2;
  for(it2 = it->second.begin(); it2 != it->second.end(); it2++)
  {
    _impl->_collisionWorld->removeCollisionObject(&* (*it2)->btobject);
  }

  _impl->bodyShapeMap.erase(it);
}

/** This class allows to iterate over all the contact points in a
 *  btCollisionWorld, returning a tuple containing the two btCollisionObjects
 *  and the btManifoldPoint.  To be called after
 *  performDiscreteCollisionDetection().
 */
class IterateContactPoints
{
public:
  SP::btCollisionWorld world;

  IterateContactPoints(SP::btCollisionWorld _world)
    : world(_world) {}

  struct ContactPointTuple
  {
    const btCollisionObject* objectA;
    const btCollisionObject* objectB;
    btManifoldPoint* point;
    btPersistentManifold* manifold;
  };

  class iterator
  {
  protected:
    SP::btCollisionWorld world;
    ContactPointTuple data;
    unsigned int numManifolds;
    unsigned int numContacts;
    unsigned int manifold_index;
    unsigned int contact_index;
    iterator(SP::btCollisionWorld _world)
    {
      numManifolds = 0;
      manifold_index = -1;
      contact_index = -1;
      numContacts = 0;
      world = _world;
      numManifolds = world->getDispatcher()->getNumManifolds();
      ++(*this);
    }
  public:
    const ContactPointTuple& operator*()
    {
      return data;
    };
    const ContactPointTuple* operator->()
    {
      return &data;
    };

    iterator()
    {
      numManifolds = 0;
      manifold_index = -1;
      contact_index = -1;
      numContacts = 0;
    }

    iterator& operator++()
    {
      if(numManifolds == 0)
        return *this;
      contact_index ++;
      while(contact_index >= numContacts)
      {
        manifold_index ++;
        if(manifold_index < numManifolds)
        {
          data.manifold = world->getDispatcher()->
                          getManifoldByIndexInternal(manifold_index);
          data.objectA = data.manifold->getBody0();
          data.objectB = data.manifold->getBody1();
          numContacts = data.manifold->getNumContacts();
          contact_index = 0;
        }
        else
        {
          numManifolds = 0;
          return *this;
        }
      }
      data.point = &(data.manifold->getContactPoint(contact_index));
      return *this;
    };

    bool operator!=(const iterator &it)
    {
      if(it.numManifolds==0) return numManifolds!=0;
      return data.objectA != it.data.objectA
             || data.objectB != it.data.objectB
             || data.point != it.data.point;
    };
    friend class IterateContactPoints;
  };

  iterator begin()
  {
    return iterator(world);
  };

  iterator end()
  {
    return iterator();
  };
};

// called once for each contact point as it is destroyed
Simulation* SiconosBulletCollisionManager::gSimulation;
bool SiconosBulletCollisionManager::bulletContactClear(void* userPersistentData)
{
  /* note: stored pointer to shared_ptr! */
  SP::Interaction *p_inter = (SP::Interaction*)userPersistentData;
  assert(p_inter!=NULL && "Contact point's stored (SP::Interaction*) is null!");
  DEBUG_PRINTF("unlinking interaction %p, number %zu \n", &**p_inter, (*p_inter)->number());

  // SP::BulletR rel_bulletR(std::dynamic_pointer_cast<BulletR>((*p_inter)->relation()));
  // SP::Bullet5DR rel_bullet5DR(std::dynamic_pointer_cast<Bullet5DR>((*p_inter)->relation()));
  // SP::Bullet2dR rel_bullet2dR(std::dynamic_pointer_cast<Bullet2dR>((*p_inter)->relation()));
  // SP::Bullet2d3DR rel_bullet2d3DR(std::dynamic_pointer_cast<Bullet2d3DR>((*p_inter)->relation()));
  // if (rel_bulletR)
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



static void siconosBulletAdjustInternalEdgeContacts(btManifoldPoint& cp, const btCollisionObjectWrapper* colObj0Wrap, const btCollisionObjectWrapper* colObj1Wrap, int partId0, int index0)
{
	//btAssert(colObj0->getCollisionShape()->getShapeType() == TRIANGLE_SHAPE_PROXYTYPE);
	if (colObj0Wrap->getCollisionShape()->getShapeType() != TRIANGLE_SHAPE_PROXYTYPE)
		return;


	btTriangleInfoMap* triangleInfoMapPtr = nullptr;

	if (colObj0Wrap->getCollisionObject()->getCollisionShape()->getShapeType() == TERRAIN_SHAPE_PROXYTYPE)
	{
		btHeightfieldTerrainShape* heightfield = (btHeightfieldTerrainShape*)colObj0Wrap->getCollisionObject()->getCollisionShape();
		triangleInfoMapPtr = heightfield->getTriangleInfoMap();

		btVector3 newNormal = btVector3(0, 0, 1);

		const btTriangleShape* tri_shape = static_cast<const btTriangleShape*>(colObj0Wrap->getCollisionShape());
		btVector3 tri_normal;
		tri_shape->calcNormal(tri_normal);
		newNormal = tri_normal;
		//					cp.m_distance1 = cp.m_distance1 * newNormal.dot(cp.m_normalWorldOnB);
    btVector3 oldNormal =  	cp.m_normalWorldOnB;

    //printf("old normal %e\t%e\t%e\n", oldNormal.x(),  oldNormal.y(), oldNormal.z());
    //printf("new normal %e\t%e\t%e\n", newNormal.x(),  newNormal.y(), newNormal.z());
    //printf("cp.m_distance1 = %e\n", cp.m_distance1 );

    btScalar cosine =  oldNormal.dot(newNormal);
    //printf("cosine %e\n", cosine);
    if (cosine < 0.0)
    {
      newNormal = -1.0*tri_normal;
      cosine =  oldNormal.dot(newNormal);
    }
    //btScalar diff  = oldNormal.distance(newNormal);
    //printf("diff %e\n", diff);
    if ((1.0 - cosine) > 3e-03 ) // around 5 degrees
    {
      //printf("--------------------------------------> change edge  normal to triangle normal\n");
      cp.m_normalWorldOnB = newNormal;
    }
    else return;

		// Reproject collision point along normal. (what about cp.m_distance1?)
		cp.m_positionWorldOnB = cp.m_positionWorldOnA - cp.m_normalWorldOnB * cp.m_distance1;
		cp.m_localPointB = colObj0Wrap->getWorldTransform().invXform(cp.m_positionWorldOnB);
		return;
	}
}



bool SiconosBulletCollisionManager::bulletContactAddedCallback(btManifoldPoint& cp, const btCollisionObjectWrapper* colObj0Wrap, int partId0, int index0, const btCollisionObjectWrapper* colObj1Wrap, int partId1, int index1)
{
  //printf("--------- bulletContactAddedCallback start\n");
	//btAdjustInternalEdgeContacts(cp, colObj1Wrap, colObj0Wrap, partId1, index1);
  siconosBulletAdjustInternalEdgeContacts(cp, colObj1Wrap, colObj0Wrap, partId1, index1);
  //printf("--------- bulletContactAddedCallback end\n");

	return true;
}
SP::BulletR SiconosBulletCollisionManager::makeBulletR(SP::RigidBodyDS ds1,
    SP::SiconosShape shape1,
    SP::RigidBodyDS ds2,
    SP::SiconosShape shape2,
    const btManifoldPoint &p)
{
  return std::make_shared<BulletR>();
}
SP::Bullet5DR SiconosBulletCollisionManager::makeBullet5DR(SP::RigidBodyDS ds1,
    SP::SiconosShape shape1,
    SP::RigidBodyDS ds2,
    SP::SiconosShape shape2,
    const btManifoldPoint &p)
{
  return std::make_shared<Bullet5DR>();
}
SP::Bullet2dR SiconosBulletCollisionManager::makeBullet2dR(SP::RigidBody2dDS ds1,
    SP::SiconosShape shape1,
    SP::RigidBody2dDS ds2,
    SP::SiconosShape shape2,
    const btManifoldPoint &p)
{
  return std::make_shared<Bullet2dR>();
}
SP::Bullet2d3DR SiconosBulletCollisionManager::makeBullet2d3DR(SP::RigidBody2dDS ds1,
    SP::SiconosShape shape1,
    SP::RigidBody2dDS ds2,
    SP::SiconosShape shape2,
    const btManifoldPoint &p)
{
  return std::make_shared<Bullet2d3DR>();
}
class CollisionUpdateVisitor : public SiconosVisitor
{
public:
  using SiconosVisitor::visit;
  SiconosBulletCollisionManager_impl &impl;

  CollisionUpdateVisitor(SiconosBulletCollisionManager_impl& _impl)
    : impl(_impl) {}

  void visit(SP::RigidBodyDS bds)
  {
    if(bds->contactors())
    {
      if(impl.bodyShapeMap.find(&*bds) == impl.bodyShapeMap.end())
      {
        impl.createCollisionObjectsForBodyContactorSet(bds);
      }
      impl.updateAllShapesForDS(*bds);
    }
  }
  void visit(SP::RigidBody2dDS bds)
  {
    if(bds->contactors())
    {
      if(impl.bodyShapeMap.find(&*bds) == impl.bodyShapeMap.end())
      {
        impl.createCollisionObjectsForBodyContactorSet(bds);
      }
      impl.updateAllShapesForDS(*bds);
    }
  }

};

void SiconosBulletCollisionManager::updateInteractions(SP::Simulation simulation)
{
  DEBUG_BEGIN("SiconosBulletCollisionManager::updateInteractions(SP::Simulation simulation)\n");
#ifdef BULLET_TIMER
  CProfileManager::Start_Profile("bullet_profile.txt");
  CProfileManager::Reset();
#endif
  // -2. update collision objects from all RigidBodyDS dynamical systems
#ifdef BULLET_TIMER
  std::chrono::time_point<std::chrono::system_clock> start, end, end_old;
  start = std::chrono::system_clock::now();
#endif

  SP::SiconosVisitor updateVisitor(new CollisionUpdateVisitor(*_impl));
  simulation->nonSmoothDynamicalSystem()->visitDynamicalSystems(updateVisitor);
#ifdef BULLET_TIMER
  end = std::chrono::system_clock::now();
  int elapsed = std::chrono::duration_cast<std::chrono::milliseconds> (end-start).count();
  std::cout << "-2 : visit " << elapsed << " ms" << std::endl;
#endif
  // Clear cache automatically before collision detection if requested
  if(_options.clearOverlappingPairCache)
    clearOverlappingPairCache();

  if(! _impl->_queuedCollisionObjects.empty())
  {
    int collisionFilterMask  = 1;

    std::vector<std::pair<SP::btCollisionObject,int>>::iterator it;
    for(it = _impl->_queuedCollisionObjects.begin();
        it != _impl->_queuedCollisionObjects.end();
        ++ it)
    {
      std::pair<SP::btCollisionObject,int> p = *it;
      int collisionFilterGroup  = p.second;
      SP::btCollisionObject collisionObject = p.first;
      _impl->_collisionWorld->addCollisionObject(&*collisionObject,collisionFilterGroup,collisionFilterMask);
    }
    _impl->_queuedCollisionObjects.clear();
  }

  // -1. reset statistical counters
  resetStatistics();

#ifdef BULLET_TIMER
  end_old=end;
  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>
    (end-end_old).count();
  std::cout << "-1 : addCollisionObject " << elapsed << " ms" << std::endl;
#endif
  // 0. set up bullet callbacks
  gSimulation = &*simulation;
  gContactDestroyedCallback = this->bulletContactClear;
  gContactAddedCallback = this->bulletContactAddedCallback;

  // Important parameter controlling contact point making and breaking
  gContactBreakingThreshold = _options.contactBreakingThreshold;

  // 1. perform bullet collision detection
  _impl->_collisionWorld->performDiscreteCollisionDetection();
#ifdef BULLET_TIMER
  end_old =end;
  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>
    (end-end_old).count();
  std::cout << "1 : collisionDectection" << elapsed << " ms" << std::endl;
#endif


#ifdef BULLET_TIMER
  CProfileManager::dumpAll();
  CProfileManager::Stop_Profile();
#endif

  DEBUG_PRINT("SiconosBulletCollisionManager :: iterating contact points:\n");
  //getchar();
  // 2. deleted contact points have been removed from the graph during the
  //    bullet collision detection callbacks

  // 3. for each contact point, if there is no interaction, create one
  IterateContactPoints t(_impl->_collisionWorld);
  IterateContactPoints::iterator it, itend=t.end();
  DEBUG_EXPR_WE(
    int num_contact_points =0;
    for(it=t.begin(); it!=itend; ++it)  num_contact_points++;
    std::cout << "Number of contacts points detected by bullet: " << num_contact_points << std::endl; );

  for(it=t.begin(); it!=itend; ++it)
  {
    DEBUG_PRINTF("\n\n\nSiconosBulletCollisionManager ::   -- %p, %p, %p\n", it->objectA, it->objectB, it->point);

    // Get the RigidBodyDS and SiconosShape pointers

    const BodyBulletShapeRecord *pairA, *pairB;
    pairA = reinterpret_cast<const BodyBulletShapeRecord*>(it->objectA->getUserPointer());
    pairB = reinterpret_cast<const BodyBulletShapeRecord*>(it->objectB->getUserPointer());
    assert(pairA && pairB && "btCollisionObject had a null user pointer!");

    // The first pair will always be the non-static object
    // As a consequence, if there is a static body, it is always associated with second pair pairB
    bool flip = false;
    if(pairB->ds && !pairA->ds)
    {
      pairA = reinterpret_cast<const BodyBulletShapeRecord*>(it->objectB->getUserPointer());
      pairB = reinterpret_cast<const BodyBulletShapeRecord*>(it->objectA->getUserPointer());
      flip = true;
    }
    DEBUG_PRINTF("SiconosBulletCollisionManager :: flip = %i \n", flip);
    // If both collision objects belong to the same body (or no body),
    // no interaction is created.
    if(pairA->ds == pairB->ds)
      continue;

    // If the two bodies are already connected by another type of
    // relation (e.g. EqualityCondition == they have a joint between
    // them), then don't create contact constraints, because it leads
    // to an ill-conditioned problem.

    DEBUG_EXPR_WE(
      if (pairA->ds && pairB->ds)
      {
        DEBUG_PRINTF("SiconosBulletCollisionManager ::   -- ds1 :  %zu,  ds2: %zu\n",
                     pairA->ds->number(),
                     pairB->ds->number());
      }
      if (pairA->ds && pairB->staticBody)
      {
        DEBUG_PRINTF("SiconosBulletCollisionManager ::   -- ds1 :  %zu  staticbody: %i\n",
                     pairA->ds->number(),
                     pairB->staticBody->number);
      }
      );

    DEBUG_PRINTF("SiconosBulletCollisionManager :: _with_equality_constraints  -- %i\n", _with_equality_constraints);


    if(_with_equality_constraints && pairA->ds && pairB->ds)
    {
      InteractionsGraph::VIterator ui, uiend;
      SP::InteractionsGraph indexSet0 = simulation->nonSmoothDynamicalSystem()->topology()->indexSet0();
      bool match = false;
      for(std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
      {
        SP::Interaction inter(indexSet0->bundle(*ui));
        SP::SecondOrderDS ds1(std::dynamic_pointer_cast<SecondOrderDS>(
                                indexSet0->properties(*ui).source));
        SP::SecondOrderDS ds2(std::dynamic_pointer_cast<SecondOrderDS>(
                                indexSet0->properties(*ui).target));
        if(ds1 && ds2 && (((&*ds1==&*pairA->ds) && (&*ds2==&*pairB->ds))
                          || ((&*ds1==&*pairB->ds) && (&*ds2==&*pairA->ds))))
        {
          SP::BulletR br(std::dynamic_pointer_cast<BulletR>(inter->relation()));
          DEBUG_EXPR(std::cout << "br" << br << std::endl;);
          if(!br)
          {
            DEBUG_PRINT("Only match on non-BulletR interactions, i.e. non-contact relations\n");
            SP::NewtonEulerJointR jr(
              std::dynamic_pointer_cast<NewtonEulerJointR>(inter->relation()));

            /* If it is a joint, check the joint self-collide property */
            if(jr && !jr->allowSelfCollide())
              match = true;

            /* If any non-contact relation is found, both bodies must
             * allow self-collide */
            // We need to check for other type of dynamical systems.
            SP::RigidBodyDS rbdsA =  std::static_pointer_cast<RigidBodyDS>(pairA->ds);
            SP::RigidBodyDS rbdsB =  std::static_pointer_cast<RigidBodyDS>(pairB->ds);
            if(!rbdsA->allowSelfCollide() || !rbdsB->allowSelfCollide())
              match = true;
          }
          if(match) break;
        }
      }
      if(match)
        continue;
    }
    DEBUG_PRINTF("SiconosBulletCollisionManager :: it->point->m_userPersistentData  %p \n", it->point->m_userPersistentData);
    if(it->point->m_userPersistentData)
    {
      /* interaction already exists */
      DEBUG_PRINT("SiconosBulletCollisionManager :: interaction already exists \n");
      SP::Interaction *p_inter =
        (SP::Interaction*)it->point->m_userPersistentData;


      SP::BulletR rel_bulletR(std::dynamic_pointer_cast<BulletR>((*p_inter)->relation()));
      SP::Bullet5DR rel_bullet5DR(std::dynamic_pointer_cast<Bullet5DR>((*p_inter)->relation()));
      SP::Bullet2dR rel_bullet2dR(std::dynamic_pointer_cast<Bullet2dR>((*p_inter)->relation()));
      SP::Bullet2d3DR rel_bullet2d3DR(std::dynamic_pointer_cast<Bullet2d3DR>((*p_inter)->relation()));

      if(rel_bulletR || rel_bullet5DR)
      {
        DEBUG_PRINT("SiconosBulletCollisionManager :: BulletR case || rel_bullet5DR\n");
        // We need to check for other type of dynamical systems.
        SP::RigidBodyDS rbdsA =  std::static_pointer_cast<RigidBodyDS>(pairA->ds);
        SP::RigidBodyDS rbdsB =  std::static_pointer_cast<RigidBodyDS>(pairB->ds);

        /* update the relation */
        SP::BulletR rel(std::static_pointer_cast<BulletR>((*p_inter)->relation()));
        rel->updateContactPointsFromManifoldPoint(*it->manifold, *it->point,
            flip, _options.worldScale,
            rbdsA,
            rbdsB ? rbdsB
            : SP::NewtonEulerDS());
      }
      else if(rel_bullet2dR)
      {
        DEBUG_PRINT("SiconosBulletCollisionManager :: Bullet2dR case");
        // We need to check for other type of dynamical systems.
        SP::RigidBody2dDS rbdsA =  std::static_pointer_cast<RigidBody2dDS>(pairA->ds);
        SP::RigidBody2dDS rbdsB =  std::static_pointer_cast<RigidBody2dDS>(pairB->ds);

        /* update the relation */
        rel_bullet2dR->updateContactPointsFromManifoldPoint(*it->manifold, *it->point,
            flip, _options.worldScale,
            rbdsA,
            rbdsB ? rbdsB
            : SP::RigidBody2dDS());
      }
      else if(rel_bullet2d3DR)
      {
        DEBUG_PRINT("SiconosBulletCollisionManager :: Bullet2d3DR case");
        // We need to check for other type of dynamical systems.
        SP::RigidBody2dDS rbdsA =  std::static_pointer_cast<RigidBody2dDS>(pairA->ds);
        SP::RigidBody2dDS rbdsB =  std::static_pointer_cast<RigidBody2dDS>(pairB->ds);

        /* update the relation */
        rel_bullet2d3DR->updateContactPointsFromManifoldPoint(*it->manifold, *it->point,
            flip, _options.worldScale,
            rbdsA,
            rbdsB ? rbdsB
            : SP::RigidBody2dDS());
      }

      else
      {
        THROW_EXCEPTION("Unknown relation type");
      }


      _stats.existing_interactions_processed ++;
    }
    else
    {
      /* new interaction */
      DEBUG_PRINT("SiconosBulletCollisionManager :: New interaction\n");
      SP::Interaction inter;

      int g1 = pairA->contactor->collision_group;
      int g2 = pairB->contactor->collision_group;
      SP::NonSmoothLaw nslaw = nonSmoothLaw(g1,g2);

      /* test nslaw type and then deduce the type of relation to be created */
      SP::NewtonImpactFrictionNSL nslaw_NewtonImpactFrictionNSL(std::dynamic_pointer_cast<NewtonImpactFrictionNSL>(nslaw));
      SP::NewtonImpactRollingFrictionNSL nslaw_NewtonImpactRollingFrictionNSL(std::dynamic_pointer_cast<NewtonImpactRollingFrictionNSL>(nslaw));

      // DEBUG_EXPR(std::cout << nslaw_NewtonImpactFrictionNSL << std::endl;);
      // DEBUG_EXPR(std::cout << nslaw_NewtonImpactRollingFrictionNSL << std::endl;);

      // we assume that this test checks if  we deal with 3D problem with RigidBodies
      // Clearly, it will not be sufficient with meshed FE bodies.
      if(nslaw && nslaw_NewtonImpactFrictionNSL)
      {
        if(nslaw->size() == 3)
        {
          DEBUG_PRINT("Creation of a relation for 3D frictional contact\n");
          SP::RigidBodyDS rbdsA =  std::static_pointer_cast<RigidBodyDS>(pairA->ds);
          SP::RigidBodyDS rbdsB =  std::static_pointer_cast<RigidBodyDS>(pairB->ds);

          SP::BulletR rel(makeBulletR(rbdsA, pairA->sshape,
                                      rbdsB, pairB->sshape,
                                      *it->point));

          if(!rel) continue;

          // Fill in extra contact information
          rel->bodyShapeRecordA = createSPtrBodyBulletShapeRecord(*const_cast<BodyBulletShapeRecord*>(pairA));
          rel->bodyShapeRecordB = createSPtrBodyBulletShapeRecord(*const_cast<BodyBulletShapeRecord*>(pairB));
          rel->btObject[0] = pairA->btobject;
          rel->btObject[1] = pairB->btobject;

          // TODO cast down btshape from BodyShapeRecord-derived classes
          // rel->btShape[0] = pairA->btshape;
          // rel->btShape[1] = pairB->btshape;

          rel->updateContactPointsFromManifoldPoint(*it->manifold, *it->point,
              flip, _options.worldScale,
              rbdsA ? rbdsA : SP::NewtonEulerDS(),
              rbdsB ? rbdsB : SP::NewtonEulerDS());

          // We wish to be sure that no Interactions are created without
          // sufficient warning before contact.  TODO: Replace with exception or
          // flag.
          if(rel->distance() < 0.0)
          {
            DEBUG_PRINTF("SiconosBulletCollisionManager :: Interactions must be created with positive "
                         "distance (%f).\n", rel->distance());
            _stats.interaction_warnings ++;
          }

          inter = std::make_shared<Interaction>(nslaw, rel);
          _stats.new_interactions_created ++;
        }
        else if(nslaw && nslaw->size() == 2)
        {
          DEBUG_PRINT("Creation of a relation for 2D frictional contact\n");
          SP::RigidBody2dDS rbdsA =  std::static_pointer_cast<RigidBody2dDS>(pairA->ds);
          SP::RigidBody2dDS rbdsB =  std::static_pointer_cast<RigidBody2dDS>(pairB->ds);

          SP::Bullet2dR rel(makeBullet2dR(rbdsA, pairA->sshape,
                                          rbdsB, pairB->sshape,
                                          *it->point));

          if(!rel) continue;

           // Fill in extra contact information
          rel->bodyShapeRecordA = createSPtrBodyBulletShapeRecord(*const_cast<BodyBulletShapeRecord*>(pairA));
          rel->bodyShapeRecordB = createSPtrBodyBulletShapeRecord(*const_cast<BodyBulletShapeRecord*>(pairB));
          rel->btObject[0] = pairA->btobject;
          rel->btObject[1] = pairB->btobject;

          // TODO cast down btshape from BodyShapeRecord-derived classes
          // rel->btShape[0] = pairA->btshape;
          // rel->btShape[1] = pairB->btshape;

          rel->updateContactPointsFromManifoldPoint(*it->manifold, *it->point,
              flip, _options.worldScale,
              rbdsA ? rbdsA : SP::RigidBody2dDS(),
              rbdsB ? rbdsB : SP::RigidBody2dDS());

          // We wish to be sure that no Interactions are created without
          // sufficient warning before contact.  TODO: Replace with exception or
          // flag.
          if(rel->distance() < 0.0)
          {
            DEBUG_PRINTF("SiconosBulletCollisionManager :: Interactions must be created with positive "
                         "distance (%f).\n", rel->distance());
            _stats.interaction_warnings ++;
          }
          DEBUG_PRINT("SiconosBulletCollisionManager :: create 2d interaction\n");
          inter = std::make_shared<Interaction>(nslaw, rel);
          _stats.new_interactions_created ++;
        }

      }
      else if(nslaw && nslaw_NewtonImpactRollingFrictionNSL)
      {
        if(nslaw && nslaw->size() == 5)
        {
          DEBUG_PRINT("Creation of a relation for 3D Rolling frictional contact\n");
          SP::RigidBodyDS rbdsA =  std::static_pointer_cast<RigidBodyDS>(pairA->ds);
          SP::RigidBodyDS rbdsB =  std::static_pointer_cast<RigidBodyDS>(pairB->ds);

          SP::Bullet5DR rel(makeBullet5DR(rbdsA, pairA->sshape,
                                          rbdsB, pairB->sshape,
                                          *it->point));

          if(!rel) continue;

          // Fill in extra contact information
          rel->bodyShapeRecordA = createSPtrBodyBulletShapeRecord(*const_cast<BodyBulletShapeRecord*>(pairA));
          rel->bodyShapeRecordB = createSPtrBodyBulletShapeRecord(*const_cast<BodyBulletShapeRecord*>(pairB));
          rel->btObject[0] = pairA->btobject;
          rel->btObject[1] = pairB->btobject;

          // TODO cast down btshape from BodyShapeRecord-derived classes
          // rel->btShape[0] = pairA->btshape;
          // rel->btShape[1] = pairB->btshape;

          rel->updateContactPointsFromManifoldPoint(*it->manifold, *it->point,
              flip, _options.worldScale,
              rbdsA ? rbdsA : SP::NewtonEulerDS(),
              rbdsB ? rbdsB : SP::NewtonEulerDS());

          // We wish to be sure that no Interactions are created without
          // sufficient warning before contact.  TODO: Replace with exception or
          // flag.
          if(rel->distance() < 0.0)
          {
            DEBUG_PRINTF("Interactions must be created with positive "
                         "distance (%f).\n", rel->distance());
            _stats.interaction_warnings ++;
          }

          inter = std::make_shared<Interaction>(nslaw, rel);
          _stats.new_interactions_created ++;
        }
        else if(nslaw && nslaw->size() == 3)
        {
          DEBUG_PRINT("Creation of a relation for 2D rolling frictional contact\n");
          SP::RigidBody2dDS rbdsA =  std::static_pointer_cast<RigidBody2dDS>(pairA->ds);
          SP::RigidBody2dDS rbdsB =  std::static_pointer_cast<RigidBody2dDS>(pairB->ds);

          SP::Bullet2d3DR rel(makeBullet2d3DR(rbdsA, pairA->sshape,
                                              rbdsB, pairB->sshape,
                                              *it->point));

          if(!rel) continue;

          // Fill in extra contact information
          rel->bodyShapeRecordA = createSPtrBodyBulletShapeRecord(*const_cast<BodyBulletShapeRecord*>(pairA));
          rel->bodyShapeRecordB = createSPtrBodyBulletShapeRecord(*const_cast<BodyBulletShapeRecord*>(pairB));
          rel->btObject[0] = pairA->btobject;
          rel->btObject[1] = pairB->btobject;

          // TODO cast down btshape from BodyShapeRecord-derived classes
          // rel->btShape[0] = pairA->btshape;
          // rel->btShape[1] = pairB->btshape;

          // TODO cast down btshape from BodyShapeRecord-derived classes
          // rel->btShape[0] = pairA->btshape;
          // rel->btShape[1] = pairB->btshape;

          rel->updateContactPointsFromManifoldPoint(*it->manifold, *it->point,
              flip, _options.worldScale,
              rbdsA ? rbdsA : SP::RigidBody2dDS(),
              rbdsB ? rbdsB : SP::RigidBody2dDS());

          // We wish to be sure that no Interactions are created without
          // sufficient warning before contact.  TODO: Replace with exception or
          // flag.
          if(rel->distance() < 0.0)
          {
            DEBUG_PRINTF("SiconosBulletCollisionManager :: Interactions must be created with positive "
                         "distance (%f).\n", rel->distance());
            _stats.interaction_warnings ++;
          }
          DEBUG_PRINT("SiconosBulletCollisionManager :: create 2d interaction\n");
          inter = std::make_shared<Interaction>(nslaw, rel);
          _stats.new_interactions_created ++;
        }
      }
      else
      {
        if(nslaw && nslaw->size() == 1)
        {
          SP::Bullet1DR rel(
            std::make_shared<Bullet1DR>(
              createSPtrbtManifoldPoint(*it->point)));
          inter = std::make_shared<Interaction>(nslaw, rel);
        }
      }

      if(inter)
      {
        /* store interaction in the contact point data, it will be freed by the
         * Bullet callback gContactDestroyedCallback */
        /* note: storing pointer to shared_ptr! */
        it->point->m_userPersistentData = (void*)(new SP::Interaction(inter));
        DEBUG_PRINT("SiconosBulletCollisionManager :: link the interaction\n");
        /* link bodies by the new interaction */
        simulation->link(inter, pairA->ds, pairB->ds);
      }
    }
    //getchar();
  }
#ifdef BULLET_TIMER
  end_old =end;
  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>
    (end-end_old).count();
  std::cout << "2 : creation of interaction " << elapsed << " ms" << std::endl;
#endif
  DEBUG_END("SiconosBulletCollisionManager::updateInteractions(SP::Simulation simulation)\n");
}

void SiconosBulletCollisionManager::clearOverlappingPairCache()
{
  if(!_impl->_collisionWorld) return;

  BodyShapeMap::iterator it;
  btOverlappingPairCache *pairCache =
    _impl->_collisionWorld->getBroadphase()->getOverlappingPairCache();

  for(it = _impl->bodyShapeMap.begin(); it != _impl->bodyShapeMap.end(); it++)
  {
    std::vector< std::shared_ptr<BodyBulletShapeRecord> >::iterator rec;
    for(rec = it->second.begin(); rec != it->second.end(); rec++)
    {
      if((*rec)->btobject)
      {
        pairCache-> cleanProxyFromPairs((*rec)->btobject->getBroadphaseHandle(),
                                        &*_impl->_dispatcher);
      }
    }
  }
}

/** Implement comparison for std::sort(). */
static
bool cmpQueryResult(const SP::SiconosCollisionQueryResult &a,
                    const SP::SiconosCollisionQueryResult &b)
{
  return a->distance < b->distance;
}

std::vector<SP::SiconosCollisionQueryResult>
SiconosBulletCollisionManager::lineIntersectionQuery(const SiconosVector& start,
    const SiconosVector& end,
    bool closestOnly,
    bool sorted)
{
  std::vector<SP::SiconosCollisionQueryResult> result_list;

  btVector3 btstart(start(0), start(1), start(2));
  btVector3 btend(end(0), end(1), end(2));

  // Return at most one object
  if(closestOnly)
  {
    btCollisionWorld::ClosestRayResultCallback rayResult(btstart, btend);
    _impl->_collisionWorld->rayTest(btstart, btend, rayResult);

    if(rayResult.hasHit())
    {
      const BodyShapeRecord *rec =
        reinterpret_cast<const BodyShapeRecord*>(
          rayResult.m_collisionObject->getUserPointer());

      if(rec)
      {
        SP::SiconosCollisionQueryResult result(
          std::make_shared<SiconosCollisionQueryResult>());
        result->point.resize(3);
        result->point.setValue(0, rayResult.m_hitPointWorld.getX());
        result->point.setValue(1, rayResult.m_hitPointWorld.getY());
        result->point.setValue(2, rayResult.m_hitPointWorld.getZ());
        result->distance = (result->point - start).norm2();
        result->body = rec->ds; // note: may be null for static contactors
        result->shape = rec->sshape;
        result->contactor = rec->contactor;
        result_list.push_back(result);
      }
      else
      {
        DEBUG_PRINT("SiconosBulletCollisionManager :: BodyShapeRecord found by intersection was null.");
      }
    }
  }

  // Return more than one object, Bullet provides a different
  // interface
  else
  {
    btCollisionWorld::AllHitsRayResultCallback rayResult(btstart, btend);
    _impl->_collisionWorld->rayTest(btstart, btend, rayResult);

    if(rayResult.hasHit())
    {
      for(int i=0; i < rayResult.m_collisionObjects.size(); i++)
      {
        const BodyShapeRecord *rec =
          reinterpret_cast<const BodyShapeRecord*>(
            rayResult.m_collisionObjects[i]->getUserPointer());

        if(rec)
        {
          SP::SiconosCollisionQueryResult result(
            std::make_shared<SiconosCollisionQueryResult>());
          result->point.resize(3);
          result->point.setValue(0, rayResult.m_hitPointWorld[i].getX());
          result->point.setValue(1, rayResult.m_hitPointWorld[i].getY());
          result->point.setValue(2, rayResult.m_hitPointWorld[i].getZ());
          result->distance = (result->point - start).norm2();
          result->body = rec->ds; // note: null for static contactors
          result->shape = rec->sshape;
          result->contactor = rec->contactor;
          result_list.push_back(result);
        }
        else
        {
          DEBUG_PRINT("Siconos warning: BodyShapeRecord found by intersection was null.");
        }
      }
    }
  }

  if(sorted && result_list.size() > 1)
    std::sort(result_list.begin(), result_list.end(), cmpQueryResult);

  return result_list;
}
