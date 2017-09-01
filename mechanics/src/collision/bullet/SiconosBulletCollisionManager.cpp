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

/*! \file SiconosBulletCollisionManager.cpp
  \brief Definition of a Bullet-based interaction handler for contact
  detection.
*/

#include <algorithm>
#include <MechanicsFwd.hpp>
#include <BulletSiconosFwd.hpp>

#undef SICONOS_VISITABLES
#define SICONOS_VISITABLES()                    \
  KERNEL_CLASSES()                              \
  MECHANICS_CLASSES()                           \
  REGISTER(BodyShapeRecord)                     \
  REGISTER(BodyBoxRecord)                       \
  REGISTER(BodySphereRecord)                    \
  REGISTER(BodyCHRecord)                        \
  REGISTER(BodyPlaneRecord)                     \
  REGISTER(BodyCylinderRecord)                  \
  REGISTER(BodyMeshRecord)                      \
  REGISTER(BodyHeightRecord)

DEFINE_SPTR(UpdateShapeVisitor)

#include "SiconosBulletCollisionManager.hpp"
#include "BodyDS.hpp"
#include "BulletR.hpp"
#include "BulletFrom1DLocalFrameR.hpp"

#include <map>
#include <limits>
#include <boost/format.hpp>
#include <boost/make_shared.hpp>
#include <CxxStd.hpp>

#include <Model.hpp>
#include <Relation.hpp>
#include <Simulation.hpp>
#include <NonSmoothDynamicalSystem.hpp>
#include <SimulationTypeDef.hpp>
#include <NonSmoothLaw.hpp>
#include <OneStepIntegrator.hpp>
#include <NewtonImpactFrictionNSL.hpp>
#include <FrictionContact.hpp>
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
#include <BulletCollision/CollisionShapes/btConvexHullShape.h>
#include <BulletCollision/CollisionShapes/btBvhTriangleMeshShape.h>
#include <BulletCollision/CollisionShapes/btTriangleIndexVertexArray.h>
#include <BulletCollision/CollisionShapes/btHeightfieldTerrainShape.h>
#include <LinearMath/btConvexHullComputer.h>

#include <LinearMath/btQuaternion.h>
#include <LinearMath/btVector3.h>

#if defined(__clang__)
#pragma clang diagnostic pop
#elif !(__INTEL_COMPILER || __APPLE__ )
#pragma GCC diagnostic pop
#endif

//#define DEBUG_MESSAGES 1
#include <debug.h>

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
class btSiconosMeshData : public btBvhTriangleMeshShape {
public:
  btSiconosMeshData(btStridingMeshInterface*i, bool b)
    : btBvhTriangleMeshShape(i,b), btScalarVertices(0) {}
  ~btSiconosMeshData() { if (btScalarVertices) delete[] btScalarVertices; }
  btScalar* btScalarVertices;
  SP::btTriangleIndexVertexArray btTriData;
};
#define BTMESHSHAPE btSiconosMeshData

// Similarly, place to store height matrix
class btSiconosHeightData : public btHeightfieldTerrainShape {
public:
  btSiconosHeightData(int width, std11::shared_ptr< std::vector<btScalar> > data,
                      btScalar min_height, btScalar max_height)
    : btHeightfieldTerrainShape(width, data->size()/width, data->data(),
                                0, // scale ignored for PHY_FLOAT
                                min_height, max_height,
                                2, PHY_FLOAT, false), // up = z, flip = false
      _data(data), _min_height(min_height), _max_height(max_height) {}
  std11::shared_ptr< std::vector<btScalar> > _data;
  btScalar _min_height, _max_height;
};
#define BTHEIGHTSHAPE btSiconosHeightData

// Default Bullet options
SiconosBulletOptions::SiconosBulletOptions()
  : contactBreakingThreshold(0.02)
  , contactProcessingThreshold(0.03)
  , worldScale(1.0)
  , useAxisSweep3(false)
  , clearOverlappingPairCache(false)
  , perturbationIterations(3)
  , minimumPointsPerturbationThreshold(3)
{
}

// We need to maintain a record associating each body with a shape,
// contactor, and collision object for each shape type.  We also need
// to access generic shape stuff (group, margin) by a pointer from the
// collision callback, so we need a record base class.
class BodyShapeRecord
{
public:
  BodyShapeRecord(SP::SiconosVector b, SP::BodyDS d, SP::SiconosShape sh,
                  SP::btCollisionObject btobj, SP::SiconosContactor con)
    : base(b), ds(d), sshape(sh), btobject(btobj), contactor(con),
      shape_version(sh->version()) {}
  virtual ~BodyShapeRecord() {}

  SP::SiconosVector base;
  SP::BodyDS ds;
  SP::SiconosShape sshape;
  SP::btCollisionObject btobject;
  SP::SiconosContactor contactor;
  unsigned int shape_version;

  VIRTUAL_ACCEPT_VISITORS();
};

template <typename SICONOSSHAPE, typename BULLETSHAPE>
class BodyShapeRecordT : BodyShapeRecord
{
public:
  BodyShapeRecordT(SP::SiconosVector base, SP::BodyDS ds,
                   SICONOSSHAPE sh, BULLETSHAPE btsh,
                   SP::btCollisionObject btobj, SP::SiconosContactor con)
    : BodyShapeRecord(base, ds, sh, btobj, con), shape(sh), btshape(btsh) {}
  SICONOSSHAPE shape;
  BULLETSHAPE btshape;
};

#define SHAPE_RECORD(X,SICONOSSHAPE,BULLETSHAPE)                        \
  class X : public BodyShapeRecord,                                     \
            public std11::enable_shared_from_this<X> {                  \
  public:                                                               \
    X(SP::SiconosVector base, SP::BodyDS ds,                            \
      SICONOSSHAPE sh, BULLETSHAPE btsh,                                \
      SP::btCollisionObject btobj, SP::SiconosContactor con)            \
      : BodyShapeRecord(base, ds, sh, btobj, con), shape(sh), btshape(btsh) {} \
    SICONOSSHAPE shape;                                                 \
    BULLETSHAPE btshape;                                                \
    ACCEPT_VISITORS();                                                  \
  };

// Body-Shape record types
SHAPE_RECORD(BodyBoxRecord, SP::SiconosBox, SP::BTBOXSHAPE);
SHAPE_RECORD(BodySphereRecord, SP::SiconosSphere, SP::BTSPHERESHAPE);
SHAPE_RECORD(BodyCHRecord, SP::SiconosConvexHull, SP::BTCHSHAPE);
SHAPE_RECORD(BodyPlaneRecord, SP::SiconosPlane, SP::BTPLANESHAPE);
SHAPE_RECORD(BodyCylinderRecord, SP::SiconosCylinder, SP::BTCYLSHAPE);
SHAPE_RECORD(BodyMeshRecord, SP::SiconosMesh, std11::shared_ptr<btSiconosMeshData>);
SHAPE_RECORD(BodyHeightRecord, SP::SiconosHeightMap, std11::shared_ptr<btSiconosHeightData>);

typedef std::map<const BodyDS*, std::vector<std11::shared_ptr<BodyShapeRecord> > >
  BodyShapeMap;

/** For associating static contactor sets and their offsets.
 * Pointer to this is what is returned as the opaque and unique
 * StaticContactorSetID so that they can be removed. */
class StaticContactorSetRecord
{
public:
  SP::SiconosContactorSet contactorSet;
  SP::SiconosVector base;
};
namespace SP {
  typedef std11::shared_ptr<StaticContactorSetRecord> StaticContactorSetRecord;
};

class CollisionUpdater;

class SiconosBulletCollisionManager_impl
{
protected:
  SP::btCollisionWorld _collisionWorld;
  SP::btDefaultCollisionConfiguration _collisionConfiguration;
  SP::btCollisionDispatcher _dispatcher;
  SP::btBroadphaseInterface _broadphase;

  /* Static contactor sets may be repeated with different positions,
   * thus each one is assocated with a list of base positions and
   * collision objects. */
  std::map< StaticContactorSetRecord*, SP::StaticContactorSetRecord >
    _staticContactorSetRecords;

  /* During iteration over DSs for position updates we need to access
   * btCollisionObject, so need a map DS->btXShape. */
  BodyShapeMap bodyShapeMap;

  SP::Simulation _simulation;

  /* Create collision objects for each shape type */
  void createCollisionObject(const SP::SiconosVector base,
                             const SP::BodyDS ds,
                             const SP::SiconosPlane plane,
                             const SP::SiconosContactor contactor);
  void createCollisionObject(const SP::SiconosVector base,
                             const SP::BodyDS ds,
                             const SP::SiconosSphere sphere,
                             const SP::SiconosContactor contactor);
  void createCollisionObject(const SP::SiconosVector base,
                             const SP::BodyDS ds,
                             const SP::SiconosBox box,
                             const SP::SiconosContactor contactor);
  void createCollisionObject(const SP::SiconosVector base,
                             const SP::BodyDS ds,
                             const SP::SiconosCylinder cyl,
                             const SP::SiconosContactor contactor);
  void createCollisionObject(const SP::SiconosVector base,
                             const SP::BodyDS ds,
                             const SP::SiconosConvexHull ch,
                             const SP::SiconosContactor contactor);
  void createCollisionObject(const SP::SiconosVector base,
                             const SP::BodyDS ds,
                             const SP::SiconosMesh mesh,
                             const SP::SiconosContactor contactor);
  void createCollisionObject(const SP::SiconosVector base,
                             const SP::BodyDS ds,
                             const SP::SiconosHeightMap height,
                             const SP::SiconosContactor contactor);

  /* Call the above functions for each shape associated with a body or contactor. */
  void createCollisionObjectsForBodyContactorSet(
    const SP::BodyDS ds,
    const SP::SiconosVector base = SP::SiconosVector(),
    const SP::SiconosContactorSet contactor = SP::SiconosContactorSet());

  /* A helper function used to initialise new shapes, generic to the
   * shape type */
  template<typename ST, typename BT, typename BR>
  void createCollisionObjectHelper(SP::SiconosVector base, SP::BodyDS ds,
                                   ST& shape, BT& btshape, BodyShapeMap& bodyShapeMap,
                                   SP::SiconosContactor contactor);

  void updateShape(BodySphereRecord &record);
  void updateShape(BodyPlaneRecord &record);
  void updateShape(BodyBoxRecord &record);
  void updateShape(BodyCylinderRecord &record);
  void updateShape(BodyCHRecord &record);
  void updateShape(BodyMeshRecord &record);
  void updateShape(BodyHeightRecord &record);

  void updateAllShapesForDS(const BodyDS &bds);
  void updateShapePosition(const BodyShapeRecord &record);

  /* Helper to apply an offset transform to a position and return as a
   * btTransform */
  btTransform offsetTransform(const SiconosVector& position,
                              const SiconosVector& offset);

  /** Helper to set the inertia of a NewtonEulerDS based on a
   * btCollisionShape */
  void updateContactorInertia(SP::NewtonEulerDS ds, SP::btCollisionShape btshape);

  SiconosBulletOptions &_options;

  std::vector<SP::btCollisionObject> _queuedCollisionObjects;

public:
  SiconosBulletCollisionManager_impl(SiconosBulletOptions &op) : _options(op) {}
  ~SiconosBulletCollisionManager_impl() {}

  friend class SiconosBulletCollisionManager;
  friend class CollisionUpdateVisitor;
  friend class CreateCollisionObjectShapeVisitor;
  friend class UpdateShapeVisitor;
};

SiconosCollisionManager::StaticContactorSetID
SiconosBulletCollisionManager::insertStaticContactorSet(
  SP::SiconosContactorSet cs,
  SP::SiconosVector position)
{
  SP::StaticContactorSetRecord rec(std11::make_shared<StaticContactorSetRecord>());
  rec->contactorSet = cs;
  if (!position)
  { // Default at center
    position = std11::make_shared<SiconosVector>(7);
    position->zero();
    (*position)(3) = 1.0;
  }
  rec->base = position;
  impl->createCollisionObjectsForBodyContactorSet(SP::BodyDS(), rec->base, cs);
  impl->_staticContactorSetRecords[&*rec] = rec;
  return static_cast<SiconosBulletCollisionManager::StaticContactorSetID>(&*rec);
}

bool SiconosBulletCollisionManager::removeStaticContactorSet(StaticContactorSetID id)
{
  StaticContactorSetRecord *recptr = static_cast<StaticContactorSetRecord *>(id);
  if (impl->_staticContactorSetRecords.find(recptr)
      == impl->_staticContactorSetRecords.end())
    return false;

  SP::StaticContactorSetRecord rec(impl->_staticContactorSetRecords[recptr]);
  // TODO
  assert(0 && "removeStaticContactorSet not implemented.");
  return false;
}

void SiconosBulletCollisionManager::initialize_impl()
{
  impl.reset(new SiconosBulletCollisionManager_impl(_options));
  impl->_collisionConfiguration.reset(
    new btDefaultCollisionConfiguration());

  if (_options.perturbationIterations > 0
      || _options.minimumPointsPerturbationThreshold > 0)
  {
    impl->_collisionConfiguration->setConvexConvexMultipointIterations(
      _options.perturbationIterations,
      _options.minimumPointsPerturbationThreshold);
    impl->_collisionConfiguration->setPlaneConvexMultipointIterations(
      _options.perturbationIterations,
      _options.minimumPointsPerturbationThreshold);
  }

  impl->_dispatcher.reset(
    new btCollisionDispatcher(&*impl->_collisionConfiguration));

  if (_options.useAxisSweep3)
    impl->_broadphase.reset(new btAxisSweep3(btVector3(), btVector3()));
  else
    impl->_broadphase.reset(new btDbvtBroadphase());

  impl->_collisionWorld.reset(
    new btCollisionWorld(&*impl->_dispatcher, &*impl->_broadphase,
                         &*impl->_collisionConfiguration));

  btGImpactCollisionAlgorithm::registerAlgorithm(&*impl->_dispatcher);
  impl->_collisionWorld->getDispatchInfo().m_useContinuous = false;
}

SiconosBulletCollisionManager::SiconosBulletCollisionManager()
  : SiconosCollisionManager()
{
  initialize_impl();
}

SiconosBulletCollisionManager::SiconosBulletCollisionManager(const SiconosBulletOptions &options)
  : _options(options)
{
  initialize_impl();
}

SiconosBulletCollisionManager::~SiconosBulletCollisionManager()
{
  // unlink() will be called on all remaining
  // contact points when world is destroyed

  // must be the first de-allocated, otherwise segfault
  impl->_collisionWorld.reset();
}

class UpdateShapeVisitor : public SiconosVisitor
{
public:
  using SiconosVisitor::visit;
  SiconosBulletCollisionManager_impl &impl;

  UpdateShapeVisitor(SiconosBulletCollisionManager_impl &_impl)
    : impl(_impl) {}

  void visit(std11::shared_ptr<BodyPlaneRecord> record)
    { impl.updateShape(*record); }
  void visit(std11::shared_ptr<BodySphereRecord> record)
    { impl.updateShape(*record); }
  void visit(std11::shared_ptr<BodyBoxRecord> record)
    { impl.updateShape(*record); }
  void visit(std11::shared_ptr<BodyCylinderRecord> record)
    { impl.updateShape(*record); }
  void visit(std11::shared_ptr<BodyCHRecord> record)
    { impl.updateShape(*record); }
  void visit(std11::shared_ptr<BodyMeshRecord> record)
    { impl.updateShape(*record); }
  void visit(std11::shared_ptr<BodyHeightRecord> record)
    { impl.updateShape(*record); }
};

void SiconosBulletCollisionManager_impl::updateAllShapesForDS(const BodyDS &bds)
{
  SP::UpdateShapeVisitor updateShapeVisitor(new UpdateShapeVisitor(*this));
  std::vector<std11::shared_ptr<BodyShapeRecord> >::iterator it;
  for (it = bodyShapeMap[&bds].begin(); it != bodyShapeMap[&bds].end(); it++)
    (*it)->acceptSP(updateShapeVisitor);
}

template<typename ST, typename BT, typename BR>
void SiconosBulletCollisionManager_impl::createCollisionObjectHelper(
  SP::SiconosVector base, SP::BodyDS ds, ST& shape, BT& btshape,
  BodyShapeMap& bodyShapeMap, SP::SiconosContactor contactor)
{
  assert(base && "Collision objects must have a base position.");

  // create corresponding Bullet object and shape
  SP::btCollisionObject btobject(new btCollisionObject());

  // default parameters
  btobject->setContactProcessingThreshold(_options.contactProcessingThreshold);

  // associate the shape with the object
  btobject->setCollisionShape(&*btshape);

  if (!ds)
    btobject->setCollisionFlags(btCollisionObject::CF_STATIC_OBJECT);
  else
    btobject->setCollisionFlags(btCollisionObject::CF_KINEMATIC_OBJECT);

  // put it in the world
  #ifdef QUEUE_STATIC_CONTACTORS
  if (!ds)
    _queuedCollisionObjects.push_back(btobject);
  else
    _collisionWorld->addCollisionObject(&*btobject);
  #else
  _collisionWorld->addCollisionObject(&*btobject);
  #endif

  // create a record to keep track of things
  // (for static contactor, ds=nil)
  std11::shared_ptr<BR> record(
    std11::make_shared<BR>(base, ds, shape, btshape, btobject, contactor));

  bodyShapeMap[ds ? &*ds : 0].push_back(record);

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
  return btTransform( rbase * roffset,
                      btVector3(position(0), position(1), position(2)) + rboffset );
}

void SiconosBulletCollisionManager_impl::updateShapePosition(const BodyShapeRecord &record)
{
  SiconosVector q(7);
  if (record.base)
    q = *record.base;
  else {
    q.zero();
    q(3) = 1;
  }

  DEBUG_PRINTF("updating shape position: %p(%ld) - %f, %f, %f\n",
               &*box,box.use_count(), q(0), q(1), q(2));

  btTransform t = offsetTransform(q, *record.contactor->offset);
  t.setOrigin(t.getOrigin() * _options.worldScale);
  record.btobject->setWorldTransform( t );
}

void SiconosBulletCollisionManager_impl::createCollisionObject(
  const SP::SiconosVector base,
  const SP::BodyDS ds,
  SP::SiconosSphere sphere,
  const SP::SiconosContactor contactor)
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
  createCollisionObjectHelper<SP::SiconosSphere, SP::BTSPHERESHAPE, BodySphereRecord>
    (base, ds, sphere, btsphere, bodyShapeMap, contactor);
}

void SiconosBulletCollisionManager_impl::updateContactorInertia(
  SP::NewtonEulerDS ds, SP::btCollisionShape btshape)
{
  btVector3 localinertia;
  btshape->calculateLocalInertia(ds->scalarMass(), localinertia);
  assert( ! ((localinertia.x() == 0.0
              && localinertia.y() == 0.0
              && localinertia.z() == 0.0)
             || isinf(localinertia.x())
             || isinf(localinertia.y())
             || isinf(localinertia.z()))
          && "calculateLocalInertia() returned garbage" );
  ds->setInertia(localinertia[0], localinertia[1], localinertia[2]);
}

void SiconosBulletCollisionManager_impl::updateShape(BodySphereRecord &record)
{
  SP::SiconosSphere sphere(record.shape);
  SP::BTSPHERESHAPE btsphere(record.btshape);

  if (sphere->version() != record.shape_version)
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

    if (record.ds && record.ds->useContactorInertia())
      updateContactorInertia(record.ds, btsphere);

    if (record.btobject->getBroadphaseHandle())
    {
      _collisionWorld->updateSingleAabb(&*record.btobject);
      _collisionWorld->getBroadphase()->getOverlappingPairCache()->
        cleanProxyFromPairs(record.btobject->getBroadphaseHandle(), &*_dispatcher);
    }

    record.shape_version = sphere->version();
  }

  updateShapePosition(record);
}

void SiconosBulletCollisionManager_impl::createCollisionObject(
  const SP::SiconosVector base,
  const SP::BodyDS ds,
  SP::SiconosPlane plane,
  SP::SiconosContactor contactor)
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
  const btScalar pts[] = {
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
  createCollisionObjectHelper<SP::SiconosPlane, SP::BTPLANESHAPE, BodyPlaneRecord>
    (base, ds, plane, btplane, bodyShapeMap, contactor);
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
  record.btobject->setWorldTransform( t );
}

void SiconosBulletCollisionManager_impl::createCollisionObject(
  const SP::SiconosVector base,
  const SP::BodyDS ds,
  SP::SiconosBox box,
  SP::SiconosContactor contactor)
{
  const btScalar half = 0.5;

#ifdef USE_CONVEXHULL_FOR_BOX
  const btScalar pts[] = {
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
  createCollisionObjectHelper<SP::SiconosBox, SP::BTBOXSHAPE, BodyBoxRecord>
    (base, ds, box, btbox, bodyShapeMap, contactor);
}

void SiconosBulletCollisionManager_impl::updateShape(BodyBoxRecord &record)
{
  SP::SiconosBox box(record.shape);
  SP::BTBOXSHAPE btbox(record.btshape);

  // Update shape parameters
  if (box->version() != record.shape_version)
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

    if (record.ds && record.ds->useContactorInertia())
      updateContactorInertia(record.ds, btbox);

    if (record.btobject->getBroadphaseHandle())
    {
      _collisionWorld->updateSingleAabb(&*record.btobject);
      _collisionWorld->getBroadphase()->getOverlappingPairCache()->
        cleanProxyFromPairs(record.btobject->getBroadphaseHandle(), &*_dispatcher);
    }

    record.shape_version = box->version();
  }

  updateShapePosition(record);
}

void SiconosBulletCollisionManager_impl::createCollisionObject(
  const SP::SiconosVector base,
  const SP::BodyDS ds,
  SP::SiconosCylinder cylinder,
  SP::SiconosContactor contactor)
{
  SP::BTCYLSHAPE btcylinder(new BTCYLSHAPE(btVector3(1.0, 1.0, 1.0)));

  // initialization
  createCollisionObjectHelper<SP::SiconosCylinder, SP::BTCYLSHAPE, BodyCylinderRecord>
    (base, ds, cylinder, btcylinder, bodyShapeMap, contactor);
}

void SiconosBulletCollisionManager_impl::updateShape(BodyCylinderRecord &record)
{
  SP::SiconosCylinder cyl(record.shape);
  SP::BTCYLSHAPE btcyl(record.btshape);

  // Update shape parameters
  if (cyl->version() != record.shape_version)
  {
    // Bullet cylinder has an inside margin, so we add the outside
    // margin explicitly.
    double m = cyl->outsideMargin();

    double radius = (cyl->radius() + m) * _options.worldScale;
    double length = (cyl->length()/2 + m) * _options.worldScale;

    assert(radius > 0 && length > 0);

    btcyl->setLocalScaling(btVector3(radius, length, radius));
    btcyl->setMargin((cyl->insideMargin() + cyl->outsideMargin()) * _options.worldScale);

    if (record.ds && record.ds->useContactorInertia())
      updateContactorInertia(record.ds, btcyl);

    if (record.btobject->getBroadphaseHandle())
    {
      _collisionWorld->updateSingleAabb(&*record.btobject);
      _collisionWorld->getBroadphase()->getOverlappingPairCache()->
        cleanProxyFromPairs(record.btobject->getBroadphaseHandle(), &*_dispatcher);
    }

    record.shape_version = cyl->version();
  }

  updateShapePosition(record);
}

void SiconosBulletCollisionManager_impl::createCollisionObject(
  const SP::SiconosVector base,
  const SP::BodyDS ds,
  SP::SiconosConvexHull ch,
  SP::SiconosContactor contactor)
{
  if (!ch->vertices())
    throw SiconosException("No vertices matrix specified for convex hull.");

  if (ch->vertices()->size(1) != 3)
    throw SiconosException("Convex hull vertices matrix must have 3 columns.");

  // Copy and scale the points
  int rows = ch->vertices()->size(0);
  std::vector<btScalar> pts;
  pts.resize(rows*3);
  for (int r=0; r < rows; r++) {
    pts[r*3+0] = (*ch->vertices())(r, 0) * _options.worldScale;
    pts[r*3+1] = (*ch->vertices())(r, 1) * _options.worldScale;
    pts[r*3+2] = (*ch->vertices())(r, 2) * _options.worldScale;
  }

  SP::btConvexHullShape btch;
  if (ch->insideMargin() == 0)
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
    if (shrunkBy < 0) {
      // TODO: Warning
      // "insideMargin is too large, convex hull would be too small.";
      btch.reset(new btConvexHullShape(&pts[0], rows, sizeof(btScalar)*3));
      ch->setInsideMargin(0);
    }
    else {
      btch.reset(new btConvexHullShape);
      for (int i=0; i < shrinkCH.vertices.size(); i++) {
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
  createCollisionObjectHelper<SP::SiconosConvexHull, SP::BTCHSHAPE, BodyCHRecord>
    (base, ds, ch, btch, bodyShapeMap, contactor);
}

void SiconosBulletCollisionManager_impl::updateShape(BodyCHRecord &record)
{
  SP::SiconosConvexHull ch(record.shape);
  SP::BTCHSHAPE btch(record.btshape);

  // Update shape parameters
  if (ch->version() != record.shape_version)
  {
    // TODO
    //btbox->setLocalScaling(btVector3(sx, sy, sz));
    btch->setMargin((ch->insideMargin() + ch->outsideMargin()) * _options.worldScale);

    if (record.ds && record.ds->useContactorInertia())
      updateContactorInertia(record.ds, btch);

    if (record.btobject->getBroadphaseHandle())
    {
      _collisionWorld->updateSingleAabb(&*record.btobject);
      _collisionWorld->getBroadphase()->getOverlappingPairCache()->
        cleanProxyFromPairs(record.btobject->getBroadphaseHandle(), &*_dispatcher);
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
    std11::make_shared<btTriangleIndexVertexArray>(
      mesh->indexes()->size()/3,
      (int*)mesh->indexes()->data(),
      sizeof(int)*3,
      mesh->vertices()->size(1),
      mesh->vertices()->getArray(),
      sizeof(btScalar)*3));

  return std::make_pair(bttris, (btScalar*)NULL);
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
  for (unsigned int i=0; i < numVertices; i++)
  {
    vertices[i*3+0] = (*mesh->vertices())(0,i);
    vertices[i*3+1] = (*mesh->vertices())(1,i);
    vertices[i*3+2] = (*mesh->vertices())(2,i);
  }
  SP::btTriangleIndexVertexArray bttris(
    std11::make_shared<btTriangleIndexVertexArray>(
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
  const SP::BodyDS ds,
  SP::SiconosMesh mesh,
  SP::SiconosContactor contactor)
{
  if (!mesh->indexes())
    throw SiconosException("No indexes matrix specified for mesh.");

  if ((mesh->indexes()->size() % 3) != 0)
    throw SiconosException("Mesh indexes size must be divisible by 3.");

  if (!mesh->vertices())
    throw SiconosException("No vertices matrix specified for mesh.");

  if (mesh->vertices()->size(0) != 3)
    throw SiconosException("Convex hull vertices matrix must have 3 columns.");

  // Create Bullet triangle list, either by copying on non-copying method
  // TODO: worldScale on vertices
  std::pair<SP::btTriangleIndexVertexArray, btScalar*> datapair(
    make_bt_vertex_array(mesh, (btScalar)0, (*mesh->vertices())(0,0)));
  SP::btTriangleIndexVertexArray bttris(datapair.first);

  // Create Bullet mesh object
  SP::BTMESHSHAPE btmesh(std11::make_shared<BTMESHSHAPE>(&*bttris, true));

  // Hold on to the data since Bullet does not make a copy
  btmesh->btTriData = bttris;
  btmesh->btScalarVertices = datapair.second;

  // initialization
  createCollisionObjectHelper<SP::SiconosMesh, SP::BTMESHSHAPE, BodyMeshRecord>
    (base, ds, mesh, btmesh, bodyShapeMap, contactor);
}

void SiconosBulletCollisionManager_impl::updateShape(BodyMeshRecord &record)
{
  SP::SiconosMesh mesh(record.shape);
  SP::BTMESHSHAPE btmesh(record.btshape);

  // Update shape parameters
  if (mesh->version() != record.shape_version)
  {
    // btBvhTriangleMeshShape supports only outsideMargin.
    // TODO: support insideMargin, scale the points by their normals.
    btmesh->setMargin(mesh->outsideMargin() * _options.worldScale);
    btmesh->recalcLocalAabb();

    // TODO: Calculate inertia from a mesh.  For now we leave it at
    // identity, the user can provide an inertia if desired.

    // if (record.ds && record.ds->useContactorInertia())
    //   updateContactorInertia(record.ds, btmesh);

    if (record.btobject->getBroadphaseHandle())
    {
      _collisionWorld->updateSingleAabb(&*record.btobject);
      _collisionWorld->getBroadphase()->getOverlappingPairCache()->
        cleanProxyFromPairs(record.btobject->getBroadphaseHandle(), &*_dispatcher);
    }

    record.shape_version = mesh->version();
  }

  updateShapePosition(record);
}

void SiconosBulletCollisionManager_impl::createCollisionObject(
  const SP::SiconosVector base,
  const SP::BodyDS ds,
  SP::SiconosHeightMap heightmap,
  SP::SiconosContactor contactor)
{
  if (!heightmap->height_data())
    throw SiconosException("No height matrix specified for heightmap.");

  SP::SiconosMatrix data = heightmap->height_data();

  if (!data || data->size(0) < 2 || data->size(1) < 2)
    throw SiconosException("Height matrix does not have sufficient dimensions "
                           "to represent a plane.");

  // Create heightfield data for Bullet.  Make a copy in case data
  // type btScalar is different from SiconosMatrix.
  // Calculate min and max value at the same time.
  double vmin = std::numeric_limits<double>::infinity();
  double vmax = -vmin;
  std11::shared_ptr< std::vector<btScalar> > heightfield(
    std11::make_shared< std::vector<btScalar> >());

  heightfield->resize(data->size(0) * data->size(1));

  for (unsigned int i=0; i < data->size(0); i++) {
    for (unsigned int j=0; j < data->size(1); j++) {
      double v = data->getValue(i,j);
      (*heightfield)[j*data->size(0)+i] = v;
      if (v > vmax) vmax = v;
      if (v < vmin) vmin = v;
    }
  }

  // Create Bullet height object
  SP::BTHEIGHTSHAPE btheight(std11::make_shared<BTHEIGHTSHAPE>(
    data->size(0), heightfield, vmin, vmax));

  // initialization
  createCollisionObjectHelper<SP::SiconosHeightMap, SP::BTHEIGHTSHAPE, BodyHeightRecord>
    (base, ds, heightmap, btheight, bodyShapeMap, contactor);
}

void SiconosBulletCollisionManager_impl::updateShape(BodyHeightRecord &record)
{
  SP::SiconosHeightMap height(record.shape);
  SP::BTHEIGHTSHAPE btheight(record.btshape);

  // Update shape parameters
  if (height->version() != record.shape_version)
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

    if (record.btobject->getBroadphaseHandle())
    {
      _collisionWorld->updateSingleAabb(&*record.btobject);
      _collisionWorld->getBroadphase()->getOverlappingPairCache()->
        cleanProxyFromPairs(record.btobject->getBroadphaseHandle(), &*_dispatcher);
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
  if (record.base)
    q = *record.base;
  else {
    q.zero();
    q(3) = 1;
  }
  t = offsetTransform(q, o);

  DEBUG_PRINTF("updating shape position: %p(%ld) - %f, %f, %f\n",
               &*box,box.use_count(), q(0), q(1), q(2));

  t.setOrigin(t.getOrigin() * _options.worldScale);
  record.btobject->setWorldTransform( t );
}

class CreateCollisionObjectShapeVisitor : public SiconosVisitor
{
public:
  using SiconosVisitor::visit;
  SiconosBulletCollisionManager_impl &impl;
  const SP::BodyDS ds;
  const SP::SiconosVector base;
  SP::SiconosContactor contactor;

  CreateCollisionObjectShapeVisitor(SiconosBulletCollisionManager_impl &_impl,
                                    const SP::BodyDS _ds,
                                    const SP::SiconosVector _base)
    : impl(_impl), ds(_ds), base(_base) {}

  void visit(SP::SiconosPlane shape)
    { impl.createCollisionObject(base, ds, shape, contactor); }
  void visit(SP::SiconosSphere shape)
    { impl.createCollisionObject(base, ds, shape, contactor); }
  void visit(SP::SiconosBox shape)
    { impl.createCollisionObject(base, ds, shape, contactor); }
  void visit(SP::SiconosCylinder shape)
    { impl.createCollisionObject(base, ds, shape, contactor); }
  void visit(SP::SiconosConvexHull shape)
    { impl.createCollisionObject(base, ds, shape, contactor); }
  void visit(SP::SiconosMesh shape)
    { impl.createCollisionObject(base, ds, shape, contactor); }
  void visit(SP::SiconosHeightMap shape)
    { impl.createCollisionObject(base, ds, shape, contactor); }
};

void SiconosBulletCollisionManager_impl::createCollisionObjectsForBodyContactorSet(
  const SP::BodyDS ds,
  SP::SiconosVector base,
  SP::SiconosContactorSet contactors)
{
  // ensure consistency between ds and base and contactor -- if they
  // can be taken from the DS, ensure nothing else is provided.
  assert (!(ds && base) && "Provide either DS or base, but not both!");
  assert (!(ds && contactors) && "Provide either DS or contactor set, but not both!");

  // of course, we need at least one of the two combinations
  assert (ds || (base && contactors));

  SP::SiconosContactorSet con(contactors);
  if (ds) {
    con = ds->contactors();
    base = ds->q();
  }
  if (!con) return;

  std11::shared_ptr<CreateCollisionObjectShapeVisitor>
    ccosv(new CreateCollisionObjectShapeVisitor(*this, ds, base));

  /* Call createCollisionObject for each shape type using the visitor
   * defined above */
  std::vector< SP::SiconosContactor >::const_iterator it;
  for (it=con->begin(); it!=con->end(); it++)
  {
    ccosv->contactor = *it;
    ccosv->contactor->shape->acceptSP(ccosv);
  }
}

void SiconosBulletCollisionManager::removeBody(const SP::BodyDS& body)
{
  BodyShapeMap::iterator it( impl->bodyShapeMap.find(&*body) );
  if (it == impl->bodyShapeMap.end())
    return;

  std::vector<std11::shared_ptr<BodyShapeRecord> >::iterator it2;
  for (it2 = it->second.begin(); it2 != it->second.end(); it2++)
  {
    printf("removing a collision object\n");
    impl->_collisionWorld->removeCollisionObject( &* (*it2)->btobject );
  }

  impl->bodyShapeMap.erase(it);
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

  class iterator {
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
    const ContactPointTuple& operator*() {
      return data;
    };
    const ContactPointTuple* operator->() {
      return &data;
    };

    iterator()
      { numManifolds = 0;
        manifold_index = -1;
        contact_index = -1;
        numContacts = 0; }

    iterator& operator++() {
      if (numManifolds == 0)
        return *this;
      contact_index ++;
      while (contact_index >= numContacts)
      {
        manifold_index ++;
        if (manifold_index < numManifolds)
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

    bool operator!=(const iterator &it) {
      if (it.numManifolds==0) return numManifolds!=0;
      return data.objectA != it.data.objectA
        || data.objectB != it.data.objectB
        || data.point != it.data.point;
    };
    friend class IterateContactPoints;
  };

  iterator begin() {
    return iterator(world);
  };

  iterator end() {
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
  DEBUG_PRINTF("unlinking interaction %p\n", &**p_inter);
  std11::static_pointer_cast<BulletR>((*p_inter)->relation())->preDelete();
  gSimulation->unlink(*p_inter);
  delete p_inter;
  return false;
}

SP::BulletR SiconosBulletCollisionManager::makeBulletR(SP::BodyDS ds1,
                                                       SP::SiconosShape shape1,
                                                       SP::BodyDS ds2,
                                                       SP::SiconosShape shape2,
                                                       const btManifoldPoint &p,
                                                       bool flip,
                                                       double y_correction_A,
                                                       double y_correction_B,
                                                       double scaling)
{
  return std11::make_shared<BulletR>(p, ds1 ? ds1->q() : SP::SiconosVector(),
                                     ds2 ? ds2->q() : SP::SiconosVector(),
                                     flip, y_correction_A,
                                     y_correction_B, scaling);
}

class CollisionUpdateVisitor : public SiconosVisitor
{
public:
  using SiconosVisitor::visit;
  SiconosBulletCollisionManager_impl &impl;

  CollisionUpdateVisitor(SiconosBulletCollisionManager_impl& _impl)
    : impl(_impl) {}

  void visit(SP::BodyDS bds)
  {
    if (bds->contactors()) {
      if (impl.bodyShapeMap.find(&*bds) == impl.bodyShapeMap.end())
      {
        impl.createCollisionObjectsForBodyContactorSet(bds);
      }
      impl.updateAllShapesForDS(*bds);
    }
  }
};

void SiconosBulletCollisionManager::updateInteractions(SP::Simulation simulation)
{
  // -2. update collision objects from all BodyDS dynamical systems
  SP::SiconosVisitor updateVisitor(new CollisionUpdateVisitor(*impl));
  simulation->nonSmoothDynamicalSystem()->visitDynamicalSystems(updateVisitor);

  // Clear cache automatically before collision detection if requested
  if (_options.clearOverlappingPairCache)
    clearOverlappingPairCache();

  if (! impl->_queuedCollisionObjects.empty())
  {

    std::vector<SP::btCollisionObject>::iterator it;
    for (it = impl->_queuedCollisionObjects.begin();
         it != impl->_queuedCollisionObjects.end();
         ++ it)
    {
      impl->_collisionWorld->addCollisionObject(&**it);
    }
    impl->_queuedCollisionObjects.clear();
  }

  // -1. reset statistical counters
  resetStatistics();

  // 0. set up bullet callbacks
  gSimulation = &*simulation;
  gContactDestroyedCallback = this->bulletContactClear;

  // Important parameter controlling contact point making and breaking
  gContactBreakingThreshold = _options.contactBreakingThreshold;

  // 1. perform bullet collision detection
  impl->_collisionWorld->performDiscreteCollisionDetection();

  // 2. deleted contact points have been removed from the graph during the
  //    bullet collision detection callbacks

  // 3. for each contact point, if there is no interaction, create one
  IterateContactPoints t(impl->_collisionWorld);
  IterateContactPoints::iterator it, itend=t.end();
  DEBUG_PRINT("iterating contact points:\n");

  for (it=t.begin(); it!=itend; ++it)
  {
    DEBUG_PRINTF("  -- %p, %p, %p\n", it->objectA, it->objectB, it->point);

    // Get the BodyDS and SiconosShape pointers

    const BodyShapeRecord *pairA, *pairB;
    pairA = reinterpret_cast<const BodyShapeRecord*>(it->objectA->getUserPointer());
    pairB = reinterpret_cast<const BodyShapeRecord*>(it->objectB->getUserPointer());
    assert(pairA && pairB && "btCollisionObject had a null user pointer!");

    // The first pair will always be the non-static object
    bool flip = false;
    if (pairB->ds && !pairA->ds) {
      pairA = reinterpret_cast<const BodyShapeRecord*>(it->objectB->getUserPointer());
      pairB = reinterpret_cast<const BodyShapeRecord*>(it->objectA->getUserPointer());
      flip = true;
    }

    // If both collision objects belong to the same body (or no body),
    // no interaction is created.
    if (pairA->ds == pairB->ds)
      continue;

    // If the two bodies are already connected by another type of
    // relation (e.g. EqualityCondition == they have a joint between
    // them), then don't create contact constraints, because it leads
    // to an ill-conditioned problem.
    if (pairA->ds && pairB->ds)
    {
      InteractionsGraph::VIterator ui, uiend;
      SP::InteractionsGraph indexSet0 = simulation->nonSmoothDynamicalSystem()->topology()->indexSet0();
      bool match = false;
      for (std11::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
      {
        SP::Interaction inter( indexSet0->bundle(*ui) );
        SP::BodyDS ds1( std11::dynamic_pointer_cast<BodyDS>(
                          indexSet0->properties(*ui).source) );
        SP::BodyDS ds2( std11::dynamic_pointer_cast<BodyDS>(
                          indexSet0->properties(*ui).target) );
        if (ds1 && ds2 && (((&*ds1==&*pairA->ds) && (&*ds2==&*pairB->ds))
                           || ((&*ds1==&*pairB->ds) && (&*ds2==&*pairA->ds))))
        {
          // Only match on non-BulletR interactions, i.e. non-contact relations
          SP::BulletR br ( std11::dynamic_pointer_cast<BulletR>(inter->relation()) );
          if (!br) {
            SP::NewtonEulerJointR jr (
              std11::dynamic_pointer_cast<NewtonEulerJointR>(inter->relation()) );

            /* If it is a joint, check the joint self-collide property */
            if (jr && !jr->allowSelfCollide())
              match = true;

            /* If any non-contact relation is found, both bodies must
             * allow self-collide */
            if (!pairA->ds->allowSelfCollide() || !pairB->ds->allowSelfCollide())
              match = true;
          }
          if (match) break;
        }
      }
      if (match)
        continue;
    }

    if (it->point->m_userPersistentData)
    {
      /* interaction already exists */
      SP::Interaction *p_inter =
        (SP::Interaction*)it->point->m_userPersistentData;

      /* update the relation */
      SP::BulletR rel(std11::static_pointer_cast<BulletR>((*p_inter)->relation()));
      rel->updateContactPoints(*it->point, pairA->ds,
                               pairB->ds ? pairB->ds : SP::NewtonEulerDS());

      _stats.existing_interactions_processed ++;
    }
    else
    {
      /* new interaction */
      SP::Interaction inter;

      int g1 = pairA->contactor->collision_group;
      int g2 = pairB->contactor->collision_group;
      SP::NonSmoothLaw nslaw = nonSmoothLaw(g1,g2);

      if (nslaw && nslaw->size() == 3)
      {
        // Remove the added outside margin as a correction factor in Relation
        double combined_margin =
          pairA->sshape->outsideMargin() + pairB->sshape->outsideMargin();

        SP::BulletR rel(makeBulletR(pairA->ds, pairA->sshape,
                                    pairB->ds, pairB->sshape,
                                    *it->point,
                                    flip,
                                    pairA->sshape->outsideMargin(),
                                    pairB->sshape->outsideMargin(),
                                    1.0 / _options.worldScale));

        if (!rel) continue;
        rel->updateContactPoints(*it->point,
                                 pairA->ds ? pairA->ds : SP::NewtonEulerDS(),
                                 pairB->ds ? pairB->ds : SP::NewtonEulerDS());

        // We wish to be sure that no Interactions are created without
        // sufficient warning before contact.  TODO: Replace with exception or
        // flag.
        if ((it->point->getDistance() + combined_margin) < 0.0) {
          DEBUG_PRINTF(stderr, "Interactions must be created with positive distance (%f).\n",
                       (it->point->getDistance() + combined_margin)/_options.worldScale);
          _stats.interaction_warnings ++;
        }

        inter.reset(new Interaction(nslaw, rel/*4 * i + z*/));
        _stats.new_interactions_created ++;
      }
      else
      {
        if (nslaw && nslaw->size() == 1)
        {
          SP::BulletFrom1DLocalFrameR rel(
            new BulletFrom1DLocalFrameR(createSPtrbtManifoldPoint(*it->point)));
          inter.reset(new Interaction(nslaw, rel /*4 * i + z*/));
        }
      }

      if (inter)
      {
        /* store interaction in the contact point data, it will be freed by the
         * Bullet callback gContactDestroyedCallback */
        /* note: storing pointer to shared_ptr! */
        it->point->m_userPersistentData = (void*)(new SP::Interaction(inter));

        /* link bodies by the new interaction */
        simulation->link(inter, pairA->ds, pairB->ds);
      }
    }
  }
}

void SiconosBulletCollisionManager::clearOverlappingPairCache()
{
  if (!impl->_collisionWorld) return;

  BodyShapeMap::iterator it;
  btOverlappingPairCache *pairCache =
    impl->_collisionWorld->getBroadphase()->getOverlappingPairCache();

  for (it = impl->bodyShapeMap.begin(); it != impl->bodyShapeMap.end(); it++)
  {
    std::vector< std11::shared_ptr<BodyShapeRecord> >::iterator rec;
    for (rec = it->second.begin(); rec != it->second.end(); rec++) {
      if ((*rec)->btobject) {
        pairCache-> cleanProxyFromPairs((*rec)->btobject->getBroadphaseHandle(),
                                        &*impl->_dispatcher);
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
  if (closestOnly)
  {
    btCollisionWorld::ClosestRayResultCallback rayResult(btstart, btend);
    impl->_collisionWorld->rayTest(btstart, btend, rayResult);

    if (rayResult.hasHit())
    {
      const BodyShapeRecord *rec =
        reinterpret_cast<const BodyShapeRecord*>(
          rayResult.m_collisionObject->getUserPointer());

      if (rec) {
        SP::SiconosCollisionQueryResult result(
          std11::make_shared<SiconosCollisionQueryResult>());
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
      else {
        DEBUG_PRINTF("Siconos warning: BodyShapeRecord found by intersection was null.");
      }
    }
  }

  // Return more than one object, Bullet provides a different
  // interface
  else
  {
    btCollisionWorld::AllHitsRayResultCallback rayResult(btstart, btend);
    impl->_collisionWorld->rayTest(btstart, btend, rayResult);

    if (rayResult.hasHit())
    {
      for (int i=0; i < rayResult.m_collisionObjects.size(); i++)
      {
        const BodyShapeRecord *rec =
          reinterpret_cast<const BodyShapeRecord*>(
            rayResult.m_collisionObjects[i]->getUserPointer());

        if (rec) {
          SP::SiconosCollisionQueryResult result(
            std11::make_shared<SiconosCollisionQueryResult>());
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
        else {
          DEBUG_PRINTF("Siconos warning: BodyShapeRecord found by intersection was null.");
        }
      }
    }
  }

  if (sorted && result_list.size() > 1)
    std::sort(result_list.begin(), result_list.end(), cmpQueryResult);

  return result_list;
}
