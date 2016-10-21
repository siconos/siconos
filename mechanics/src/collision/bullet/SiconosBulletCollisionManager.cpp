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

/*! \file SiconosBulletCollisionManager.cpp
  \brief Definition of a Bullet-based interaction handler for contact
  detection.
*/

#include <MechanicsFwd.hpp>

#include "BulletSiconosFwd.hpp"
#include "SiconosBulletCollisionManager.hpp"
#include "BodyDS.hpp"
#include "BulletR.hpp"
#include "BulletFrom1DLocalFrameR.hpp"

#include <map>
#include <boost/format.hpp>
#include <boost/make_shared.hpp>

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

#include <Question.hpp>

#include <BulletCollision/CollisionDispatch/btCollisionWorld.h>
#include <BulletCollision/Gimpact/btGImpactCollisionAlgorithm.h>
#include <BulletCollision/CollisionDispatch/btDefaultCollisionConfiguration.h>
#include <BulletCollision/BroadphaseCollision/btDbvtBroadphase.h>
#include <BulletCollision/BroadphaseCollision/btAxisSweep3.h>

#include <BulletCollision/CollisionShapes/btStaticPlaneShape.h>
#include <BulletCollision/CollisionShapes/btSphereShape.h>
#include <BulletCollision/CollisionShapes/btBoxShape.h>
#include <BulletCollision/CollisionShapes/btConvexHullShape.h>
#include <LinearMath/btConvexHullComputer.h>

#include <LinearMath/btQuaternion.h>
#include <LinearMath/btVector3.h>

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

// We need to maintain a record associating each body with a shape,
// contactor, and collision object for each shape type.  We also need
// to access generic shape stuff (group, margin) by a pointer from the
// collision callback, so we need a record base class.
struct BodyShapeRecord
{
  BodyShapeRecord(SP::SiconosVector b, SP::BodyDS d, SP::SiconosShape sh,
                  SP::btCollisionObject btobj, SP::SiconosContactor con)
    : base(b), ds(d), sshape(sh), btobject(btobj), contactor(con),
      shape_version(sh->version()) {}
  SP::SiconosVector base;
  SP::BodyDS ds;
  SP::SiconosShape sshape;
  SP::btCollisionObject btobject;
  SP::SiconosContactor contactor;
  unsigned int shape_version;
};

template <typename SICONOSSHAPE, typename BULLETSHAPE>
struct BodyShapeRecordT : BodyShapeRecord
{
  BodyShapeRecordT(SP::SiconosVector base, SP::BodyDS ds,
                   SICONOSSHAPE sh, BULLETSHAPE btsh,
                   SP::btCollisionObject btobj, SP::SiconosContactor con)
    : BodyShapeRecord(base, ds, sh, btobj, con), shape(sh), btshape(btsh) {}
  SICONOSSHAPE shape;
  BULLETSHAPE btshape;
};

// Body-Shape map types
typedef BodyShapeRecordT<SP::SiconosBox, SP::BTBOXSHAPE> BodyBoxRecord;
typedef std::map<const BodyDS*, std::vector<std11::shared_ptr<BodyBoxRecord> > >
  BodyBoxMap;

typedef BodyShapeRecordT<SP::SiconosSphere, SP::BTSPHERESHAPE> BodySphereRecord;
typedef std::map<const BodyDS*, std::vector<std11::shared_ptr<BodySphereRecord> > >
  BodySphereMap;

typedef BodyShapeRecordT<SP::SiconosConvexHull, SP::BTCHSHAPE> BodyCHRecord;
typedef std::map<const BodyDS*, std::vector<std11::shared_ptr<BodyCHRecord> > >
  BodyCHMap;

typedef BodyShapeRecordT<SP::SiconosPlane, SP::BTPLANESHAPE> BodyPlaneRecord;
typedef std::map<const BodyDS*, std::vector<std11::shared_ptr<BodyPlaneRecord> > >
  BodyPlaneMap;

/** For associating static contactor sets and their offsets.
 * Pointer to this is what is returned as the opaque and unique
 * StaticContactorSetID so that they can be removed. */
struct StaticContactorSetRecord
{
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
   * btCollisionObject, so need a map DS->btXShape. We don't use an
   * SP::BodyDS because we need to use it from a const visitor. */
  BodyBoxMap bodyBoxMap;
  BodyCHMap bodyCHMap;
  BodySphereMap bodySphereMap;
  BodyPlaneMap bodyPlaneMap;

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
                             const SP::SiconosConvexHull ch,
                             const SP::SiconosContactor contactor);

  /* Call the above functions for each shape associated with a body or contactor. */
  void createCollisionObjectsForBodyContactorSet(
    const SP::BodyDS ds,
    const SP::SiconosVector base = SP::SiconosVector(),
    const SP::SiconosContactorSet contactor = SP::SiconosContactorSet());

  /* A helper function used to initialise new shapes, generic to the
   * shape type */
  template<typename ST, typename BT, typename BR, typename BSM>
  void createCollisionObjectHelper(SP::SiconosVector base, SP::BodyDS ds,
                                   ST& shape, BT& btshape, BSM& bodyShapeMap,
                                   SP::SiconosContactor contactor);

  void updateShape(BodyBoxRecord &record);
  void updateShape(BodySphereRecord &record);
  void updateShape(BodyCHRecord &record);
  void updateShape(BodyPlaneRecord &record);

  void updateAllShapesForDS(const BodyDS &bds);
  void updateShapePosition(const BodyShapeRecord &record);

  /* Helper to apply an offset transform to a position and return as a
   * btTransform */
  btTransform offsetTransform(const SiconosVector& position,
                              const SiconosVector& offset);

  SiconosBulletOptions &_options;

  std::vector<SP::btCollisionObject> _queuedCollisionObjects;

public:
  SiconosBulletCollisionManager_impl(SiconosBulletOptions &op) : _options(op) {}
  ~SiconosBulletCollisionManager_impl() {}

  friend class SiconosBulletCollisionManager;
  friend class CollisionUpdateVisitor;
  friend class CreateCollisionObjectShapeVisitor;
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

  if (_options.useMultipointIterations)
  {
    impl->_collisionConfiguration->setConvexConvexMultipointIterations();
    impl->_collisionConfiguration->setPlaneConvexMultipointIterations();
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
  impl->_collisionWorld->getDispatchInfo().m_convexConservativeDistanceThreshold = 0.1f;
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

void SiconosBulletCollisionManager_impl::updateAllShapesForDS(const BodyDS &bds)
{
  std::vector<std11::shared_ptr<BodyPlaneRecord> >::iterator itp;
  for (itp = bodyPlaneMap[&bds].begin(); itp != bodyPlaneMap[&bds].end(); itp++)
    updateShape(**itp);

  std::vector<std11::shared_ptr<BodySphereRecord> >::iterator its;
  for (its = bodySphereMap[&bds].begin(); its != bodySphereMap[&bds].end(); its++)
    updateShape(**its);

  std::vector<std11::shared_ptr<BodyBoxRecord> >::iterator itb;
  for (itb = bodyBoxMap[&bds].begin(); itb != bodyBoxMap[&bds].end(); itb++)
    updateShape(**itb);

  std::vector<std11::shared_ptr<BodyCHRecord> >::iterator itc;
  for (itc = bodyCHMap[&bds].begin(); itc != bodyCHMap[&bds].end(); itc++)
    updateShape(**itc);
}

template<typename ST, typename BT, typename BR, typename BSM>
void SiconosBulletCollisionManager_impl::createCollisionObjectHelper(
  SP::SiconosVector base, SP::BodyDS ds, ST& shape, BT& btshape,
  BSM& bodyShapeMap, SP::SiconosContactor contactor)
{
  assert(base && "Collision objects must have a base position.");

  // create corresponding Bullet object and shape
  SP::btCollisionObject btobject(new btCollisionObject());

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
  btQuaternion r(rbase * roffset);

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

  record.btobject->setWorldTransform( offsetTransform(q, *record.contactor->offset) );
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
  createCollisionObjectHelper<SP::SiconosSphere, SP::BTSPHERESHAPE,
                              BodySphereRecord, BodySphereMap>
    (base, ds, sphere, btsphere, bodySphereMap, contactor);
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
  SP::btBoxShape btplane(
    new btBoxShape(btVector3(1000*_options.worldScale,
                             1000*_options.worldScale,
                             1000*_options.worldScale)));
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
  createCollisionObjectHelper<SP::SiconosPlane, SP::BTPLANESHAPE,
                              BodyPlaneRecord, BodyPlaneMap>
    (base, ds, plane, btplane, bodyPlaneMap, contactor);
}

void SiconosBulletCollisionManager_impl::updateShape(BodyPlaneRecord& record)
{
  SP::SiconosPlane plane(record.shape);
  SP::BTPLANESHAPE btplane(record.btshape);

  // Update object parameters
  SiconosVector &q = *record.base;
#ifdef USE_BOX_FOR_PLANE
  btTransform tr(btQuaternion(q(4), q(5), q(6), q(3)),
    btVector3(q(0)*_options.worldScale,
              q(1)*_options.worldScale,
              (q(2) + plane->outsideMargin() - 1000)*_options.worldScale));
#else
#ifdef USE_CONVEXHULL_FOR_PLANE
  btTransform tr(btQuaternion(q(4), q(5), q(6), q(3)),
    btVector3(q(0)*_options.worldScale,
              q(1)*_options.worldScale,
              (q(2) - plane->insideMargin())*_options.worldScale));
#else
  btTransform tr(btQuaternion(q(4), q(5), q(6), q(3)),
    btVector3(q(0)*_options.worldScale,
              q(1)*_options.worldScale,
              (q(2) - plane->insideMargin())*_options.worldScale));
#endif
#endif

  // Note, we do not use generic updateShapePosition for plane
}

void SiconosBulletCollisionManager_impl::createCollisionObject(
  const SP::SiconosVector base,
  const SP::BodyDS ds,
  SP::SiconosBox box,
  SP::SiconosContactor contactor)
{
  const double half = 0.5;

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
  createCollisionObjectHelper<SP::SiconosBox, SP::BTBOXSHAPE,
                              BodyBoxRecord, BodyBoxMap>
    (base, ds, box, btbox, bodyBoxMap, contactor);
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
        btch->addPoint(v, false);
      }
      ch->setInsideMargin(shrunkBy / _options.worldScale);
    }
  }

  // Add external margin and recalc bounding box
  btch->setMargin((ch->insideMargin() + ch->outsideMargin()) * _options.worldScale);
  btch->recalcLocalAabb();

  // initialization
  createCollisionObjectHelper<SP::SiconosConvexHull, SP::BTCHSHAPE,
                              BodyCHRecord, BodyCHMap>
    (base, ds, ch, btch, bodyCHMap, contactor);
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

struct CreateCollisionObjectShapeVisitor : public SiconosVisitor
{
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
  void visit(SP::SiconosConvexHull shape)
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
  };

  class iterator {
  protected:
    SP::btCollisionWorld world;
    ContactPointTuple data;
    unsigned int numManifolds;
    unsigned int numContacts;
    unsigned int manifold_index;
    unsigned int contact_index;
    btPersistentManifold* contactManifold;
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
          contactManifold = world->getDispatcher()->
            getManifoldByIndexInternal(manifold_index);
          data.objectA = contactManifold->getBody0();
          data.objectB = contactManifold->getBody1();
          numContacts = contactManifold->getNumContacts();
          contact_index = 0;
        }
        else
        {
          numManifolds = 0;
          return *this;
        }
      }
      data.point = &(contactManifold->getContactPoint(contact_index));
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
  gSimulation->unlink(*p_inter);
  delete p_inter;
  return false;
}

struct CollisionUpdateVisitor : public SiconosVisitor
{
  using SiconosVisitor::visit;
  SiconosBulletCollisionManager_impl &impl;

  CollisionUpdateVisitor(SiconosBulletCollisionManager_impl& _impl)
    : impl(_impl) {}

  void visit(SP::BodyDS bds)
  {
    if (bds->contactors()) {
      BodyBoxMap::iterator it = impl.bodyBoxMap.find(&*bds);
      if (impl.bodyBoxMap.find(&*bds) == impl.bodyBoxMap.end()
          && impl.bodyCHMap.find(&*bds) == impl.bodyCHMap.end()
          && impl.bodySphereMap.find(&*bds) == impl.bodySphereMap.end()
          && impl.bodyPlaneMap.find(&*bds) == impl.bodyPlaneMap.end())
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
  gContactBreakingThreshold = _options.breakingThreshold;

  // 1. perform bullet collision detection
  impl->_collisionWorld->performDiscreteCollisionDetection();

  // 2. deleted contact points have been removed from the graph during the
  //    bullet collision detection callbacks

  // 3. for each contact point, if there is no interaction, create one
  IterateContactPoints t(impl->_collisionWorld);
  IterateContactPoints::iterator it, itend=t.end();
  DEBUG_PRINT("iterating contact points:\n");
  int n_points=0;
  int late_interaction=0;
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

    if (it->point->m_userPersistentData)
    {
      // do what's needed for an interaction already present
      SP::Interaction *p_inter =
        (SP::Interaction*)it->point->m_userPersistentData;
      // (note: nothing for now!)
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

        SP::BulletR rel(new BulletR(createSPtrbtManifoldPoint(*it->point),
                                    flip,
                                    pairA->sshape->outsideMargin(),
                                    pairB->sshape->outsideMargin(),
                                    1.0 / _options.worldScale));

        // We wish to be sure that no Interactions are created without
        // sufficient warning before contact.  TODO: Replace with exception or
        // flag.
        if ((it->point->getDistance() + combined_margin) < 0.0) {
          DEBUG_PRINTF(stderr, "Interactions must be created with positive distance (%f).\n",
                       (it->point->getDistance() + combined_margin)/_options.worldScale);
          _stats.interaction_warnings ++;
        }

        inter.reset(new Interaction(3, nslaw, rel, 0 /*4 * i + z*/));
        _stats.new_interactions_created ++;
      }
      else
      {
        if (nslaw && nslaw->size() == 1)
        {
          SP::BulletFrom1DLocalFrameR rel(
            new BulletFrom1DLocalFrameR(createSPtrbtManifoldPoint(*it->point)));
          inter.reset(new Interaction(1, nslaw, rel, 0 /*4 * i + z*/));
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
