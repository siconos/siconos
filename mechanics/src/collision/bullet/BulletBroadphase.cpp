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

#include <MechanicsFwd.hpp>

#include "BulletSiconosFwd.hpp"
#include "BulletBroadphase.hpp"
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

// We need to maintain a 3-way record associating each body shape pair
// for each shape type, and maintaining other memory we need too.  We
// also need to access generic shape stuff (group, margin) by a
// pointer, so we need a base class.
struct BodyShapeRecord
{
  BodyShapeRecord(SP::BodyDS d, SP::SiconosShape sh,
                  SP::btCollisionObject btobj, SP::SiconosVector off)
    : ds(d), sshape(sh), btobject(btobj), offset(off) {}
  SP::BodyDS ds;
  SP::SiconosShape sshape;
  SP::btCollisionObject btobject;
  SP::SiconosVector offset;
};

template <typename SICONOSSHAPE, typename BULLETSHAPE>
struct BodyShapeRecordT : BodyShapeRecord
{
  BodyShapeRecordT(SP::BodyDS d, SICONOSSHAPE sh, BULLETSHAPE btsh,
                   SP::btCollisionObject btobj, SP::SiconosVector off)
    : BodyShapeRecord(d, sh, btobj, off), shape(sh), btshape(btsh) {}
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

class CollisionUpdater;

class BulletBroadphase_impl
{
protected:
  SP::btCollisionWorld _collisionWorld;
  SP::btDefaultCollisionConfiguration _collisionConfiguration;
  SP::btCollisionDispatcher _dispatcher;
  SP::btBroadphaseInterface _broadphase;

  std::vector<SP::SiconosContactor> staticContactors;

  /* Work-around: Adding static contactors before adding the initial
   * bodies changes the contact points that are generated by Bullet!
   * Therefore we need to "queue" these contactors and add them
   * *after* the initial bodies have been added.  It's a bit ugly,
   * but seems to be the only way to maintain same order of calls to
   * addCollisionObject(). */
  std::vector<SP::SiconosContactor> queuedStaticContactors;
  void insertQueuedContactors(BulletBroadphase &broad);

  // Non-smooth laws
  std::map<std::pair<int,int>, SP::NonSmoothLaw> nslaws;

  /* During iteration over DSs for position updates we need to access
   * btCollisionObject, so need a map DS->btXShape. We don't use an
   * SP::BodyDS because we need to use it from a const visitor. */
  BodyBoxMap bodyBoxMap;
  BodyCHMap bodyCHMap;
  BodySphereMap bodySphereMap;
  BodyPlaneMap bodyPlaneMap;

  SP::Simulation _simulation;

  /* Create collision objects for each shape type */
  void createCollisionObject(const SP::BodyDS ds,
                             const SP::SiconosPlane plane,
                             const SP::SiconosVector offset);
  void createCollisionObject(const SP::BodyDS ds,
                             const SP::SiconosSphere sphere,
                             const SP::SiconosVector offset);
  void createCollisionObject(const SP::BodyDS ds,
                             const SP::SiconosBox box,
                             const SP::SiconosVector offset);
  void createCollisionObject(const SP::BodyDS ds,
                             const SP::SiconosConvexHull ch,
                             const SP::SiconosVector offset);

  /* Call the above functions for each shape associated with a body or contactor. */
  void createCollisionObjectsForBodyContactor(
    const SP::BodyDS ds, const SP::SiconosContactor contactor = SP::SiconosContactor());

  /* A helper function used to initialise new shapes, generic to the
   * shape type */
  template<typename ST, typename BT, typename BR, typename BSM>
  void createCollisionObjectHelper(SP::BodyDS ds, ST& shape, BT& btshape,
                                   BSM& bodyShapeMap, SP::SiconosVector offset);

  void updateShape(const BodyBoxRecord &record);
  void updateShape(const BodySphereRecord &record);
  void updateShape(const BodyCHRecord &record);
  void updateShape(const BodyPlaneRecord &record);

  void updateAllShapesForDS(const BodyDS &bds);
  void updateShapePosition(const BodyShapeRecord &record);

  /* Helper to apply an offset transform to a position and return as a
   * btTransform */
  btTransform offsetTransform(const SiconosVector& position,
                              const SiconosVector& offset);

  BulletOptions &_options;

public:
  BulletBroadphase_impl(BulletOptions &op) : _options(op) {}
  ~BulletBroadphase_impl() {}

  friend class BulletBroadphase;
  friend class CollisionUpdater;
};

void BulletBroadphase_impl::insertQueuedContactors(BulletBroadphase &broad)
{
  std::vector<SP::SiconosContactor>::iterator con;
  for (con = queuedStaticContactors.begin();
       con != queuedStaticContactors.end();
       con ++)
  {
    createCollisionObjectsForBodyContactor(NULL, *con);
    staticContactors.push_back(*con);
  }
  queuedStaticContactors.clear();
}

void BulletBroadphase::initialize_impl()
{
  impl.reset(new BulletBroadphase_impl(_options));
  impl->_collisionConfiguration.reset(
    new btDefaultCollisionConfiguration());

  impl->_collisionConfiguration->setConvexConvexMultipointIterations();
  impl->_collisionConfiguration->setPlaneConvexMultipointIterations();

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

BulletBroadphase::BulletBroadphase()
{
  initialize_impl();
}

BulletBroadphase::BulletBroadphase(const BulletOptions &options)
  : _options(options)
{
  initialize_impl();
}

BulletBroadphase::~BulletBroadphase()
{
  // unlink() will be called on all remaining
  // contact points when world is destroyed
  gBulletBroadphase = this;

  // must be the first de-allocated, otherwise segfault
  impl->_collisionWorld.reset();
}

void BulletBroadphase::insertStaticContactor(SP::SiconosContactor contactor)
{
  /* Work-around: Instead of adding them directly, we have to queue
   * them to be added after the initial bodies, otherwise Bullet
   * behaves differently. */
  impl->queuedStaticContactors.push_back(contactor);
}

void BulletBroadphase_impl::updateAllShapesForDS(const BodyDS &bds)
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
void BulletBroadphase_impl::createCollisionObjectHelper(
  SP::BodyDS ds, ST& shape, BT& btshape, BSM& bodyShapeMap, SP::SiconosVector offset)
{
  // create corresponding Bullet object and shape
  SP::btCollisionObject btobject(new btCollisionObject());

  // associate the shape with the object
  btobject->setCollisionShape(&*btshape);

  // put it in the world
  _collisionWorld->addCollisionObject(&*btobject);

  // create a record to keep track of things
  // (for static contactor, ds=nil)
  std11::shared_ptr<BR> record(
    std11::make_shared<BR>(ds, shape, btshape, btobject, offset));

  bodyShapeMap[ds ? &*ds : 0].push_back(record);

  assert(record->btobject);
  assert(record->sshape);
  assert(record->shape);
  assert(record->btshape);
  assert(record->offset);
  assert(record->offset->size() == 7);

  // Allow Bullet to report colliding DSs.  We need to access it from
  // the collision callback as the record base class so down-cast it.
  btobject->setUserPointer(
    reinterpret_cast<void*>(
      static_cast<BodyShapeRecord*>(&*record)));

  // initial parameter update
  updateShape(*record);
}

btTransform BulletBroadphase_impl::offsetTransform(const SiconosVector& position,
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

void BulletBroadphase_impl::updateShapePosition(const BodyShapeRecord &record)
{
  SiconosVector q(7);
  if (record.ds)
    q = *record.ds->q();
  else {
    q.zero();
    q(3) = 1;
  }

  DEBUG_PRINTF("updating shape position: %p(%ld) - %f, %f, %f\n",
               &*box,box.use_count(), q(0), q(1), q(2));

  record.btobject->setWorldTransform( offsetTransform(q, *record.offset) );
}

void BulletBroadphase_impl::createCollisionObject(const SP::BodyDS ds,
                                                  SP::SiconosSphere sphere,
                                                  SP::SiconosVector offset)
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
    (ds, sphere, btsphere, bodySphereMap, offset);
}

void BulletBroadphase_impl::updateShape(const BodySphereRecord &record)
{
  SP::SiconosSphere sphere(record.shape);
  SP::BTSPHERESHAPE btsphere(record.btshape);

  double r = (sphere->radius() + sphere->outsideMargin()) * _options.worldScale;

  // Update shape parameters
#ifdef USE_CONVEXHULL_FOR_SPHERE
  btsphere->setMargin(r);
#else
  btsphere->setLocalScaling(btVector3(r, r, r));

  // btSphereShape has an internal margin
  btsphere->setMargin(sphere->insideMargin() * _options.worldScale);
#endif

  updateShapePosition(record);
}

void BulletBroadphase_impl::createCollisionObject(const SP::BodyDS ds,
                                                  SP::SiconosPlane plane,
                                                  SP::SiconosVector offset)
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
    (ds, plane, btplane, bodyPlaneMap, offset);
}

void BulletBroadphase_impl::updateShape(const BodyPlaneRecord& record)
{
  SP::SiconosPlane plane(record.shape);
  SP::BTPLANESHAPE btplane(record.btshape);

  // Update object parameters
  SiconosVector &q = *record.ds->q();
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

void BulletBroadphase_impl::createCollisionObject(const SP::BodyDS ds,
                                                  SP::SiconosBox box,
                                                  SP::SiconosVector offset)
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
    (ds, box, btbox, bodyBoxMap, offset);
}

void BulletBroadphase_impl::updateShape(const BodyBoxRecord &record)
{
  SP::SiconosBox box(record.shape);
  SP::BTBOXSHAPE btbox(record.btshape);

  // Update shape parameters
#ifdef USE_CONVEXHULL_FOR_BOX
  double m = -box->insideMargin();
#else
  double m = box->outsideMargin();
#endif

  double sx = ((*box->dimensions())(0) + m*2) * _options.worldScale;
  double sy = ((*box->dimensions())(1) + m*2) * _options.worldScale;
  double sz = ((*box->dimensions())(2) + m*2) * _options.worldScale;

  btbox->setLocalScaling(btVector3(sx, sy, sz));
  btbox->setMargin((box->insideMargin() + box->outsideMargin()) * _options.worldScale);

  updateShapePosition(record);
}

void BulletBroadphase_impl::createCollisionObject(const SP::BodyDS ds,
                                                  SP::SiconosConvexHull ch,
                                                  SP::SiconosVector offset)
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
    (ds, ch, btch, bodyCHMap, offset);
}

void BulletBroadphase_impl::updateShape(const BodyCHRecord &record)
{
  SP::SiconosConvexHull ch(record.shape);
  SP::BTCHSHAPE btch(record.btshape);

  // Update shape parameters

  // TODO
  //btbox->setLocalScaling(btVector3(sx, sy, sz));
  btch->setMargin((ch->insideMargin() + ch->outsideMargin()) * _options.worldScale);

  updateShapePosition(record);
}

void BulletBroadphase_impl::createCollisionObjectsForBodyContactor(
  const SP::BodyDS ds, SP::SiconosContactor contactor)
{
  SP::SiconosContactor con(contactor);
  if (ds) con = ds->contactor();
  if (!con) return;

  std::vector<std::pair<SP::SiconosPlane, SP::SiconosVector> >::const_iterator itp;
  for (itp=con->planes().begin(); itp!=con->planes().end(); itp++)
    createCollisionObject(ds, itp->first, itp->second);

  std::vector<std::pair<SP::SiconosSphere, SP::SiconosVector> >::const_iterator its;
  for (its=con->spheres().begin(); its!=con->spheres().end(); its++)
    createCollisionObject(ds, its->first, its->second);

  std::vector<std::pair<SP::SiconosBox, SP::SiconosVector> >::const_iterator itb;
  for (itb=con->boxes().begin(); itb!=con->boxes().end(); itb++)
    createCollisionObject(ds, itb->first, itb->second);

  std::vector<std::pair<SP::SiconosConvexHull, SP::SiconosVector> >::const_iterator itc;
  for (itc=con->convexhulls().begin(); itc!=con->convexhulls().end(); itc++)
    createCollisionObject(ds, itc->first, itc->second);
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
BulletBroadphase* BulletBroadphase::gBulletBroadphase = NULL;
bool BulletBroadphase::bulletContactClear(void* userPersistentData)
{
  /* note: stored pointer to shared_ptr! */
  SP::Interaction *p_inter = (SP::Interaction*)userPersistentData;
  assert(p_inter!=NULL && "Contact point's stored (SP::Interaction*) is null!");
  DEBUG_PRINTF("unlinking interaction %p\n", &**p_inter);
  gBulletBroadphase->unlink(*p_inter);
  delete p_inter;
  return false;
}

void BulletBroadphase::updateInteractions(SP::Simulation simulation)
{
  resetStatistics();

  // Update static contactors, initial bodies have already been
  // added via the visitor.
  impl->insertQueuedContactors(*this);

  // 0. set up bullet callbacks
  gBulletBroadphase = this;
  gContactDestroyedCallback = this->bulletContactClear;

  // TODO: This must be either configured dynamically or made available to the
  // user.
  gContactBreakingThreshold = _options.breakingThreshold;

  // 1. perform bullet collision detection
  impl->_collisionWorld->performDiscreteCollisionDetection();
  gBulletBroadphase = 0;

  if (!simulation)
    return;

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

      int g1 = pairA->sshape->group();
      int g2 = pairB->sshape->group();
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
        link(inter, pairA->ds, pairB->ds);
      }
    }
  }
}

struct CollisionUpdater : public SiconosVisitor
{
  using SiconosVisitor::visit;
  BulletBroadphase_impl &impl;
  BulletBroadphase &broad;

  CollisionUpdater(BulletBroadphase &_broad, BulletBroadphase_impl &_impl)
    : broad(_broad), impl(_impl) {}

  void visit(SP::BodyDS bds)
  {
    if (bds->contactor()) {
      BodyBoxMap::iterator it = impl.bodyBoxMap.find(&*bds);
      if (impl.bodyBoxMap.find(&*bds) == impl.bodyBoxMap.end()
          && impl.bodyCHMap.find(&*bds) == impl.bodyCHMap.end()
          && impl.bodySphereMap.find(&*bds) == impl.bodySphereMap.end()
          && impl.bodyPlaneMap.find(&*bds) == impl.bodyPlaneMap.end())
      {
        impl.createCollisionObjectsForBodyContactor(bds);
      }
      impl.updateAllShapesForDS(*bds);
    }
  }
};

SP::SiconosVisitor BulletBroadphase::getDynamicalSystemsVisitor(SP::Simulation simulation)
{
  return SP::SiconosVisitor(new CollisionUpdater(*this, *impl));
}

void BulletBroadphase::insertNonSmoothLaw(SP::NonSmoothLaw nslaw,
                                          int group1, int group2)
{
  impl->nslaws[std::pair<int,int>(group1,group2)] = nslaw;
}

SP::NonSmoothLaw BulletBroadphase::nonSmoothLaw(int group1, int group2)
{
  try {
    return impl->nslaws.at(std::pair<int,int>(group1,group2));
  } catch (const std::out_of_range &) {
    return SP::NonSmoothLaw();
  }
}
