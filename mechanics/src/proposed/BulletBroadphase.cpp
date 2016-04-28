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

#include <Model.hpp>
#include <Relation.hpp>
#include <Simulation.hpp>
#include <NonSmoothDynamicalSystem.hpp>
#include <SimulationTypeDef.hpp>
#include <NonSmoothLaw.hpp>
#include <OneStepIntegrator.hpp>
#include <NewtonImpactFrictionNSL.hpp>
#include <FrictionContact.hpp>

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

//#define DEBUG_MESSAGES 1
#include <debug.h>

// We can replace the primitives by alternative implementations.  To date,
// everything works (under test conditions, very tentative) except
// btStaticPlaneShape, so we replace it with a large box.

#define USE_CONVEXHULL_FOR_BOX 1
#define USE_CONVEXHULL_FOR_SPHERE 1
// #define USE_BOX_FOR_PLANE 1
#define USE_CONVEXHULL_FOR_PLANE 1

class BulletBroadphase_impl : public SiconosShapeHandler
{
protected:
  SP::btCollisionWorld _collisionWorld;
  SP::btDefaultCollisionConfiguration _collisionConfiguration;
  SP::btCollisionDispatcher _dispatcher;
  SP::btBroadphaseInterface _broadphase;

  SP::SiconosContactor currentContactor;
  SP::BodyDS currentBodyDS;

  std::vector<SP::SiconosPlane> dirtyPlanes;
  std::vector<SP::SiconosSphere> dirtySpheres;
  std::vector<SP::SiconosBox> dirtyBoxes;

  std::map<SP::SiconosShape, SP::btCollisionObject> objectMap;
#ifdef USE_BOX_FOR_PLANE
  std::map<SP::SiconosPlane, SP::btBoxShape> planeMap;
#else
#ifdef USE_CONVEXHULL_FOR_PLANE
  std::map<SP::SiconosPlane, SP::btConvexHullShape> planeMap;
#else
  std::map<SP::SiconosPlane, SP::btStaticPlaneShape> planeMap;
#endif
#endif

#ifdef USE_CONVEXHULL_FOR_SPHERE
  std::map<SP::SiconosSphere, SP::btConvexHullShape> sphereMap;
#else
  std::map<SP::SiconosSphere, SP::btSphereShape> sphereMap;
#endif

#ifdef USE_CONVEXHULL_FOR_BOX
  std::map<SP::SiconosBox, SP::btConvexHullShape> boxMap;
#else
  std::map<SP::SiconosBox, SP::btBoxShape> boxMap;
#endif

  std::map<Interaction*, bool> orphanedInteractions;

  SP::NonSmoothLaw nslaw;

public:
  BulletBroadphase_impl() {}
  ~BulletBroadphase_impl() {}

  virtual void onChanged(SP::SiconosPlane plane);
  virtual void onChanged(SP::SiconosSphere sphere);
  virtual void onChanged(SP::SiconosBox box);

  friend class BulletBroadphase;
};

struct BodyShapePair
{
  BodyShapePair(SP::BodyDS d, SP::SiconosShape sh) : ds(d), shape(sh) {}
  SP::BodyDS ds;
  SP::SiconosShape shape;
};

void BulletBroadphase_impl::onChanged(SP::SiconosPlane plane)
{
  dirtyPlanes.push_back(plane);
}

void BulletBroadphase_impl::onChanged(SP::SiconosSphere sphere)
{
  dirtySpheres.push_back(sphere);
}

void BulletBroadphase_impl::onChanged(SP::SiconosBox box)
{
  dirtyBoxes.push_back(box);
}

BulletBroadphase::BulletBroadphase(const BulletOptions &_options)
  : options(_options)
{
  impl.reset(new BulletBroadphase_impl());
  impl->_collisionConfiguration.reset(
    new btDefaultCollisionConfiguration());

  impl->_collisionConfiguration->setConvexConvexMultipointIterations();
  impl->_collisionConfiguration->setPlaneConvexMultipointIterations();

  impl->_dispatcher.reset(
    new btCollisionDispatcher(&*impl->_collisionConfiguration));

  if (options.useAxisSweep3)
    impl->_broadphase.reset(new btAxisSweep3(btVector3(), btVector3()));
  else
    impl->_broadphase.reset(new btDbvtBroadphase());

  impl->_collisionWorld.reset(
    new btCollisionWorld(&*impl->_dispatcher, &*impl->_broadphase,
                         &*impl->_collisionConfiguration));

  btGImpactCollisionAlgorithm::registerAlgorithm(&*impl->_dispatcher);
  impl->_collisionWorld->getDispatchInfo().m_useContinuous = false;
  impl->_collisionWorld->getDispatchInfo().m_convexConservativeDistanceThreshold = 0.1f;

  impl->nslaw.reset(new NewtonImpactFrictionNSL(0.8, 0., 0.0, 3));
}

BulletBroadphase::~BulletBroadphase()
{
  // unlink() will be called on all remaining
  // contact points when world is destroyed
  gBulletBroadphase = this;

  // must be the first de-allocated, otherwise segfault
  impl->_collisionWorld.reset();

  // free memory used to point back to DS
  std::map<SP::SiconosShape, SP::btCollisionObject>::iterator it;
  for (it = impl->objectMap.begin(); it != impl->objectMap.end(); ++it)
  {
    BodyShapePair *pair =
      static_cast<BodyShapePair*>(it->second->getUserPointer());
    if (pair) {
      delete pair;
      it->second->setUserPointer(0);
    }
  }
}

void BulletBroadphase::buildGraph(SP::Model model)
{
  // required later in performBroadphase().
  _model = model;

  DynamicalSystemsGraph& dsg =
    *(model->nonSmoothDynamicalSystem()->dynamicalSystems());
  DynamicalSystemsGraph::VIterator dsi, dsiend;
  std11::tie(dsi, dsiend) = dsg.vertices();

  for (; dsi != dsiend; ++dsi)
  {
    SP::DynamicalSystem ds(dsg.bundle(*dsi));
    ds->accept(*this);
  }
}

void BulletBroadphase::buildGraph(std::vector<SP::BodyDS> bodies)
{
  std::vector<SP::BodyDS>::iterator it;
  for (it=bodies.begin(); it!=bodies.end(); ++it)
    (*it)->accept(*this);
}

void BulletBroadphase::buildGraph(SP::SiconosContactor contactor)
{
  impl->currentBodyDS.reset();
  impl->currentContactor = contactor;

  std::vector<SP::SiconosShape>::const_iterator it;
  for (it=contactor->shapes().begin();
       it!=contactor->shapes().end();
       it++)
  {
    (*it)->acceptSP(shared_from_this());
  }
}

template<typename ST, typename BT>
void BulletBroadphase::visit_helper(ST& shape, BT& btshape,
                                    std::map<ST,BT>& shapemap)
{
  // create corresponding Bullet object and shape
  SP::btCollisionObject btobject(new btCollisionObject());

  // track association (and to keep a reference)
  impl->objectMap[shape] = SP::btCollisionObject(btobject);

  // associate the shape with the object
  btobject->setCollisionShape(&*btshape);

  // allow Bullet to report colliding DSs
  btobject->setUserPointer(
    static_cast<void*>(new BodyShapePair(impl->currentBodyDS, shape)));

  // track association (and to keep a reference)
  shapemap[shape] = btshape;

  // put it in the world
  impl->_collisionWorld->addCollisionObject(&*btobject);

  // install handler for updates
  shape->setHandler(impl);

  // initial update of the shape properties
  update(shape);
}

void BulletBroadphase::visit(SP::SiconosSphere sphere)
{
  DEBUG_PRINTF("contactor: %p, ", impl->currentContactor);
  DEBUG_PRINTF("sphere: %p(%ld)\n",
         &*sphere,sphere.use_count());

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
  visit_helper(sphere, btsphere, impl->sphereMap);
}

void BulletBroadphase::update(SP::SiconosSphere sphere)
{
  const SP::SiconosVector pos = sphere->position();
  DEBUG_PRINTF("updating sphere: %p(%ld) - %f, %f, %f (r=%0.2f)\n",
         &*sphere,sphere.use_count(),
         (*pos)(0), (*pos)(1), (*pos)(2), sphere->radius());

  double r = (sphere->radius() + sphere->outsideMargin()) * options.worldScale;

  // Update shape parameters
#ifdef USE_CONVEXHULL_FOR_SPHERE
  SP::btConvexHullShape btsphere(impl->sphereMap[sphere]);
  assert(btsphere
         && "BulletBroadphase::update(), sphere not found in sphereMap.");
  btsphere->setMargin(r);
#else
  SP::btSphereShape btsphere(impl->sphereMap[sphere]);
  assert(btsphere
         && "BulletBroadphase::update(), sphere not found in sphereMap.");

  btsphere->setLocalScaling(btVector3(r, r, r));

  // btSphereShape has an internal margin
  btsphere->setMargin(sphere->insideMargin() * options.worldScale);
#endif

  // Update object parameters
  SP::btCollisionObject btobject(impl->objectMap[sphere]);
  assert(btobject
         && "BulletBroadphase::update(), sphere not found in objectMap.");

  SiconosVector &q = *sphere->position();
  btTransform tr(btQuaternion(q(4), q(5), q(6), q(3)),
                 btVector3(q(0)*options.worldScale,
                           q(1)*options.worldScale,
                           q(2)*options.worldScale));
  btobject->setWorldTransform(tr);
}

void BulletBroadphase::visit(SP::SiconosPlane plane)
{
  DEBUG_PRINTF("contactor: %p, ", impl->currentContactor);
  DEBUG_PRINTF("plane: %p(%ld)\n",
         &*plane,plane.use_count());

  // create the initial plane with default parameters
#ifdef USE_BOX_FOR_PLANE
  SP::btBoxShape btplane(
    new btBoxShape(btVector3(1000*options.worldScale,
                             1000*options.worldScale,
                             1000*options.worldScale)));
#else
#ifdef USE_CONVEXHULL_FOR_PLANE
  btScalar h = 1000 * options.worldScale;
  const btScalar pts[] = {
    h, h, 0,
    h, -h, 0,
    -h, -h, 0,
    -h, h, 0,
  };
  SP::btConvexHullShape btplane(
    new btConvexHullShape(pts, 4, sizeof(pts[0])*3));

  // We ignore the desired internal margin for plane and just use a large one.
  plane->setInsideMargin(1000 * options.worldScale);

  // External margin
  btplane->setMargin((plane->insideMargin() + plane->outsideMargin())
                     * options.worldScale);
#else
  SP::btStaticPlaneShape btplane(
    new btStaticPlaneShape(btVector3(0, 0, 1), 0.0));
  btplane->setMargin(plane->outsideMargin() * options.worldScale);
#endif
#endif

  // initialization
  visit_helper(plane, btplane, impl->planeMap);
}

void BulletBroadphase::update(SP::SiconosPlane plane)
{
  DEBUG_PRINTF("updating plane: %p(%ld)\n",
         &*plane,plane.use_count());

  // Update object parameters
  SP::btCollisionObject btobject(impl->objectMap[plane]);
  assert(btobject
         && "BulletBroadphase::update(), plane not found in objectMap.");

  SiconosVector &q = *plane->position();
#ifdef USE_BOX_FOR_PLANE
  btTransform tr(btQuaternion(q(4), q(5), q(6), q(3)),
    btVector3(q(0)*options.worldScale,
              q(1)*options.worldScale,
              (q(2) + plane->outsideMargin() - 1000)*options.worldScale));
#else
#ifdef USE_CONVEXHULL_FOR_PLANE
  btTransform tr(btQuaternion(q(4), q(5), q(6), q(3)),
    btVector3(q(0)*options.worldScale,
              q(1)*options.worldScale,
              (q(2) - plane->insideMargin())*options.worldScale));
#else
  btTransform tr(btQuaternion(q(4), q(5), q(6), q(3)),
                 btVector3(q(0)*options.worldScale,
                           q(1)*options.worldScale,
              (q(2) - plane->insideMargin())*options.worldScale));
#endif
#endif
  btobject->setWorldTransform(tr);
}

void BulletBroadphase::visit(SP::SiconosBox box)
{
  DEBUG_PRINTF("contactor: %p, ", impl->currentContactor);
  DEBUG_PRINTF("box: %p(%ld)\n",
         &*box,box.use_count());

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
  btbox->setMargin(box->insideMargin() * options.worldScale);
#endif

  // initialization
  visit_helper(box, btbox, impl->boxMap);
}

void BulletBroadphase::update(SP::SiconosBox box)
{
  const SP::SiconosVector pos = box->position();
  DEBUG_PRINTF("updating box: %p(%ld) - %f, %f, %f\n",
         &*box,box.use_count(),
         (*pos)(0), (*pos)(1), (*pos)(2));

  // Update shape parameters
#ifdef USE_CONVEXHULL_FOR_BOX
  SP::btConvexHullShape btbox(impl->boxMap[box]);
#else
  SP::btBoxShape btbox(impl->boxMap[box]);
#endif
  assert(btbox
         && "BulletBroadphase::update(), box not found in boxMap.");

#ifdef USE_CONVEXHULL_FOR_BOX
  double m = -box->insideMargin();
#else
  double m = box->outsideMargin();
#endif

  double sx = ((*box->dimensions())(0) + m*2) * options.worldScale;
  double sy = ((*box->dimensions())(1) + m*2) * options.worldScale;
  double sz = ((*box->dimensions())(2) + m*2) * options.worldScale;

  btbox->setLocalScaling(btVector3(sx, sy, sz));
  btbox->setMargin((box->insideMargin() + box->outsideMargin()) * options.worldScale);

  // Update object parameters
  SP::btCollisionObject btobject(impl->objectMap[box]);
  assert(btobject
         && "BulletBroadphase::update(), box not found in objectMap.");

  SiconosVector &q = *box->position();
  btTransform tr(btQuaternion(q(4), q(5), q(6), q(3)),
                 btVector3(q(0)*options.worldScale,
                           q(1)*options.worldScale,
                           q(2)*options.worldScale));

  btobject->setWorldTransform(tr);
}

void BulletBroadphase::visit(const BodyDS &bds)
{
  SP::SiconosContactor contactor = bds.contactor();

  // (Stripping constness!  Only because SP::BodyDS visitor does not work.)
  SP::BodyDS spbds( ((BodyDS*)(&bds))->shared_from_this() );

  impl->currentBodyDS = spbds;
  impl->currentContactor = contactor;

  std::vector<SP::SiconosShape>::const_iterator it;
  for (it=contactor->shapes().begin();
       it!=contactor->shapes().end();
       it++)
  {
    (*it)->acceptSP(shared_from_this());
    (*it)->setPosition(bds.q());
  }
}

void BulletBroadphase::updateGraph()
{
  if (!impl->dirtyPlanes.empty())
  {
    std::vector<SP::SiconosPlane>::iterator it;
    for (it=impl->dirtyPlanes.begin();
         it!=impl->dirtyPlanes.end(); it++)
    {
      update(*it);
    }
    impl->dirtyPlanes.clear();
  }

  if (!impl->dirtySpheres.empty())
  {
    std::vector<SP::SiconosSphere>::iterator it;
    for (it=impl->dirtySpheres.begin();
         it!=impl->dirtySpheres.end(); it++)
    {
      update(*it);
    }
    impl->dirtySpheres.clear();
  }

  if (!impl->dirtyBoxes.empty())
  {
    std::vector<SP::SiconosBox>::iterator it;
    for (it=impl->dirtyBoxes.begin();
         it!=impl->dirtyBoxes.end(); it++)
    {
      update(*it);
    }
    impl->dirtyBoxes.clear();
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
BulletBroadphase* BulletBroadphase::gBulletBroadphase = NULL;
bool BulletBroadphase::bulletContactClear(void* userPersistentData)
{
  /* note: stored pointer to shared_ptr! */
  SP::Interaction *p_inter = (SP::Interaction*)userPersistentData;
  assert(p_inter!=NULL && "Contact point's stored (SP::Interaction*) is null!");
  DEBUG_PRINTF("unlinking interaction %p\n", &**p_inter);
  gBulletBroadphase->unlink(*p_inter);
  delete p_inter;
}

void BulletBroadphase::performBroadphase()
{
  // 0. set up bullet callbacks
  gBulletBroadphase = this;
  gContactDestroyedCallback = this->bulletContactClear;

  // TODO: This must be either configured dynamically or made available to the
  // user.
  gContactBreakingThreshold = options.breakingThreshold;

  // 1. perform bullet collision detection
  impl->orphanedInteractions.clear();
  impl->_collisionWorld->performDiscreteCollisionDetection();
  gBulletBroadphase = 0;

  if (!model())
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
    const BodyShapePair *pairA, *pairB;
    pairA = static_cast<BodyShapePair*>(it->objectA->getUserPointer());
    pairB = static_cast<BodyShapePair*>(it->objectB->getUserPointer());
    assert(pairA && pairB && "btCollisionObject had a null user pointer!");

    // The first pair will always be the non-static object
    if (pairB->ds && !pairA->ds) {
      const BodyShapePair *tmp = pairA;
      pairA = pairB;
      pairB = tmp;
    }

    // If both bodies are static, no interaction is created.
    if (!(pairA->ds || pairB->ds))
      continue;

    if (it->point->m_userPersistentData)
    {
      // do what's needed for an interaction already present
      SP::Interaction *p_inter =
        (SP::Interaction*)it->point->m_userPersistentData;
      // (note: nothing for now!)
      stats.existing_interactions_processed ++;
    }
    else
    {
      /* new interaction */
      SP::Interaction inter;
      if (impl->nslaw->size() == 3)
      {
        // Remove the added outside margin as a correction factor in Relation
        double combined_margin =
          pairA->shape->outsideMargin() + pairB->shape->outsideMargin();

        SP::BulletR rel(new BulletR(createSPtrbtManifoldPoint(*it->point),
                                    pairA->shape->outsideMargin(),
                                    pairB->shape->outsideMargin(),
                                    1.0 / options.worldScale));

        // We wish to be sure that no Interactions are created without
        // sufficient warning before contact.  TODO: Replace with exception or
        // flag.
        if ((it->point->getDistance() + combined_margin) < 0.0) {
          fprintf(stderr, "Interactions must be created with positive distance (%f).\n",
                  (it->point->getDistance() + combined_margin)/options.worldScale);
          stats.interaction_warnings ++;
        }

        inter.reset(new Interaction(3, impl->nslaw, rel, 0 /*4 * i + z*/));
        stats.new_interactions_created ++;
      }
      else
      {
        if (impl->nslaw->size() == 1)
        {
          SP::BulletFrom1DLocalFrameR rel(
            new BulletFrom1DLocalFrameR(createSPtrbtManifoldPoint(*it->point)));
          inter.reset(new Interaction(1, impl->nslaw, rel, 0 /*4 * i + z*/));
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
    // if (n_points++ < 4) {
    //   printf(", %f", it->point->getDistance());
    // }
  }

  // while (n_points++ < 4)
  //   printf(", 0");
  // printf(", %d\n", late_interaction);

  /* Update non smooth problem */
  model()->simulation()->initOSNS();
}
