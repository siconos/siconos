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

#include "BulletSpaceFilter.hpp"
#include "BulletSpaceFilter_impl.hpp"
#include "SpaceFilter_impl.hpp"

#include <Model.hpp>
#include <Simulation.hpp>
#include <NonSmoothDynamicalSystem.hpp>
#include <SimulationTypeDef.hpp>
#include <NonSmoothLaw.hpp>
#include <OneStepIntegrator.hpp>
#include <Topology.hpp>

#include <BulletCollision/Gimpact/btGImpactCollisionAlgorithm.h>
#include <BulletCollision/CollisionDispatch/btDefaultCollisionConfiguration.h>
#include <BulletCollision/BroadphaseCollision/btDbvtBroadphase.h>

#include <Question.hpp>
#include "BulletDS.hpp"
#include "BulletDS_impl.hpp"

#include <RuntimeException.hpp>
#include <boost/format.hpp>
#include <boost/typeof/typeof.hpp>

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES 1
#include <debug.h>


#include "MoreauJeanOSI.hpp"
#include "MoreauJeanGOSI.hpp"
#include "EulerMoreauOSI.hpp"
#include "SchatzmanPaoliOSI.hpp"

#define VISITOR_CLASSES() \
  REGISTER(MoreauJeanOSI) \
  REGISTER(MoreauJeanGOSI) \
  REGISTER(EulerMoreauOSI) \
  REGISTER(SchatzmanPaoliOSI)

#include <VisitorMaker.hpp>

using namespace Experimental;

extern btScalar gContactBreakingThreshold;
typedef bool (*ContactProcessedCallback)(btManifoldPoint& cp, void* body0, void* body1);
typedef bool (*ContactDestroyedCallback)(void* userPersistentData);
extern ContactDestroyedCallback gContactDestroyedCallback;
extern ContactProcessedCallback gContactProcessedCallback;

int btScalarSize()
{
  return sizeof(btScalar);
}

std::map<Interaction*, bool> gOrphanedInteractions;

bool contactClear(void* userPersistentData);
bool contactClear(void* userPersistentData)
{
  if (true)
  {
    DEBUG_PRINTF("contactClear : Interaction : %p\n",   static_cast<Interaction *>(userPersistentData));
    gOrphanedInteractions[static_cast<Interaction *>(userPersistentData)] = true;
  }

  return true;
}

bool contactProcess(btManifoldPoint& cp, void *body0, void *body1);
bool contactProcess(btManifoldPoint& cp, void *body0, void *body1)
{
  if (true)
  {
    DEBUG_PRINTF("gContactProcessedCallback : process contactPoint : %p, distance %g\n", &cp, cp.getDistance());
  }
  return true;
}

#include "BulletDS.hpp"

struct ForPosition : public Question<SP::SiconosVector>
{
  ANSWER(BulletDS, q());
};

// /* initializeIterationMatrixW is not a member of OneStepIntegrator (should it be ?),
//    so we visit some integrators which provide this initialization.
// */

// /* first, a generic visitor is defined. */
// struct CallInitW : public SiconosVisitor
// {
//   double time;
//   SP::DynamicalSystem ds;
//   DynamicalSystemsGraph::VDescriptor dsv;

//   template<typename T>
//   void operator()(const T& osi)
//   {
//     const_cast<T*>(&osi)->initializeIterationMatrixW(this->time, this->ds, this->dsv);
//   }
// };

// /* the visit is made on classes which provide the function initializeIterationMatrixW */
// typedef Visitor < Classes < MoreauJeanOSI,
//                             MoreauJeanGOSI,
//                             EulerMoreauOSI,
//                             SchatzmanPaoliOSI >,
//                   CallInitW >::Make InitW;


/* initializeIterationMatrixW is not a member of OneStepIntegrator (should it be ?),
   so we visit some integrators which provide this initialization.
*/

/* first, a generic visitor is defined. */
struct CallInitDS : public SiconosVisitor
{
  double time;
  SP::DynamicalSystem ds;
  SP::Model m ; 
  template<typename T>
  void operator()(const T& osi)
  {
    const_cast<T*>(&osi)->initializeDynamicalSystem(*(this->m), this->time, this->ds);
  }
};

/* the visit is made on classes which provide the function initializeIterationMatrixW */
// typedef Visitor < Classes < MoreauJeanOSI >,
//                   CallInitDS >::Make InitDynamicalSystem;

typedef Visitor < Classes < MoreauJeanOSI,
                            MoreauJeanGOSI,
                            EulerMoreauOSI,
                            SchatzmanPaoliOSI >,
                  CallInitDS >::Make InitDynamicalSystem;

BulletSpaceFilter::BulletSpaceFilter(SP::Model model) :
  SpaceFilter(),
  _dynamicCollisionsObjectsInserted(false),
  _staticCollisionsObjectsInserted(false),
  _closeContactsThreshold(0.)
{

  _model = model;
  _nslaws.reset(new NSLawMatrix());
  _staticObjects.reset(new StaticObjects());

  _collisionConfiguration.reset(new btDefaultCollisionConfiguration());

  /* not done by default anymore */
  /* one have to get collision configuration and call explicitely
   * these methods */

/*  _collisionConfiguration->setConvexConvexMultipointIterations();
    _collisionConfiguration->setPlaneConvexMultipointIterations();*/

  _dispatcher.reset(new btCollisionDispatcher(&*_collisionConfiguration));

  _broadphase.reset(new btDbvtBroadphase());

  _collisionWorld.reset(new btCollisionWorld(&*_dispatcher, &*_broadphase,
                                             &*_collisionConfiguration));

  btGImpactCollisionAlgorithm::registerAlgorithm(&*_dispatcher);

  _collisionWorld->getDispatchInfo().m_useContinuous = false;

  //gContactBreakingThreshold = 10.;

  gContactProcessedCallback = &contactProcess;
  gContactDestroyedCallback = &contactClear;

}

BulletSpaceFilter::~BulletSpaceFilter()
{
    // btCollisionObjects contain pointers to _broadphase (so-called
    // "broadphase handles") that are referenced during the
    // _collisionWorld destructor, therefore if _broadphase is
    // destroyed too early we get a segfault.  Avoid this by
    // explicitly destroying _collisionWorld before automatic
    // shared_ptr destruction.
    _collisionWorld.reset();
}

void BulletSpaceFilter::setCollisionConfiguration(
  SP::btDefaultCollisionConfiguration collisionConfig)
{
  _collisionConfiguration = collisionConfig;
  _dispatcher.reset(new btCollisionDispatcher(&*_collisionConfiguration));
  _collisionWorld.reset(new btCollisionWorld(&*_dispatcher, &*_broadphase, &*_collisionConfiguration));
  _dynamicCollisionsObjectsInserted = false;
  _staticCollisionsObjectsInserted = false;
}

void BulletSpaceFilter::buildInteractions(double time)
{
  DEBUG_PRINT("-----start build interaction\n");

  DEBUG_PRINT("----- insert dynamic collision objects if needed\n");

  if (! _dynamicCollisionsObjectsInserted ||  ! _staticCollisionsObjectsInserted)
  {
    _collisionWorld.reset();

    _dispatcher.reset(new btCollisionDispatcher(&*_collisionConfiguration));
    _collisionWorld.reset(new btCollisionWorld(&*_dispatcher, &*_broadphase,
                                               &*_collisionConfiguration));
    // btGImpactCollisionAlgorithm::registerAlgorithm(&*_dispatcher);
    // _collisionWorld->getDispatchInfo().m_useContinuous = false;
    _dynamicCollisionsObjectsInserted = false;
    _staticCollisionsObjectsInserted = false;
  }

  if (! _dynamicCollisionsObjectsInserted)
  {
    DynamicalSystemsGraph& dsg = *(_model->nonSmoothDynamicalSystem()->dynamicalSystems());
    DynamicalSystemsGraph::VIterator dsi, dsiend;
    std11::tie(dsi, dsiend) = dsg.vertices();
    for (; dsi != dsiend; ++dsi)
    {
      CollisionObjects& collisionObjects = *ask<ForCollisionObjects>(*dsg.bundle(*dsi));

      for (CollisionObjects::iterator ico = collisionObjects.begin();
           ico != collisionObjects.end(); ++ico)
      {
        _collisionWorld->addCollisionObject(const_cast<btCollisionObject*>((*ico).first));
        DEBUG_PRINTF("add dynamic collision object %p\n", const_cast<btCollisionObject*>((*ico).first));
      }
    }

    _dynamicCollisionsObjectsInserted = true;
  }

  if (! _staticCollisionsObjectsInserted)
  {
    for(StaticObjects::iterator
          ic = _staticObjects->begin(); ic != _staticObjects->end(); ++ic)
    {
      _collisionWorld->addCollisionObject(const_cast<btCollisionObject*>((*ic).first));
      DEBUG_PRINTF("add static collision object %p\n", const_cast<btCollisionObject*>((*ic).first));
    }

    _staticCollisionsObjectsInserted = true;
  }

  //  1. perform bullet collision detection
  DEBUG_PRINT("-----  1. perform bullet collision detection\n");
  gOrphanedInteractions.clear();
  _collisionWorld->performDiscreteCollisionDetection();



  // 2. collect old contact points from Siconos graph
  DEBUG_PRINT("-----  2. collect old contact points from Siconos graph\n");
  std::map<btManifoldPoint*, bool> contactPoints;

  std::map<Interaction*, bool> activeInteractions;

  SP::InteractionsGraph indexSet0 = model()->nonSmoothDynamicalSystem()->topology()->indexSet(0);
  InteractionsGraph::VIterator ui0, ui0end, v0next;
  std11::tie(ui0, ui0end) = indexSet0->vertices();
  for (v0next = ui0 ;
       ui0 != ui0end; ui0 = v0next)
  {
    ++v0next;  // trick to iterate on a dynamic bgl graph
    SP::Interaction inter0 = indexSet0->bundle(*ui0);

    if (gOrphanedInteractions.find(&*inter0) != gOrphanedInteractions.end())
    {
      model()->nonSmoothDynamicalSystem()->removeInteraction(inter0);
    }

    else
    {
      SP::btManifoldPoint cp = ask<ForContactPoint>(*(inter0->relation()));
      if(cp)
      {
        DEBUG_PRINTF("Contact point in interaction : %p\n", &*cp);

        contactPoints[&*cp] = false;
      }
    }
  }
  unsigned int numManifolds =
    _collisionWorld->getDispatcher()->getNumManifolds();

  DEBUG_PRINTF("number of manifolds : %d\n", numManifolds);

  for (unsigned int i = 0; i < numManifolds; ++i)
  {
    btPersistentManifold* contactManifold =
      _collisionWorld->getDispatcher()->getManifoldByIndexInternal(i);

    DEBUG_PRINTF("processing manifold %d : %p\n", i, contactManifold);

    const btCollisionObject* obA = contactManifold->getBody0();
    const btCollisionObject* obB = contactManifold->getBody1();

    //    contactManifold->refreshContactPoints(obA->getWorldTransform(),obB->getWorldTransform());

    unsigned int numContacts = contactManifold->getNumContacts();

    if ( (obA->getUserPointer() && obA->getUserPointer() != obB->getUserPointer()) ||
         (obB->getUserPointer() && obA->getUserPointer() != obB->getUserPointer()) )
    {

      for (unsigned int z = 0; z < numContacts; ++z)
      {

        SP::btManifoldPoint cpoint(createSPtrbtManifoldPoint(contactManifold->getContactPoint(z)));
        DEBUG_PRINTF("manifold %d, contact %d, &contact %p, lifetime %d\n", i, z, &*cpoint, cpoint->getLifeTime());


        // should no be mixed with something else that use UserPointer!
        SP::BulletDS dsa;
        SP::BulletDS dsb;
        unsigned long int gid1, gid2;

        if(obA->getUserPointer())
        {
          dsa = static_cast<BulletDS*>(obA->getUserPointer())->shared_ptr();

          assert(dsa->collisionObjects()->find(contactManifold->getBody0()) !=
                 dsa->collisionObjects()->end());
          gid1 = boost::get<2>((*((*dsa->collisionObjects()).find(obA))).second);
        }
        else
        {
          gid1 = (*_staticObjects->find(obA)).second.second;
        }

        SP::NonSmoothLaw nslaw;
        if (obB->getUserPointer())
        {
          dsb = static_cast<BulletDS*>(obB->getUserPointer())->shared_ptr();

          assert(dsb->collisionObjects()->find(obB) != dsb->collisionObjects()->end());

          gid2 = boost::get<2>((*((*dsb->collisionObjects()).find(obB))).second);
        }

        else
        {
          gid2 = (*_staticObjects->find(obB)).second.second;
        }


        DEBUG_PRINTF("collision between group %ld and %ld\n", gid1, gid2);

        try {
          nslaw = (*_nslaws)(gid1, gid2);
        } catch (ublas::bad_index &e) {
          DEBUG_PRINTF("Warning: NonSmoothLaw for groups %u and %u not found!\n",
                       gid1, gid2);
        }

        if (nslaw)
        {
          std::map<btManifoldPoint*, bool>::iterator itc;
          itc = contactPoints.find(&*cpoint);

          DEBUG_EXPR(if (itc == contactPoints.end())
                     {
                       DEBUG_PRINT("contact point not found\n");
                       for(std::map<btManifoldPoint*, bool>::iterator itd=contactPoints.begin();
                           itd != contactPoints.end(); ++itd)
                       {
                         DEBUG_PRINTF("-->%p != %p\n", &*cpoint, &*(*itd).first);
                       }
                     });


          bool flip = false;
          if (itc == contactPoints.end() || !cpoint->m_userPersistentData)
          {
            /* new interaction */

            SP::Interaction inter;
            if (nslaw->size() == 3)
            {
              /* if objectB is the only DS, (A is static), then flip
               * the contact points and normal otherwise the relation
               * is to the wrong side */
              flip = !dsa && dsb;
              SP::BulletR rel(new BulletR(*cpoint,
                                          flip ? dsb->q() : dsa->q(),
                                          (flip ? (dsa?dsa->q():SP::SiconosVector())
                                                : (dsb?dsb->q():SP::SiconosVector())),
                                          flip));
              rel->setContactPoint(cpoint);
              inter.reset(new Interaction(nslaw, rel));//, 4 * i + z));
            }
            else
            {
              if (nslaw->size() == 1)
              {
              SP::BulletFrom1DLocalFrameR rel(new BulletFrom1DLocalFrameR(cpoint));
              inter.reset(new Interaction(nslaw, rel));//, 4 * i + z));
              }
            }

            if (dsa != dsb)
            {
              DEBUG_PRINTF("LINK obA:%p obB:%p inter:%p\n", obA, obB, &*inter);
            }

            if (dsa && !dsb)
            {
                cpoint->m_userPersistentData = &*inter;
                model()->simulation()->link(inter, dsa);
            }
            else if (!dsa && dsb)
            {
                cpoint->m_userPersistentData = &*inter;
                model()->simulation()->link(inter, dsb);
            }
            else if (dsa && dsb && (dsa != dsb))
            {
                cpoint->m_userPersistentData = &*inter;
                model()->simulation()->link(inter, dsa, dsb);
            }
          }

          if (cpoint->m_userPersistentData)
          {
            Interaction *inter = static_cast<Interaction *>(cpoint->m_userPersistentData);
            activeInteractions[inter] = true;
            DEBUG_PRINTF("Interaction %p = true\n", static_cast<Interaction *>(cpoint->m_userPersistentData));
            DEBUG_PRINTF("cpoint %p  = true\n", &*cpoint);
            SP::BulletR rel(std11::static_pointer_cast<BulletR>(inter->relation()));
            rel->updateContactPoints(*cpoint,
                                     flip ? dsb : dsa,
                                     flip ? (dsa ? dsa : SP::NewtonEulerDS())
                                          : (dsb ? dsb : SP::NewtonEulerDS()));
          }
          else
          {
            assert(false);
            DEBUG_PRINT("cpoint->m_userPersistentData is empty\n");
          }

          contactPoints[&*cpoint] = true;
          DEBUG_PRINTF("cpoint %p  = true\n", &*cpoint);
        }
      }
    }
  }
  // 4. remove old contact points
  std11::tie(ui0, ui0end) = indexSet0->vertices();
  for (v0next = ui0 ;
       ui0 != ui0end; ui0 = v0next)
  {
    ++v0next;  // trick to iterate on a dynamic bgl graph
    SP::Interaction inter0 = indexSet0->bundle(*ui0);
    SP::btManifoldPoint cp = ask<ForContactPoint>(*(inter0->relation()));
    if (cp)
    {
      if (!contactPoints[&*cp])
      {

        //      assert (!contactPoints[&*ask<ForContactPoint>(*(inter0->relation()))]);

        DEBUG_PRINTF("remove contact %p, lifetime %d\n",
                     &*ask<ForContactPoint>(*(inter0->relation())),
                     ask<ForContactPoint>(*(inter0->relation()))->getLifeTime());
        model()->nonSmoothDynamicalSystem()->removeInteraction(inter0);
      }
    }
  }

  DEBUG_PRINT("-----end build interaction\n");

  model()->simulation()->initOSNS();

}

void BulletSpaceFilter::addStaticObject(SP::btCollisionObject co, unsigned int id)
{
  (*_staticObjects)[&*co]= std::pair<SP::btCollisionObject, int>(co, id);
}

void BulletSpaceFilter::addDynamicObject(SP::BulletDS ds,
                                         SP::Simulation simulation,
                                         SP::OneStepIntegrator osi)
{
  if (!osi && simulation->oneStepIntegrators()
      && simulation->oneStepIntegrators()->size() > 0)
  {
      osi = *(simulation->oneStepIntegrators()->begin());
  }
  else if (!osi) {
      RuntimeException::selfThrow(
          "BulletSpaceFilter::addDynamicObject -- no OSI found");
  }

  /* Insert the new DS into the OSI, model, and simulation. */
  this->model()->nonSmoothDynamicalSystem()->insertDynamicalSystem(ds);

  /* Associate/initialize the OSI */
  simulation->prepareIntegratorForDS(osi, ds, this->model(), simulation->nextTime());

  /* Partially re-initialize the simulation. */
  simulation->initialize(this->model(), false);

  /* Re-create the world from scratch */
#if 0
  _collisionWorld.reset();

  _dispatcher.reset(new btCollisionDispatcher(&*_collisionConfiguration));
  _collisionWorld.reset(new btCollisionWorld(&*_dispatcher, &*_broadphase,
                                             &*_collisionConfiguration));
  // btGImpactCollisionAlgorithm::registerAlgorithm(&*_dispatcher);
  // _collisionWorld->getDispatchInfo().m_useContinuous = false;
  _dynamicCollisionsObjectsInserted = false;
  _staticCollisionsObjectsInserted = false;
#else
  for (CollisionObjects::iterator ico = ds->collisionObjects()->begin();
       ico != ds->collisionObjects()->end(); ++ico)
  {
    btCollisionObject *ob = const_cast<btCollisionObject*>((*ico).first);
    _collisionWorld->addCollisionObject(ob);
    _collisionWorld->updateSingleAabb(ob);
    _collisionWorld->getBroadphase()->getOverlappingPairCache()->
        cleanProxyFromPairs(ob->getBroadphaseHandle(), &*_dispatcher);
  }
#endif
}
