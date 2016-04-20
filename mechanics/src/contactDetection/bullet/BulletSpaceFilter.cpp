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

#include "BulletSpaceFilter.hpp"
#include "BulletSpaceFilter_impl.hpp"
#include "SpaceFilter_impl.hpp"

#include <Model.hpp>
#include <Simulation.hpp>
#include <NonSmoothDynamicalSystem.hpp>
#include <SimulationTypeDef.hpp>
#include <NonSmoothLaw.hpp>
#include <OneStepIntegrator.hpp>

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

extern btScalar gContactBreakingThreshold;
typedef bool (*ContactProcessedCallback)(btManifoldPoint& cp, void* body0, void* body1);
typedef bool (*ContactDestroyedCallback)(void* userPersistentData);
extern ContactDestroyedCallback gContactDestroyedCallback;
extern ContactProcessedCallback gContactProcessedCallback;


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

        nslaw = (*_nslaws)(gid1, gid2);

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


          if (itc == contactPoints.end() || !cpoint->m_userPersistentData)
          {
            /* new interaction */

            SP::Interaction inter;
            if (nslaw->size() == 3)
            {
              SP::BulletR rel(new BulletR(cpoint));
              inter.reset(new Interaction(3, nslaw, rel, 4 * i + z));
            }
            else
            {
              if (nslaw->size() == 1)
              {
              SP::BulletFrom1DLocalFrameR rel(new BulletFrom1DLocalFrameR(cpoint));
              inter.reset(new Interaction(1, nslaw, rel, 4 * i + z));
              }
            }

            if (obB->getUserPointer())
            {
              SP::BulletDS dsb(static_cast<BulletDS*>(obB->getUserPointer())->shared_ptr());

              if (dsa != dsb)
              {
                DEBUG_PRINTF("LINK obA:%p obB:%p inter:%p\n", obA, obB, &*inter);
                assert(inter);

                cpoint->m_userPersistentData = &*inter;
                link(inter, dsa, dsb);
              }
              /* else collision shapes belong to the same object do nothing */
            }
            else
            {
              DEBUG_PRINTF("LINK obA:%p inter :%p\n", obA, &*inter);
              assert(inter);

              cpoint->m_userPersistentData = &*inter;
              link(inter, dsa);
            }
          }

          if (cpoint->m_userPersistentData)
          {
            activeInteractions[static_cast<Interaction *>(cpoint->m_userPersistentData)] = true;
            DEBUG_PRINTF("Interaction %p = true\n", static_cast<Interaction *>(cpoint->m_userPersistentData));
            DEBUG_PRINTF("cpoint %p  = true\n", &*cpoint);
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
  //osi->insertDynamicalSystem(ds);
  this->model()->nonSmoothDynamicalSystem()->insertDynamicalSystem(ds);
  simulation->addInOSIMap(ds, osi);

  /* Initialize the DS at the current time */
  ds->initialize(simulation->nextTime(), osi->getSizeMem());

  /* Partially re-initialize the simulation. */
  simulation->initialize(this->model(), false);

  /* Re-create the world from scratch */
  _dispatcher.reset(new btCollisionDispatcher(&*_collisionConfiguration));
  _collisionWorld.reset(new btCollisionWorld(&*_dispatcher, &*_broadphase,
                                             &*_collisionConfiguration));
  _dynamicCollisionsObjectsInserted = false;
  _staticCollisionsObjectsInserted = false;
}
