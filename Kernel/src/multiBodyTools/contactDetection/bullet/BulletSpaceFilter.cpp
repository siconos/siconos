/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
#include "BulletCollision/Gimpact/btGImpactCollisionAlgorithm.h"

//#define DEBUG_MESSAGES 1
#include <debug.h>

extern btScalar gContactBreakingThreshold;
typedef bool (*ContactProcessedCallback)(btManifoldPoint& cp, void* body0, void* body1);
typedef bool (*ContactDestroyedCallback)(void* userPersistentData);
extern ContactDestroyedCallback gContactDestroyedCallback;
extern ContactProcessedCallback gContactProcessedCallback;

std::vector<int> gOrphanedInteractions;

bool contactClear(void* userPersistentData)
{
  //  DEBUG_PRINTF("gContactDestroyedCallback : push Interaction number %d, distance %g\n", static_cast<Interaction *>(userPersistentData)->number(),  static_cast<Interaction *>(userPersistentData)->y(0)->getValue(0));
  gOrphanedInteractions.push_back(static_cast<Interaction *>(userPersistentData)->number());
  return true;
}

bool contactProcess(btManifoldPoint& cp, void *body0, void *body1)
{
  if (true)
  {
    DEBUG_PRINTF("gContactProcessedCallback : process contactPoint : %p, distance %g\n", &cp, cp.getDistance());
  }
  return true;
}


struct ForPosition : public Question<SP::SiconosVector>
{
  using SiconosVisitor::visit;

  ANSWER(BulletDS, q());
};

BulletSpaceFilter::BulletSpaceFilter(SP::Model model,
                                     SP::NonSmoothLaw nslaw,
                                     SP::btVector3 aabbMin,
                                     SP::btVector3 aabbMax) :
  SpaceFilter(),
  _worldAabbMin(aabbMin),
  _worldAabbMax(aabbMax)
{

  _model = model;
  _nslaw = nslaw;

  _staticObjects.reset(new std::vector<SP::btCollisionObject>());
  _staticShapes.reset(new std::vector<SP::btCollisionShape>());

  _collisionConfiguration.reset(new btDefaultCollisionConfiguration());

  _collisionConfiguration->setConvexConvexMultipointIterations();

  _dispatcher.reset(new btCollisionDispatcher(&*_collisionConfiguration));

  _broadphase.reset(new BulletBroadPhase(*_worldAabbMin, *_worldAabbMax));

  _collisionWorld.reset(new btCollisionWorld(&*_dispatcher, &*_broadphase, &*_collisionConfiguration));

  btGImpactCollisionAlgorithm::registerAlgorithm(&*_dispatcher);

  _collisionWorld->getDispatchInfo().m_useContinuous = false;

  DynamicalSystemsGraph& dsg = *(model->nonSmoothDynamicalSystem()->dynamicalSystems());
  DynamicalSystemsGraph::VIterator dsi, dsiend;
  boost::tie(dsi, dsiend) = dsg.vertices();
  for (; dsi != dsiend; ++dsi)
  {
    _collisionWorld->addCollisionObject(&*(ask<ForCollisionObject>(*(dsg.bundle(*dsi)))));
  }

  //gContactBreakingThreshold = 10.;

  gContactProcessedCallback = &contactProcess;
  gContactDestroyedCallback = &contactClear;

}


void BulletSpaceFilter::buildInteractions(double time)
{

  DEBUG_PRINT("-----start build interaction\n");

  // 1. perform bullet collision detection
  _collisionWorld->performDiscreteCollisionDetection();

  // 2. collect old contact points from Siconos graph
  std::map<btManifoldPoint*, bool> contactPoints;

  std::map<int, bool> activeInteractions;

  SP::InteractionsGraph indexSet0 = model()->nonSmoothDynamicalSystem()->topology()->indexSet(0);
  InteractionsGraph::VIterator ui0, ui0end, v0next;
  boost::tie(ui0, ui0end) = indexSet0->vertices();
  for (v0next = ui0 ;
       ui0 != ui0end; ui0 = v0next)
  {
    Interaction& inter0 = *(indexSet0->bundle(*ui0));
    ++v0next;  // trick to iterate on a dynamic bgl graph
    contactPoints[&*ask<ForContactPoint>(*(inter0.relation()))] = false;
  }

  unsigned int numManifolds =
    _collisionWorld->getDispatcher()->getNumManifolds();

  for (unsigned int i = 0; i < numManifolds; ++i)
  {
    btPersistentManifold* contactManifold =
      _collisionWorld->getDispatcher()->getManifoldByIndexInternal(i);

    btCollisionObject* obA =
      static_cast<btCollisionObject*>(contactManifold->getBody0());
    btCollisionObject* obB =
      static_cast<btCollisionObject*>(contactManifold->getBody1());

    //    contactManifold->refreshContactPoints(obA->getWorldTransform(),obB->getWorldTransform());

    unsigned int numContacts = contactManifold->getNumContacts();

    if (obA->getUserPointer())
    {

      btVector3 center;
      btScalar radius;

      obA->getCollisionShape()->getBoundingSphere(center, radius);
      double contactThreshold = radius * .1;

      /* closed (=10% radius) contact points elimination */
      unsigned int zone[4];
      unsigned int maxZone = 0;
      unsigned int bestContact[4];
      double minDistance[4];
      for (unsigned int j = 0; j < numContacts; ++j)
      {
        btManifoldPoint& cpointj = contactManifold->getContactPoint(j);
        zone[j] = maxZone;
        bestContact[zone[j]] = j;

        btVector3 posaj = cpointj.getPositionWorldOnA();
        btVector3 posbj = cpointj.getPositionWorldOnB();

        for (unsigned int k = 0; k < j; ++k)
        {

          btManifoldPoint& cpointk = contactManifold->getContactPoint(k);

          btVector3 posak = cpointk.getPositionWorldOnA();
          btVector3 posbk = cpointk.getPositionWorldOnB();
          double da = (posak - posaj).dot(posak - posaj);
          double db = (posbk - posbj).dot(posbk - posbj);

          DEBUG_PRINTF("j : %d, k : %d, da : %g,  db %g\n", j, k, da, db);

          if (da < contactThreshold || db < contactThreshold)
          {
            zone[j] = zone[k];
            bestContact[zone[j]] = bestContact[zone[k]];
          }
        }
        if (zone[j] == maxZone) ++maxZone;
      }

      assert(maxZone <= 4);

      DEBUG_PRINTF("maxZone : %d\n", maxZone);
      for (unsigned int z = 0; z < maxZone; ++z)
      {
        DEBUG_PRINTF("z=%d, bestContact[z]=%d, getContactPoint(bestContact[z]).getDistance()=%g\n", z, bestContact[z], contactManifold->getContactPoint(bestContact[z]).getDistance());
        minDistance[z] = contactManifold->getContactPoint(bestContact[z]).getDistance();
      }

      for (unsigned int j = 0; j < numContacts; ++j)
      {
        assert(zone[j] <= maxZone);

        btManifoldPoint& cpointj = contactManifold->getContactPoint(j);
        DEBUG_PRINTF("zone[j] : j = %d,  zone[j]=%d\n", j, zone[j]);
        DEBUG_PRINTF("cpointj.getDistance() = %g\n", cpointj.getDistance());
        DEBUG_PRINTF("cpointj.getDistance() > %g\n", minDistance[zone[j]]);
        if (cpointj.getDistance() < minDistance[zone[j]])
        {
          DEBUG_PRINTF("bestContact[zone[j]] = j : zone[j]=%d, j=%d\n", zone[j], j);
          bestContact[zone[j]] = j;
          DEBUG_PRINTF(" minDistance[zone[j]] = cpointj.getDistance() = %g\n",  cpointj.getDistance());
          minDistance[zone[j]] = cpointj.getDistance();
        }
      }

#ifndef NDEBUG
      for (unsigned int j = 0; j < numContacts; ++j)
      {
        assert(zone[j] < 4);
        assert(bestContact[zone[j]] < 4);

        btManifoldPoint& cpointj = contactManifold->getContactPoint(j);
        btManifoldPoint& cpointbj = contactManifold->getContactPoint(bestContact[zone[j]]);
        btVector3 posaj = cpointj.getPositionWorldOnA();
        btVector3 posbj = cpointj.getPositionWorldOnB();

        DEBUG_PRINTF("manifold %d, numContacts %d, best : %g,  contact %d : %g\n", i, numContacts, cpointbj.getDistance(), j, cpointj.getDistance());

        assert(cpointbj.getDistance() <= cpointj.getDistance());

        for (unsigned int k = 0; k < j; ++k)
        {
          btManifoldPoint& cpointk = contactManifold->getContactPoint(k);

          btVector3 posak = cpointk.getPositionWorldOnA();
          btVector3 posbk = cpointk.getPositionWorldOnB();

          if (((posak - posaj).dot(posak - posaj) < contactThreshold) ||
              (posbk - posbj).dot(posbk - posbj) < contactThreshold)
          {
            DEBUG_PRINTF("zone[j]==zone[k] ? j : %d, k : %d, zone[j : %d, zone[k] : %d\n", j, k, zone[j], zone[k]);
            //  assert (zone[j] == zone[k]);
          }
        }
      }
#endif

      for (unsigned int z = 0; z < maxZone; ++z)
      {

        SP::btManifoldPoint cpoint(createSPtrbtManifoldPoint(contactManifold->getContactPoint(bestContact[z])));
        DEBUG_PRINTF("manifold %d, contact %d, &contact %p, lifetime %d\n", i, bestContact[z], &*cpoint, cpoint->getLifeTime());

        std::map<btManifoldPoint*, bool>::iterator itc;
        itc = contactPoints.find(&*cpoint);

        if (itc == contactPoints.end())
        {
          SP::Interaction inter;
          if (_nslaw->size() == 3)
          {
            SP::BulletR rel(new BulletR(cpoint, createSPtrbtPersistentManifold(*contactManifold)));
            inter.reset(new Interaction(3, _nslaw, rel, 4 * i + z));
          }
          else
          {
            if (_nslaw->size() == 1)
            {
              SP::BulletFrom1DLocalFrameR rel(new BulletFrom1DLocalFrameR(cpoint));
              inter.reset(new Interaction(1, _nslaw, rel, 4 * i + z));
            }
          }

          // should no be mixed with something else that use UserPointer!
          assert(obA->getUserPointer());

          SP::BulletDS dsa(static_cast<BulletDS*>(obA->getUserPointer())->shared_ptr());

          if (obB->getUserPointer())
          {
            SP::BulletDS dsb(static_cast<BulletDS*>(obB->getUserPointer())->shared_ptr());

            inter->insert(dsa);
            inter->insert(dsb);
            cpoint->m_userPersistentData = &*inter;
            insertInteraction(inter);
          }
          else
          {
            inter->insert(dsa);
            cpoint->m_userPersistentData = &*inter;
            insertInteraction(inter);
          }
        }

        if (cpoint->m_userPersistentData)
        {
          activeInteractions[static_cast<Interaction *>(cpoint->m_userPersistentData)->number()] = true;
          DEBUG_PRINTF("Interaction number %d = true\n", static_cast<Interaction *>(cpoint->m_userPersistentData)->number());
          DEBUG_PRINTF("cpoint %p  = true\n", &*cpoint);
        }
        contactPoints[&*cpoint] = true;
        DEBUG_PRINTF("cpoint %p  = true\n", &*cpoint);
      }
    }
  }

  for (std::vector<int>::iterator it = gOrphanedInteractions.begin();
       it != gOrphanedInteractions.end();
       ++it)
  {
    DEBUG_PRINTF("setting contact point to false for orphaned interaction %d (was %d)\n", *it, activeInteractions[*it]);
    activeInteractions[*it] = false;
  }
  gOrphanedInteractions.clear();

  // 4. remove old contact points
  boost::tie(ui0, ui0end) = indexSet0->vertices();
  for (v0next = ui0 ;
       ui0 != ui0end; ui0 = v0next)
  {
    ++v0next;  // trick to iterate on a dynamic bgl graph
    SP::Interaction inter0 = indexSet0->bundle(*ui0);

    if (!contactPoints[&*ask<ForContactPoint>(*(inter0->relation()))])
    {

      //      assert (!contactPoints[&*ask<ForContactPoint>(*(inter0->relation()))]);

      DEBUG_PRINTF("remove contact %p, lifetime %d\n",
                   &*ask<ForContactPoint>(*(inter0->relation())),
                   ask<ForContactPoint>(*(inter0->relation()))->getLifeTime());
      model()->nonSmoothDynamicalSystem()->removeInteraction(inter0);
    }

  }

  DEBUG_PRINT("-----end build interaction\n");

  model()->simulation()->initOSNS();

}
