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

#include "BulletDS.hpp"

#include <bullet/btBulletCollisionCommon.h>
#include <BulletCollision/Gimpact/btGImpactCollisionAlgorithm.h>

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
    btVector3 posa = cp.getPositionWorldOnA();
    btVector3 posb = cp.getPositionWorldOnB();

    DEBUG_PRINTF("gContactProcessedCallback : process contactPoint %p \n  position on A : (%g,%g,%g)\n  position on B : (%g,%g,%g)\n  distance : %g\n", &cp, posa[0], posa[1], posa[2], posb[0], posb[1], posb[2], cp.getDistance());
  }
  return true;
}


struct ForPosition : public Question<SP::SiconosVector>
{
  using SiconosVisitor::visit;

  ANSWER(BulletDS, q());
};

struct ForContactNum : public Question<unsigned int>
{
  using SiconosVisitor::visit;

  ANSWER(BulletR, contactNum());
};



BulletSpaceFilter::BulletSpaceFilter(SP::Model model,
                                     SP::NonSmoothLaw nslaw) :
  SpaceFilter(),
  _dynamicCollisionsObjectsInserted(false),
  _staticCollisionsObjectsInserted(false),
  _closeContactsThreshold(0.)
{

  _model = model;
  _nslaw = nslaw;

  _staticObjects.reset(new std::vector<SP::btCollisionObject>());
  _staticShapes.reset(new std::vector<SP::btCollisionShape>());

  _collisionConfiguration.reset(new btDefaultCollisionConfiguration());

//  _collisionConfiguration->setConvexConvexMultipointIterations();
//  _collisionConfiguration->setPlaneConvexMultipointIterations();

  _dispatcher.reset(new btCollisionDispatcher(&*_collisionConfiguration));

  _broadphase.reset(new btDbvtBroadphase());

  _collisionWorld.reset(new btCollisionWorld(&*_dispatcher, &*_broadphase, &*_collisionConfiguration));

  btGImpactCollisionAlgorithm::registerAlgorithm(&*_dispatcher);

  _collisionWorld->getDispatchInfo().m_useContinuous = false;

  gContactBreakingThreshold = 1.0;
//  gContactProcessedCallback = &contactProcess;
//  gContactDestroyedCallback = &contactClear;

}

void BulletSpaceFilter::buildInteractions(double time)
{

  if (! _dynamicCollisionsObjectsInserted)
  {
    DynamicalSystemsGraph& dsg = *(_model->nonSmoothDynamicalSystem()->dynamicalSystems());
    DynamicalSystemsGraph::VIterator dsi, dsiend;
    std11::tie(dsi, dsiend) = dsg.vertices();
    for (; dsi != dsiend; ++dsi)
    {
      _collisionWorld->addCollisionObject(&*(ask<ForCollisionObject>(*(dsg.bundle(*dsi)))));
    }
    
    _dynamicCollisionsObjectsInserted = true;
  }

  if (! _staticCollisionsObjectsInserted)
  {
    for(std::vector<SP::btCollisionObject>::iterator 
          ic = _staticObjects->begin(); ic != _staticObjects->end(); ++ic)
    {
      _collisionWorld->addCollisionObject((*ic).get());
    }

    _staticCollisionsObjectsInserted = true;
  }

  DEBUG_PRINT("-----start build interaction\n");

  // 1. perform broadphase bullet collision detection
  _collisionWorld->performDiscreteCollisionDetection();

  // 2. collect old interactions from Siconos graph
  std::map<btPersistentManifold*, bool> contactManifolds;
  std::map<btManifoldPoint*, bool> contactPoints;


  SP::InteractionsGraph indexSet0 = model()->nonSmoothDynamicalSystem()->topology()->indexSet(0);
  InteractionsGraph::VIterator ui0, ui0end, v0next;
  for (std11::tie(ui0,ui0end) = indexSet0->vertices(); ui0!=ui0end ; ++ui0)
  {
    Relation& relation = *(indexSet0->bundle(*ui0)->relation());
    btPersistentManifold& manifold = *ask<ForContactManifold>(relation);
    unsigned int contactNum = ask<ForContactNum>(relation);

    contactManifolds[&manifold] = false;
    contactPoints[&manifold.getContactPoint(contactNum)] = false;
  }


  unsigned int numManifolds =
    _collisionWorld->getDispatcher()->getNumManifolds();

  for (unsigned int i = 0; i < numManifolds; ++i)
  {

    btPersistentManifold* contactManifold =
      _collisionWorld->getDispatcher()->getManifoldByIndexInternal(i);

    for (unsigned int j = 0; j < contactManifold->getNumContacts(); ++j)
    {
      btManifoldPoint& cpoint = contactManifold->getContactPoint(j);

      if (contactPoints.find(&cpoint) == contactPoints.end())
      {
        assert(_nslaw->size() == 3); // 1D not yet
        
        SP::BulletR rel(new BulletR(j, createSPtrbtPersistentManifold(*contactManifold), _collisionWorld));
        SP::Interaction inter(new Interaction(3, _nslaw, rel));

        const btCollisionObject* obA =
          static_cast<const btCollisionObject*>(contactManifold->getBody0());
        const btCollisionObject* obB =
          static_cast<const btCollisionObject*>(contactManifold->getBody1());

        assert(obA->getUserPointer());
        SP::BulletDS dsa(static_cast<BulletDS*>(obA->getUserPointer())->shared_ptr());
        
        if (obB->getUserPointer())
        {
          SP::BulletDS dsb(static_cast<BulletDS*>(obB->getUserPointer())->shared_ptr());
          link(inter, dsa, dsb);
        }
        else
        {
          link(inter, dsa);
        }
        
      }
      contactPoints[&cpoint] = true;
    }    
  }
  
  // remove old interactions
  std11::tie(ui0, ui0end) = indexSet0->vertices();
  for (v0next = ui0 ;
       ui0 != ui0end; ui0 = v0next)
  {
    SP::Interaction inter = indexSet0->bundle(*ui0);
    ++v0next;  // trick to iterate on a dynamic bgl graph

    Relation& relation = *inter->relation();
    btPersistentManifold& manifold = *ask<ForContactManifold>(relation);
    unsigned int contactNum = ask<ForContactNum>(relation);
    
    if (! contactPoints[&manifold.getContactPoint(contactNum)])
    {
      model()->nonSmoothDynamicalSystem()->removeInteraction(inter);
    }

  }
  DEBUG_PRINT("-----end build interaction\n");

  model()->simulation()->initOSNS();

}
