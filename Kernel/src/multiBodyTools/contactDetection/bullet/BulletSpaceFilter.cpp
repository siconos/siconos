/* Siconos-Kernel, Copyright INRIA 2005-2010.
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

#include "BulletR.hpp"
#include "BulletDS.hpp"


struct ForPosition : public Question<SP::SiconosVector>
{
  ANSWER(BulletDS, q());
};

BulletSpaceFilter::BulletSpaceFilter(SP::NonSmoothDynamicalSystem nsds,
                                     SP::NonSmoothLaw nslaw) : SpaceFilter()
{

  _nsds = nsds;
  _nslaw = nslaw;

  _staticObjects.reset(new std::vector<SP::btCollisionObject>());
  _staticShapes.reset(new std::vector<SP::btCollisionShape>());

  _collisionConfiguration.reset(new btDefaultCollisionConfiguration());
  _dispatcher.reset(new btCollisionDispatcher(&*_collisionConfiguration));

  _worldAabbMin.reset(new btVector3(-1000, -1000, -1000));
  _worldAabbMax.reset(new btVector3(1000, 1000, 1000));

  _broadphase.reset(new btAxisSweep3(*_worldAabbMin, *_worldAabbMax));

  _collisionWorld.reset(new btCollisionWorld(&*_dispatcher, &*_broadphase, &*_collisionConfiguration));

  SP::DynamicalSystemsGraph dsg = _nsds->dynamicalSystems();
  DynamicalSystemsGraph::VIterator dsi, dsiend;
  boost::tie(dsi, dsiend) = dsg->vertices();
  for (; dsi != dsiend; ++dsi)
  {
    _collisionWorld->addCollisionObject(&*(ask<ForCollisionObject>(*(dsg->bundle(*dsi)))));
  };


};

void BulletSpaceFilter::buildInteractions(double time)
{

  // set bullet coordinates from Siconos coordinates
  SP::DynamicalSystemsGraph dsg = _nsds->dynamicalSystems();
  DynamicalSystemsGraph::VIterator dsi, dsiend;
  boost::tie(dsi, dsiend) = dsg->vertices();
  for (; dsi != dsiend; ++dsi)
  {
    SP::btCollisionObject co = ask<ForCollisionObject>(*(dsg->bundle(*dsi)));
    SP::SiconosVector q = ask<ForPosition>(*(dsg->bundle(*dsi)));

    co->getWorldTransform().getOrigin().setX((*q)(0));
    co->getWorldTransform().getOrigin().setY((*q)(1));
    co->getWorldTransform().getOrigin().setZ((*q)(2));

    co->getWorldTransform().getBasis().setRotation(btQuaternion((*q)(4), (*q)(5),
        (*q)(6), (*q)(3)));
  };


  _collisionWorld->performDiscreteCollisionDetection();

  std::map<btManifoldPoint*, bool> contactPoints;


  SP::UnitaryRelationsGraph indexSet0 = _nsds->topology()->indexSet(0);
  UnitaryRelationsGraph::VIterator ui0, ui0end, v0next;
  boost::tie(ui0, ui0end) = indexSet0->vertices();
  for (v0next = ui0 ;
       ui0 != ui0end; ui0 = v0next)
  {
    UnitaryRelation& ur0 = *(indexSet0->bundle(*ui0));
    ++v0next;  // trick to iterate on a dynamic bgl graph
    contactPoints[&*ask<ForContactPoints>(*(ur0.interaction()->relation()))] = false;
  };

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

    unsigned int numContacts = contactManifold->getNumContacts();

    for (unsigned int j = 0; j < numContacts; ++j)
    {
      SP::btManifoldPoint cpoint(createSPtrbtManifoldPoint(contactManifold->getContactPoint(j)));

      std::map<btManifoldPoint*, bool>::iterator itc;
      itc = contactPoints.find(&*cpoint);

      if (itc == contactPoints.end())
      {
        SP::BulletR rel(new BulletR(cpoint));
        SP::Interaction inter(new Interaction(3, _nslaw, rel, 4 * i + j));

        // should no be mixed with something else that use UserPointer!
        if (obA->getUserPointer())
        {
          SP::BulletDS dsa(static_cast<BulletDS*>(obA->getUserPointer())->shared_ptr());
          inter->insert(dsa);
        }
        if (obB->getUserPointer())
        {
          SP::BulletDS dsb(static_cast<BulletDS*>(obB->getUserPointer())->shared_ptr());
          inter->insert(dsb);
        }
        _nsds->topology()->insertInteraction(inter);

      }
      contactPoints[&*cpoint] = true;
    }

    boost::tie(ui0, ui0end) = indexSet0->vertices();
    for (v0next = ui0 ;
         ui0 != ui0end; ui0 = v0next)
    {
      ++v0next;  // trick to iterate on a dynamic bgl graph
      UnitaryRelation& ur0 = *(indexSet0->bundle(*ui0));

      if (!contactPoints[&*ask<ForContactPoints>(*(ur0.interaction()->relation()))])
      {
        _nsds->removeInteraction(ur0.interaction());
      }

    }
  }

  _nsds->topology()->computeRelativeDegrees();

}
