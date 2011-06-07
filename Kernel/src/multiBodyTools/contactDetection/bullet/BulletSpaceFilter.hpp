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

#ifndef BulletSpaceFilter_hpp
#define BulletSpaceFilter_hpp

#include "BulletSiconos.hpp"
#include "SpaceFilter.hpp"
#include "BulletR.hpp"
#include "BulletRImpact.hpp"
#include "BulletDS.hpp"

#include <btBulletCollisionCommon.h>

//typedef btSimpleBroadphase BulletBroadPhase;
//typedef bt32BitAxisSweep3 BulletBroadPhase;
typedef btDbvtBroadphase BulletBroadPhase;


TYPEDEF_SPTR(BulletBroadPhase);

class BulletSpaceFilter : public SpaceFilter
{

protected:
  boost::shared_ptr<std::vector<SP::btCollisionObject> > _staticObjects;
  boost::shared_ptr<std::vector<SP::btCollisionShape> > _staticShapes;

  SP::btCollisionWorld _collisionWorld;
  SP::btDefaultCollisionConfiguration _collisionConfiguration;
  SP::btCollisionDispatcher _dispatcher;
  SP::btVector3 _worldAabbMin;
  SP::btVector3 _worldAabbMax;
  SP::BulletBroadPhase _broadphase;

public:
  BulletSpaceFilter(SP::NonSmoothDynamicalSystem nsds, SP::NonSmoothLaw nslaw,
                    SP::btVector3 aabbMin, SP::btVector3 aabbMax);

  SP::BulletBroadPhase broadphase() const
  {
    return _broadphase;
  };

  SP::btCollisionWorld collisionWorld() const
  {
    return _collisionWorld;
  };

  boost::shared_ptr<std::vector<SP::btCollisionObject> >staticObjects() const
  {
    return _staticObjects;
  };

  boost::shared_ptr<std::vector<SP::btCollisionShape> > staticShapes() const
  {
    return _staticShapes;
  };

  void buildInteractions(double);

  ACCEPT_STD_VISITORS();
};

DEFINE_SPTR(BulletSpaceFilter);


struct ForCollisionWorld : public Question<SP::btCollisionWorld>
{
  ANSWER(BulletSpaceFilter, collisionWorld());
};

struct ForStaticObjects : public Question< boost::shared_ptr<std::vector<SP::btCollisionObject> > >
{
  ANSWER(BulletSpaceFilter, staticObjects());
};

struct ForStaticShapes : public Question< boost::shared_ptr<std::vector<SP::btCollisionShape> > >
{
  ANSWER(BulletSpaceFilter, staticShapes());
};


struct ForContactPoint : public Question<SP::btManifoldPoint>
{
  ANSWER(BulletR, contactPoint());
  ANSWER(BulletRImpact, contactPoint());
};





#endif
