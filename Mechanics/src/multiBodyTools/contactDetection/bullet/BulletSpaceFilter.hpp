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

#ifndef BulletSpaceFilter_hpp
#define BulletSpaceFilter_hpp

#include "BulletSiconos.hpp"
#include "SpaceFilter.hpp"
#include "BulletR.hpp"
#include "BulletFrom1DLocalFrameR.hpp"
#include "BulletDS.hpp"

#include <bullet/btBulletCollisionCommon.h>

//typedef btSimpleBroadphase BulletBroadPhase;
typedef bt32BitAxisSweep3 BulletBroadPhase;
//typedef btDbvtBroadphase BulletBroadPhase;


TYPEDEF_SPTR(BulletBroadPhase)

class BulletSpaceFilter : public SpaceFilter
{

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(BulletSpaceFilter);

  std11::shared_ptr<std::vector<SP::btCollisionObject> > _staticObjects;
  std11::shared_ptr<std::vector<SP::btCollisionShape> > _staticShapes;

  SP::btCollisionWorld _collisionWorld;
  SP::btDefaultCollisionConfiguration _collisionConfiguration;
  SP::btCollisionDispatcher _dispatcher;
  SP::btVector3 _worldAabbMin;
  SP::btVector3 _worldAabbMax;
  SP::BulletBroadPhase _broadphase;

public:
  BulletSpaceFilter(SP::Model model, SP::NonSmoothLaw nslaw,
                    SP::btVector3 aabbMin, SP::btVector3 aabbMax);

  SP::BulletBroadPhase broadphase() const
  {
    return _broadphase;
  };

  SP::btCollisionWorld collisionWorld() const
  {
    return _collisionWorld;
  };

  std11::shared_ptr<std::vector<SP::btCollisionObject> >staticObjects() const
  {
    return _staticObjects;
  };

  std11::shared_ptr<std::vector<SP::btCollisionShape> > staticShapes() const
  {
    return _staticShapes;
  };

  void buildInteractions(double);

  ACCEPT_STD_VISITORS();
};

DEFINE_SPTR(BulletSpaceFilter)


struct ForCollisionWorld : public Question<SP::btCollisionWorld>
{
  using SiconosVisitor::visit;

  ANSWER(BulletSpaceFilter, collisionWorld());
};

struct ForStaticObjects : public Question< std11::shared_ptr<std::vector<SP::btCollisionObject> > >
{
  using SiconosVisitor::visit;

  ANSWER(BulletSpaceFilter, staticObjects());
};

struct ForStaticShapes : public Question< std11::shared_ptr<std::vector<SP::btCollisionShape> > >
{
  using SiconosVisitor::visit;

  ANSWER(BulletSpaceFilter, staticShapes());
};


struct ForContactPoint : public Question<SP::btManifoldPoint>
{
  using SiconosVisitor::visit;

  ANSWER(BulletR, contactPoint());
  ANSWER(BulletFrom1DLocalFrameR, contactPoint());
};





#endif
