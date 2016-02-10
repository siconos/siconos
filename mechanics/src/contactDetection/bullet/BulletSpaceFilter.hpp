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

#include "BulletSiconosFwd.hpp"
#include "SpaceFilter.hpp"

class BulletSpaceFilter : public SpaceFilter
{

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(BulletSpaceFilter);

  SP::StaticObjects _staticObjects;

  SP::btCollisionWorld _collisionWorld;
  SP::btDefaultCollisionConfiguration _collisionConfiguration;
  SP::btCollisionDispatcher _dispatcher;
  SP::btBroadphaseInterface _broadphase;
  bool _dynamicCollisionsObjectsInserted;
  bool _staticCollisionsObjectsInserted;

  double _closeContactsThreshold;

public:
  BulletSpaceFilter(SP::Model model);

  /** get Bullet collision configuration
      \return a pointer on a Bullet collision configuration
  */
  SP::btDefaultCollisionConfiguration collisionConfiguration() const
  {
    return _collisionConfiguration;
  }

  /** get Bullet broadphase
      \return pointer on a BulletBroadPhase
  */
  SP::btBroadphaseInterface broadphase() const
  {
    return _broadphase;
  };

  /** get bullet collision world
      \return pointer on btCollisionWorld
  */
  SP::btCollisionWorld collisionWorld() const
  {
    return _collisionWorld;
  };

  /** get static objects
      \return a std::vector of btCollisionObject
  */
  SP::StaticObjects staticObjects() const
  {
    return _staticObjects;
  };

  /** add a static object
   * \param co a btCollisionObject
   * \param id contactor id of the collision object for non smooth law selection
   */
  void addStaticObject(SP::btCollisionObject co, unsigned int id);

  /** dynamically add a bullet dynamical system
   * \param ds a BulletDS to be added
   * \param simulation a simulation to which ds shall belong
   * \param osi a OneStepIntegrator to use, or if NULL, simulation's first OSI is used.
   */
  void addDynamicObject(SP::BulletDS ds,
                        SP::Simulation simulation,
                        SP::OneStepIntegrator osi = NULL);

  /** execute the broadphase contact detection and build indexSet0
   */
  void buildInteractions(double);


  /** set close contact parameter
   *  \param threshold double value that will be multiplicated by the
   *         radius of the object bouncing box
  */
  void setCloseContactFilterParam(double threshold)
  {
    _closeContactsThreshold = threshold;
  }

  ACCEPT_STD_VISITORS();

  /** set a new collision configuration
   * \param collisionConfig the new bullet collision configuration
   */
  void setCollisionConfiguration(
    SP::btDefaultCollisionConfiguration collisionConfig);


};


#include <Question.hpp>
#include "BulletR.hpp"
#include "BulletFrom1DLocalFrameR.hpp"

struct ForCollisionWorld : public Question<SP::btCollisionWorld>
{
  ANSWER(BulletSpaceFilter, collisionWorld());
};





struct ForContactPoint : public Question<SP::btManifoldPoint>
{
  ANSWER(BulletR, contactPoint());
  ANSWER_NOUSING(BulletFrom1DLocalFrameR, contactPoint());
  void visit(const NewtonEulerR&)
  {
  }
};

#endif
