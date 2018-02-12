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
  BulletSpaceFilter(SP::NonSmoothDynamicalSystem nsds);

  /** Destructor.
   */
  virtual ~BulletSpaceFilter();

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
                        SP::OneStepIntegrator osi = std11::shared_ptr<OneStepIntegrator>());

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

int btScalarSize();

#endif
