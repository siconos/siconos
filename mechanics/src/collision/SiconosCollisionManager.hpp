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

/*! \file SiconosCollisionManager.hpp
\brief A mechanics world is a Siconos InteractionManager that supports
  static contactors and dynamic contactors attached to special
  Dynamical Systems (BodyDS, derived from NewtonEulerDS) found in the
  NonSmoothDynamicalSystem.
*/

#ifndef SiconosCollisionManager_h
#define SiconosCollisionManager_h

#include <InteractionManager.hpp>
#include <SiconosContactor.hpp>

class SiconosCollisionManager : public InteractionManager
{
protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(SiconosCollisionManager);

public:
  SiconosCollisionManager() : InteractionManager() {}
  virtual ~SiconosCollisionManager() {}

  /** An opaque handle can be used to refer to a specific static
   * contactor set previously added to the collision manager. */
  typedef void* StaticContactorSetID;

public:
  virtual StaticContactorSetID insertStaticContactorSet(
    SP::SiconosContactorSet cs, SP::SiconosVector position = SP::SiconosVector()) = 0;

  virtual bool removeStaticContactorSet(StaticContactorSetID id) = 0;
};

#endif /* SiconosCollisionManager.hpp */
