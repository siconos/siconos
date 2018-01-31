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

/*! \file SiconosBodies.hpp
  \brief SiconosBodies class - model + plans + space filter
*/
#ifndef SiconosBodies_hpp
#define SiconosBodies_hpp
#include "Simulation.hpp"
#include "MechanicsFwd.hpp"
#include <SiconosFwd.hpp>
#include <SiconosSerialization.hpp>

/** SiconosBodies : a Siconos Model, some plans and space filtering capabilities
 */

class SiconosBodies
{

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(SiconosBodies);


  SP::FMatrix _moving_plans;
  SP::SiconosMatrix _plans;
  SP::Simulation _sim;
  SP::SpaceFilter _playground;

public:

  virtual void init() = 0;

  virtual void compute();

  SP::Simulation simulation()
  {
    return _sim;
  }


  SP::FMatrix movingPlans()
  {
    return _moving_plans;
  }
  SP::SiconosMatrix plans()
  {
    return _plans;
  }


  SP::SpaceFilter spaceFilter()
  {
    return _playground;
  };

  /** destructor
   */
  virtual ~SiconosBodies() {};

};

#endif // SiconosBodies_hpp
