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

#ifndef BulletTimeStepping_hpp
#define BulletTimeStepping_hpp

#include "BulletSiconosFwd.hpp"
#include "TimeStepping.hpp"

class BulletTimeStepping : public TimeStepping
{

public:
  BulletTimeStepping(SP::TimeDiscretisation td,
		     SP::OneStepIntegrator osi = SP::OneStepIntegrator(),
		     SP::OneStepNSProblem osnspb = SP::OneStepNSProblem()) :
    TimeStepping(td, osi, osnspb) {};

  void updateWorldFromDS();
};

#endif

