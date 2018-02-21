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

#ifndef BulletTimeSteppingDirectProjection_hpp
#define BulletTimeSteppingDirectProjection_hpp

#include "MechanicsFwd.hpp"
#include "TimeSteppingDirectProjection.hpp"
#include "BulletSpaceFilter.hpp"

class BulletTimeSteppingDirectProjection : public TimeSteppingDirectProjection
{

  SP::BulletSpaceFilter _spaceFilter;

public:
  BulletTimeSteppingDirectProjection(
    SP::NonSmoothDynamicalSystem nsds,
    SP::BulletSpaceFilter sf,
    SP::TimeDiscretisation t,
    SP::OneStepIntegrator osi,
    SP::OneStepNSProblem osnspb_velo,
    SP::OneStepNSProblem osnspb_pos,
    unsigned int level = 1) :
    TimeSteppingDirectProjection(nsds,t, osi, osnspb_velo, osnspb_pos),
    _spaceFilter(sf) {};


  void updateWorldFromDS();
};

TYPEDEF_SPTR(BulletTimeSteppingDirectProjection)
#endif


