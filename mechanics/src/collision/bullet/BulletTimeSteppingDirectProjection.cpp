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

#include "BulletTimeSteppingDirectProjection.hpp"
#include "BulletSiconosFwd.hpp"
#include "BulletDS.hpp"

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunreachable-code"
#pragma clang diagnostic ignored "-Woverloaded-virtual"
#elif !(__INTEL_COMPILER || __APPLE__ )
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#endif


// looks like bullet have overloaded virtual functions ...

#include <btBulletCollisionCommon.h>

#if defined(__clang__)
#pragma clang diagnostic pop
#elif !(__INTEL_COMPILER || __APPLE__ )
#pragma GCC diagnostic pop
#endif

#ifndef DEBUG_BULLET_TIMESTEPPING
#define DEBUG_MESSAGES 1
#endif

#include <debug.h>

#include <NonSmoothDynamicalSystem.hpp>

void BulletTimeSteppingDirectProjection::updateWorldFromDS()
{
  DynamicalSystemsGraph& dsg = *_nsds->dynamicalSystems();
  DynamicalSystemsGraph::VIterator dsi, dsiend;
  std11::tie(dsi, dsiend) = dsg.vertices();

  static UpdateCollisionObjects up;

  for (; dsi != dsiend; ++dsi)
  {
    dsg.bundle(*dsi)->accept(up);
  }

}


