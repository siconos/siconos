/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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


// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES 1
#include "siconos_debug.h"

#include "ContactR.hpp"
#include <RigidBodyDS.hpp>
#include <Interaction.hpp>
#include <BlockVector.hpp>

ContactR::ContactR()
  : NewtonEuler3DR()
{
}

void ContactR::computeh(double time, const BlockVector& q0, SiconosVector& y)
{
  DEBUG_BEGIN("ContactR::computeh(...)\n");

  // Update contact points and distance if necessary
  NewtonEuler3DR::computeh(time, q0, y);

  y.setValue(0, distance());

  DEBUG_PRINTF("position on A : %g,%g,%g\n", (*pc1())(0), (*pc1())(1), (*pc1())(2));
  DEBUG_PRINTF("position on B : %g,%g,%g\n", (*pc2())(0), (*pc2())(1), (*pc2())(2));
  DEBUG_PRINTF("normal on B   : %g,%g,%g\n", (*nc())(0), (*nc())(1), (*nc())(2));

  DEBUG_END("ContactR::computeh(...)\n");
}

void ContactR::updateContactPoints(const SiconosVector& pos1,
                                   const SiconosVector& pos2,
                                   const SiconosVector& normal)
{
  // Copy relative positions
  *_relPc1 = pos1;
  *_relPc2 = pos2;

  // Update normal
  *_relNc = normal;

  assert(!((*_relNc)(0)==0 && (*_relNc)(1)==0 && (*_relNc)(2)==0)
         && "nc = 0, problems..\n");
}
