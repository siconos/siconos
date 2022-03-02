/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

#include "Contact2dR.hpp"
#include <RigidBodyDS.hpp>
#include <Interaction.hpp>
#include <BlockVector.hpp>


// void Contact2dR::computeh(const BlockVector& q, BlockVector& z, SiconosVector& y)
// {
//   DEBUG_BEGIN("Contact2dR::computeh(...)\n");

//   // Update contact points and distance if necessary
//   Lagrangian2d2DR::computeh(q, z, y);

//   y.setValue(0, distance());


//   DEBUG_PRINTF("distance : %g \n", distance());
//   DEBUG_PRINTF("position on A : %g,%g\n", (*pc1())(0), (*pc1())(1));
//   DEBUG_PRINTF("position on B : %g,%g\n", (*pc2())(0), (*pc2())(1));
//   DEBUG_PRINTF("normal on B   : %g,%g\n", (*nc())(0), (*nc())(1));

//   DEBUG_END("Contact2dR::computeh(...)\n");
// }

// void Contact2dR::updateContactPoints(const SiconosVector& pos1,
//                                      const SiconosVector& pos2,
//                                      const SiconosVector& normal)
// {
//   // Copy relative positions
//   *_relPc1 = pos1;
//   *_relPc2 = pos2;

//   // Update normal
//   *_relNc = normal;

//   assert(!((*_relNc)(0)==0 && (*_relNc)(1)==0)
//          && "nc = 0, problems..\n");
// }
void Contact2dR::updateContactPointsInAbsoluteFrame(const SiconosVector& pos1,
                                                    const SiconosVector& pos2,
                                                    const SiconosVector& normal)
{
  // Copy positions
  *_Pc1 = pos1;
  *_Pc2 = pos2;

  // Update normal
  *_Nc = normal;

  assert(!((*_Nc)(0)==0 && (*_Nc)(1)==0)
         && "nc = 0, problems..\n");
}
