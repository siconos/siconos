/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
// #include "siconos_debug.h"

#include "Contact2d3DR.hpp"

#include "SiconosVector.hpp"

// void Contact2d3DR::computeh(const siconos::algebra::BlockVector& q,
// siconos::algebra::BlockVector& z, siconos::algebra::SiconosVector& y)
// {
//   DEBUG_BEGIN("Contact2d3DR::computeh(...)\n");

//   // Update contact points and distance if necessary
//   Lagrangian2dVelocityAngularVelocityR::computeh(q, z, y);

//   y(0) = distance();

//   DEBUG_PRINTF("distance : %g \n", distance());
//   DEBUG_PRINTF("position on A : %g,%g\n", (*pc1())(0), (*pc1())(1));
//   DEBUG_PRINTF("position on B : %g,%g\n", (*pc2())(0), (*pc2())(1));
//   DEBUG_PRINTF("normal on B   : %g,%g\n", (*nc())(0), (*nc())(1));

//   DEBUG_END("Contact2d3DR::computeh(...)\n");
// }

// void Contact2d3DR::updateContactPoints(const siconos::algebra::SiconosVector& pos1,
//                                        const siconos::algebra::SiconosVector& pos2,
//                                        const siconos::algebra::SiconosVector& normal)
// {
//   // Copy relative positions
//   *_relPc1 = pos1;
//   *_relPc2 = pos2;

//   // Update normal
//   *_relNc = normal;

//   assert(!((*_relNc)(0)==0 && (*_relNc)(1)==0)
//          && "nc = 0, problems..\n");
// }
void siconos::collision::Contact2d3DR::updateContactPointsInAbsoluteFrame(
    const siconos::algebra::SiconosVector2& pos1, const siconos::algebra::SiconosVector2& pos2,
    const siconos::algebra::SiconosVector2& normal) {
  // Copy  positions
  contactPoint1_ = pos1;
  contactPoint2_ = pos2;

  // Update normal
  nc_ = normal;

  assert(!(nc_(0) == 0 && nc_(1) == 0) && "nc = 0, problems..\n");
}
