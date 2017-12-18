/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2017 INRIA.
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
#include <debug.h>

#include "ContactR.hpp"
#include <BodyDS.hpp>
#include <Interaction.hpp>
#include <BlockVector.hpp>

#include <boost/math/quaternion.hpp>

static
void copyQuatRot(boost::math::quaternion<double>& from, SiconosVector& to)
{
  to(3) = from.R_component_1();
  to(4) = from.R_component_2();
  to(5) = from.R_component_3();
  to(6) = from.R_component_4();
}

static
void copyQuatRot(const SiconosVector& from, boost::math::quaternion<double>& to)
{
  to = boost::math::quaternion<double>(from(3), from(4), from(5), from(6));
}

static
void copyQuatPos(const boost::math::quaternion<double>& from, SiconosVector& to)
{
  to(0) = from.R_component_2();
  to(1) = from.R_component_3();
  to(2) = from.R_component_4();
}

static
void copyQuatPos(const SiconosVector& from, boost::math::quaternion<double>& to)
{
  to = boost::math::quaternion<double>(0, from(0), from(1), from(2));
}

ContactR::ContactR(SP::SiconosVector q1, SP::SiconosVector q2, bool flip,
                   double y_correction_A, double y_correction_B, double scaling)
  : NewtonEulerFrom3DLocalFrameR(),
    _y_correction_A(y_correction_A),
    _y_correction_B(y_correction_B),
    _scaling(scaling),
    _flip(flip)
{
}

void ContactR::computeh(double time, BlockVector& q0, SiconosVector& y)
{
  DEBUG_BEGIN("ContactR::computeh(...)\n");

  // Update contact points and distance if necessary
  NewtonEulerFrom3DLocalFrameR::computeh(time, q0, y);

  // Since Pc1 and Pc2 may have changed, _contactDistance must be updated
  SiconosVector dpc(*_Pc2 - *_Pc1);
  double dist = dpc.norm2();

  _contactDistance = dist * (inner_prod(*_Nc, dpc) >= 0 ? -1 : 1);

  // Due to margins we add, objects are reported as closer than they really
  // are, so we correct by a factor.
  double correction = _y_correction_A + _y_correction_B;
  y.setValue(0, _contactDistance*_scaling + correction);

  DEBUG_PRINTF("distance : %g\n",  y.getValue(0));

  DEBUG_PRINTF("position on A : %g,%g,%g\n", (*pc1())(0), (*pc1())(1), (*pc1())(2));
  DEBUG_PRINTF("position on B : %g,%g,%g\n", (*pc2())(0), (*pc2())(1), (*pc2())(2));
  DEBUG_PRINTF("normal on B   : %g,%g,%g\n", (*nc())(0), (*nc())(1), (*nc())(2));

  DEBUG_END("ContactR::computeh(...)\n");
}

void ContactR::updateContactPoints(SiconosVector& pos1, SiconosVector& pos2,
                                   double distance, SiconosVector& normal,
                                   SP::NewtonEulerDS ds1,
                                   SP::NewtonEulerDS ds2)
{
  // Flip contact points if requested
  SiconosVector posa(pos1);
  SiconosVector posb(pos2);
  if (_flip) {
    SiconosVector posa = pos2;
    SiconosVector posb = pos1;
  }

  // Store distance
  _contactDistance = distance;

  // Update normal
  *_Nc = normal * (_flip ? -1 : 1);

  // Adjust contact point positions correspondingly along normal.
  // (The correction factor is split in two in case objects have
  // different margins.)
  posa = posa * _scaling + normal * _y_correction_A;
  posb = posb * _scaling - normal * _y_correction_B;

  // Update relative contact point locations.
  boost::math::quaternion<double> qq1, pq1, pp1;
  copyQuatRot(*ds1->q(), qq1);
  copyQuatPos(*ds1->q(), pq1);
  copyQuatPos(posa, pp1);

  // Unrotate q1-posa vector
  pq1 = (1.0 / qq1) * (pp1-pq1) * qq1;
  copyQuatPos(pq1, *_relPc1);

  if (ds2)
  {
    boost::math::quaternion<double> qq2, pq2, pp2;
    copyQuatRot(*ds2->q(), qq2);
    copyQuatPos(*ds2->q(), pq2);
    copyQuatPos(posb, pp2);

    // Unrotate q2-posb vector
    pq2 = (1.0 / qq2) * (pp2-pq2) * qq2;
    copyQuatPos(pq2, *_relPc2);
  }
  else
    *_relPc2 = posb;

  // Update initial contact point locations which may be modified
  // during Newton loop
  *_Pc1 = posa;
  *_Pc2 = posb;

  assert(!((*_Nc)(0)==0 && (*_Nc)(1)==0 && (*_Nc)(2)==0)
         && "nc = 0, problems..\n");
}
