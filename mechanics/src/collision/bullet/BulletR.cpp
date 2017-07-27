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


// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES 1
#include <debug.h>

#include "BulletR.hpp"
#include <BodyDS.hpp>
#include <Interaction.hpp>

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunreachable-code"
#pragma clang diagnostic ignored "-Woverloaded-virtual"
#elif !(__INTEL_COMPILER || __APPLE__ )
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#endif

#include <BulletCollision/NarrowPhaseCollision/btManifoldPoint.h>
#include <BulletCollision/CollisionDispatch/btCollisionObject.h>

#include <btBulletCollisionCommon.h>

#if defined(__clang__)
#pragma clang diagnostic pop
#elif !(__INTEL_COMPILER || __APPLE__ )
#pragma GCC diagnostic pop
#endif

static void copyBtVector3(const btVector3 &from, SiconosVector& to)
{
  to(0) = from.x();
  to(1) = from.y();
  to(2) = from.z();
}

// TODO: sppoint parameter used only by BulletSpaceFilter
BulletR::BulletR(const btManifoldPoint &point,
                 SP::SiconosVector q1, SP::SiconosVector q2,
                 bool flip,
                 double y_correction_A,
                 double y_correction_B,
                 double scaling) :
  NewtonEulerFrom3DLocalFrameR(),
  _y_correction_A(y_correction_A),
  _y_correction_B(y_correction_B),
  _scaling(scaling),
  _flip(flip)
{
}

#include <BlockVector.hpp>
void BulletR::computeh(double time, BlockVector& q0, SiconosVector& y)
{
  DEBUG_BEGIN("BulletR::computeh(...)\n");

  // Update contact points and distance if necessary
  NewtonEulerFrom3DLocalFrameR::computeh(time, q0, y);

  // Since Pc1 and Pc2 may have changed, _contactDistance must be updated
  btVector3 dpc((*_Pc2)(0) - (*_Pc1)(0),
                (*_Pc2)(1) - (*_Pc1)(1),
                (*_Pc2)(2) - (*_Pc1)(2));
  btVector3 nc((*_Nc)(0), (*_Nc)(1), (*_Nc)(2));
  double dist = dpc.length();

  _contactDistance = dist * (nc.dot(dpc) >= 0 ? -1 : 1);

  // Due to margins we add, objects are reported as closer than they really
  // are, so we correct by a factor.
  double correction = _y_correction_A + _y_correction_B;
  y.setValue(0, _contactDistance*_scaling + correction);

  DEBUG_PRINTF("distance : %g\n",  y.getValue(0));

  DEBUG_PRINTF("position on A : %g,%g,%g\n", (*pc1())(0), (*pc1())(1), (*pc1())(2));
  DEBUG_PRINTF("position on B : %g,%g,%g\n", (*pc2())(0), (*pc2())(1), (*pc2())(2));
  DEBUG_PRINTF("normal on B   : %g,%g,%g\n", (*nc())(0), (*nc())(1), (*nc())(2));

  DEBUG_END("BulletR::computeh(...)\n");
}

void BulletR::updateContactPoints(const btManifoldPoint& point,
                                  SP::NewtonEulerDS ds1, SP::NewtonEulerDS ds2)
{
  // Flip contact points if requested
  btVector3 posa = point.getPositionWorldOnA();
  btVector3 posb = point.getPositionWorldOnB();
  if (_flip) {
    posa = point.getPositionWorldOnB();
    posb = point.getPositionWorldOnA();
  }

  // Store distance
  _contactDistance = point.getDistance();

  // Update normal
  btVector3 n(point.m_normalWorldOnB  * (_flip ? -1 : 1));
  copyBtVector3(n, *_Nc);

  // Adjust contact point positions correspondingly along normal.  TODO: This
  // assumes same distance in each direction, i.e. same margin per object.
  posa = posa * _scaling + n * _y_correction_A;
  posb = posb * _scaling - n * _y_correction_B;

  // Update relative contact point locations.
  btQuaternion qq1((*ds1->q())(4), (*ds1->q())(5),
                   (*ds1->q())(6), (*ds1->q())(3));
  btQuaternion pq1(posa.x() - (*ds1->q())(0),
                   posa.y() - (*ds1->q())(1),
                   posa.z() - (*ds1->q())(2), 0);

  // Unrotate q1-posa vector
  pq1 = qq1.inverse() * pq1 * qq1;
  (*_relPc1)(0) = pq1.x();
  (*_relPc1)(1) = pq1.y();
  (*_relPc1)(2) = pq1.z();

  if (ds2)
  {
    btQuaternion qq2((*ds2->q())(4), (*ds2->q())(5), (*ds2->q())(6), (*ds2->q())(3));
    btQuaternion pq2(posb.x() - (*ds2->q())(0), posb.y() - (*ds2->q())(1), posb.z() - (*ds2->q())(2), 0);

    // Unrotate q2-posb vector
    pq2 = qq2.inverse() * pq2 * qq2;
    (*_relPc2)(0) = pq2.x();
    (*_relPc2)(1) = pq2.y();
    (*_relPc2)(2) = pq2.z();
  }
  else
    copyBtVector3(posb, *_relPc2);

  // Update initial contact point locations which may be modified
  // during Newton loop
  copyBtVector3(posa, *_Pc1);
  copyBtVector3(posb, *_Pc2);

  assert(!((*_Nc)(0)==0 && (*_Nc)(1)==0 && (*_Nc)(2)==0)
         && "nc = 0, problems..\n");
}
