/* Siconos-Kernel, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/


// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES 1
#include <debug.h>

#include "BulletR.hpp"

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

#include <Interaction.hpp>

BulletR::BulletR(SP::btManifoldPoint point, double y_correction, double scaling) :
  NewtonEulerFrom3DLocalFrameR(),
  _contactPoints(point),
  _y_correction(y_correction),
  _scaling(scaling)
{
  btVector3 posa = _contactPoints->getPositionWorldOnA();
  btVector3 posb = _contactPoints->getPositionWorldOnB();

  (*pc1())(0) = posa[0];
  (*pc1())(1) = posa[1];
  (*pc1())(2) = posa[2];
  (*pc2())(0) = posb[0];
  (*pc2())(1) = posb[1];
  (*pc2())(2) = posb[2];

  (*nc())(0) = _contactPoints->m_normalWorldOnB[0];
  (*nc())(1) = _contactPoints->m_normalWorldOnB[1];
  (*nc())(2) = _contactPoints->m_normalWorldOnB[2];

  assert(!((*nc())(0)==0 && (*nc())(1)==0 && (*nc())(2)==0)
         && "nc = 0, problems..\n");
}

void BulletR::computeh(double time, BlockVector& q0, SiconosVector& y)
{
  DEBUG_BEGIN("BulletR::computeh(...)\n");

  NewtonEulerR::computeh(time, q0, y);

  DEBUG_PRINT("start of computeh\n");

  btVector3 posa = _contactPoints->getPositionWorldOnA();
  btVector3 posb = _contactPoints->getPositionWorldOnB();

  (*pc1())(0) = posa[0]*_scaling;
  (*pc1())(1) = posa[1]*_scaling;
  (*pc1())(2) = posa[2]*_scaling;
  (*pc2())(0) = posb[0]*_scaling;
  (*pc2())(1) = posb[1]*_scaling;
  (*pc2())(2) = posb[2]*_scaling;

  {
    // Due to margins we add, objects are reported as closer than they really
    // are, so we correct by a factor.
    y.setValue(0, (_contactPoints->getDistance() + _y_correction)*_scaling);

    (*nc())(0) = _contactPoints->m_normalWorldOnB[0];
    (*nc())(1) = _contactPoints->m_normalWorldOnB[1];
    (*nc())(2) = _contactPoints->m_normalWorldOnB[2];

    // Adjust contact point positions correspondingly along normal.  TODO: This
    // assumes same distance in each direction, i.e. same margin per object.
    // (*pc1())(0) -= (*nc())(0) * _y_correction/2;
    // (*pc1())(1) -= (*nc())(1) * _y_correction/2;
    // (*pc1())(2) -= (*nc())(2) * _y_correction/2;
    // (*pc2())(0) += (*nc())(0) * _y_correction/2;
    // (*pc2())(1) += (*nc())(1) * _y_correction/2;
    // (*pc2())(2) += (*nc())(2) * _y_correction/2;
  }

  DEBUG_PRINTF("distance : %g\n",  y.getValue(0));


  DEBUG_PRINTF("position on A : %g,%g,%g\n", posa[0], posa[1], posa[2]);
  DEBUG_PRINTF("position on B : %g,%g,%g\n", posb[0], posb[1], posb[2]);
  DEBUG_PRINTF("normal on B   : %g,%g,%g\n", (*nc())(0), (*nc())(1), (*nc())(2));

  DEBUG_END("BulletR::computeh(...)\n");


}
