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

#include <Model.hpp>
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


