/* Siconos-Kernel, Copyright INRIA 2005-2010.
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

#include "BulletTimeStepping.hpp"
#include "BulletSiconos.hpp"
#include "BulletDS.hpp"
#include <btBulletCollisionCommon.h>

struct ForPosition : public Question<SP::SiconosVector>
{
  ANSWER(BulletDS, q());
};

void BulletTimeStepping::updateWorldFromDS()
{
  DynamicalSystemsGraph& dsg = *model()->nonSmoothDynamicalSystem()->dynamicalSystems();
  DynamicalSystemsGraph::VIterator dsi, dsiend;
  boost::tie(dsi, dsiend) = dsg.vertices();
  for (; dsi != dsiend; ++dsi)
  {
    btCollisionObject& co = *ask<ForCollisionObject>(*(dsg.bundle(*dsi)));
    SiconosVector& q = *ask<ForPosition>(*(dsg.bundle(*dsi)));

    co.getWorldTransform().getOrigin().setX(q(0));
    co.getWorldTransform().getOrigin().setY(q(1));
    co.getWorldTransform().getOrigin().setZ(q(2));


    assert(fabs(sqrt(pow(q(3), 2) + pow(q(4), 2) +  pow(q(5), 2) +  pow(q(6), 2)) - 1.) < 1e-10);

    co.getWorldTransform().getBasis().setRotation(btQuaternion(q(4), q(5),
        q(6), q(3)));

    co.setActivationState(ACTIVE_TAG);

  }
}


