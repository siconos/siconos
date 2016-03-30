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

#include "MechanicsFwd.hpp"
#include "BodyDS.hpp"
#include "BodyTimeStepping.hpp"
#include "SiconosContactor.hpp"

#include <Model.hpp>
#include <NonSmoothDynamicalSystem.hpp>

struct CollisionUpdater : public SiconosVisitor
{
  using SiconosVisitor::visit;

  void visit(const BodyDS& bds)
  {
    bds.contactor()->setPosition(bds.q());
  }
};

void BodyTimeStepping::updateWorldFromDS()
{
  DynamicalSystemsGraph& dsg = *model()->nonSmoothDynamicalSystem()->dynamicalSystems();
  DynamicalSystemsGraph::VIterator dsi, dsiend;
  std11::tie(dsi, dsiend) = dsg.vertices();

  CollisionUpdater up;

  for (; dsi != dsiend; ++dsi)
  {
    dsg.bundle(*dsi)->accept(up);
  }
}
