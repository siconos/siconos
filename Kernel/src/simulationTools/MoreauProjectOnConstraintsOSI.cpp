/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
#include "MoreauProjectOnConstraintsOSI.hpp"
#include "Simulation.hpp"
#include "Model.hpp"
#include "NewtonEulerDS.hpp"
#include "LagrangianDS.hpp"


//#define DEBUG_MESSAGES
//#define DEBUG_WHERE_MESSAGES
#include <debug.h>


void MoreauProjectOnConstraintsOSI::initialize()
{

  Moreau::initialize();

  ConstDSIterator itDS;
  for (itDS = OSIDynamicalSystems->begin(); itDS != OSIDynamicalSystems->end(); ++itDS)
  {
    Type::Siconos dsType = Type::value(**itDS);
    if (dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
    {
      SP::LagrangianDS d = boost::static_pointer_cast<LagrangianDS> (*itDS);
      d->allocateWorkVector(DynamicalSystem::qtmp, d->getNdof());
    }
    else if (dsType == Type::NewtonEulerDS)
    {
      SP::NewtonEulerDS d = boost::static_pointer_cast<NewtonEulerDS>(*itDS);
      d->allocateWorkVector(DynamicalSystem::qtmp, d->q()->size());
    }
    else
    {
      RuntimeException::selfThrow("MoreauProjectOnConstraintsOSI::initialize() - DS not of the right type");
    }
  }
}

