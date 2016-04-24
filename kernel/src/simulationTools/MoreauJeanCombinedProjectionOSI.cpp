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
#include "MoreauJeanCombinedProjectionOSI.hpp"
#include "Simulation.hpp"
#include "LagrangianDS.hpp"
#include "NewtonEulerDS.hpp"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
//#define DEBUG_WHERE_MESSAGES
#include <debug.h>

void MoreauJeanCombinedProjectionOSI::initialize()
{

  MoreauJeanOSI::initialize();
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if (!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    Type::Siconos dsType = Type::value(*ds);

    if (dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
    {
      SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (ds);
      d->allocateWorkVector(DynamicalSystem::qtmp, d->getNdof());
    }
    else if (dsType == Type::NewtonEulerDS)
    {
      SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS>(ds);
      d->allocateWorkVector(DynamicalSystem::qtmp, d->q()->size());
    }
    else
    {
      RuntimeException::selfThrow("MoreauJeanCombinedProjectionOSI::initialize() - DS not of the right type");
    }
  }
}





bool MoreauJeanCombinedProjectionOSI::addInteractionInIndexSet(SP::Interaction inter, unsigned int i)
{
  assert(i == 1 || i == 2);
  //double h = _simulation->timeStep();
  if (i == 1) // index set for resolution at the velocity
  {
    double y = (inter->y(0))->getValue(0); // y(0) is the position
    DEBUG_PRINTF("MoreauJeanCombinedProjectionOSI::addInteractionInIndexSet yref=%e \n", y);
#ifdef DEBUG_MESSAGES
    if (y <= 0)
      DEBUG_PRINTF("MoreauJeanCombinedProjectionOSI::addInteractionInIndexSet ACTIVATE in indexSet level = %i.\n", i);
#endif
    return (y <= 0);
  }
  else if (i == 2)  //  special index for the projection
  {
    double lambda = 0;
    lambda = (inter->lambda(1))->getValue(0); // lambda(1) is the contact impulse for MoreauJeanOSI scheme
    DEBUG_PRINTF("MoreauJeanCombinedProjectionOSI::addInteractionInIndexSet lambdaref=%e \n", lambda);
#ifdef DEBUG_MESSAGES
    if (lambda > 0)
      DEBUG_PRINTF("MoreauJeanCombinedProjectionOSI::addInteractionInIndexSet ACTIVATE in indexSet level = %i.\n", i);
#endif
    //    return (lambda > 0);
    return true;
  }
  else
  {
    return false;
  }
}


bool MoreauJeanCombinedProjectionOSI::removeInteractionInIndexSet(SP::Interaction inter, unsigned int i)
{
  assert(0);
  return(0);
}

