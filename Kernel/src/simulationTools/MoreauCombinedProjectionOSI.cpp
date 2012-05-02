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
#include "MoreauCombinedProjectionOSI.hpp"
#include "Simulation.hpp"
// #include "Model.hpp"
// #include "NonSmoothDynamicalSystem.hpp"
// #include "NewtonEulerDS.hpp"
// #include "LagrangianLinearTIDS.hpp"
// #include "FirstOrderLinearTIDS.hpp"
// #include "NewtonEulerR.hpp"
// #include "LagrangianRheonomousR.hpp"
// #include "FirstOrderLinearTIR.hpp"
// #include "FirstOrderLinearR.hpp"
// #include "NewtonImpactNSL.hpp"
// #include "MultipleImpactNSL.hpp"
// #include "NewtonImpactFrictionNSL.hpp"


#define DEBUG_MESSAGES
#define DEBUG_WHERE_MESSAGES
#include <debug.h>

using namespace std;
//using namespace RELATION;

bool MoreauCombinedProjectionOSI::addInteractionInIndexSet(SP::Interaction inter, unsigned int i)
{
  assert(i == 1 || i == 2);
  //double h = simulationLink->timeStep();
  if (i == 1) // index set for resolution at the velocity
  {
    double y = (inter->y(0))->getValue(0); // y(0) is the position
    DEBUG_PRINTF("MoreauCombinedProjectionOSI::addInteractionInIndexSet yref=%e \n", y);
    if (y <= 0)
      DEBUG_PRINTF("MoreauCombinedProjectionOSI::addInteractionInIndexSet ACTIVATE in indeSet level = %i.\n", i);
    return (y <= 0);
  }
  else if (i == 2)  //  special index for the projection
  {
    double lambda = (inter->lambda(1))->getValue(0); // lambda(1) is the contact impulse for Moreau scheme
    DEBUG_PRINTF("MoreauCombinedProjectionOSI::addInteractionInIndexSet lambdaref=%e \n", lambda);
    if (lambda > 0)
      DEBUG_PRINTF("MoreauCombinedProjectionOSI::addInteractionInIndexSet ACTIVATE in indeSet level = %i.\n", i);
    return (lambda > 0);
  }
  return(0);
}


bool MoreauCombinedProjectionOSI::removeInteractionInIndexSet(SP::Interaction inter, unsigned int i)
{
  assert(0);
  return(0);
}

