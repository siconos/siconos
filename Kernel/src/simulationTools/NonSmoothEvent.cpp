/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
#include "NonSmoothEvent.h"
#include "EventFactory.h"
#include "EventDriven.h"
#include "OneStepIntegrator.h"

using namespace std;
using namespace EventFactory;

// Default constructor
NonSmoothEvent::NonSmoothEvent(): Event(0.0, "NonSmoothEvent")
{}

NonSmoothEvent::NonSmoothEvent(double time, const std::string& name): Event(time, "NonSmoothEvent")
{}

NonSmoothEvent::~NonSmoothEvent()
{}

void NonSmoothEvent::process(Simulation* simulation)
{
  if (simulation->getType() != "EventDriven")
    RuntimeException::selfThrow("NonSmoothEvent::process failed; Simulation is not of EventDriven type.");

  if (!(simulation->getOneStepNSProblems().empty()))
  {
    EventDriven * eventDriven = static_cast<EventDriven*>(simulation);

    // Compute y[0], y[1] and update index sets. => already done during advance to event ...

    //       simulation->updateOutput(0, 1);

    //       simulation->updateIndexSets();

    // Get the required index sets ...
    UnitaryRelationsSet * indexSet1 = simulation->getIndexSetPtr(1);
    UnitaryRelationsSet * indexSet2 = simulation->getIndexSetPtr(2);

    // ---> solve impact LCP if IndexSet[1]\IndexSet[2] is not empty.
    if (!(*indexSet1 - *indexSet2).isEmpty())
    {

      // For Event-Driven algo., memories vectors are of size 2 (ie 2 blocks).
      // First block (pos 0, last in) for post-event values and last block (pos 1, first in) for
      // pre-event values.

      simulation->saveInMemory();  // To save pre-impact values

      // solve the LCP-impact => y[1],lambda[1]
      eventDriven->computeOneStepNSProblem("impact"); // solveLCPImpact();

      // compute p[1], post-impact velocity, y[1] and indexSet[2]
      eventDriven->update(1);
      // Update the corresponding index set ...
      eventDriven->updateIndexSets();

      // check that IndexSet[1]-IndexSet[2] is now empty
      //    if( !((*indexSet1-*indexSet2).isEmpty()))
      //      RuntimeException::selfThrow("NonSmoothEvent::process, error after impact-LCP solving.");
    }

    // ---> solve acceleration LCP if IndexSet[2] is not empty
    if (!((indexSet2)->isEmpty()))
    {
      // Update the state of the DS
      //    OSIIterator itOSI;
      //    for(itOSI = simulation->getOneStepIntegrators().begin(); itOSI!=simulation->getOneStepIntegrators().end() ; ++itOSI)
      //      (*itOSI)->updateState(2);

      // solve LCP-acceleration
      eventDriven->computeOneStepNSProblem("acceleration"); //solveLCPAcceleration();

      // for all index in IndexSets[2], update the index set according to y[2] and/or lambda[2] sign.
      eventDriven->updateIndexSetsWithDoubleCondition();
    }

    // Save results in memory
    simulation->saveInMemory();
  }
}

AUTO_REGISTER_EVENT("NonSmoothEvent", NonSmoothEvent);
