/* Siconos-Kernel version 2.0.0, Copyright INRIA 2005-2006.
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
#include "EventDriven.h"
using namespace std;

// Default constructor
NonSmoothEvent::NonSmoothEvent(): Event(DEFAULT_EVENT_TIME, "NonSmoothEvent")
{}

// copy constructor
NonSmoothEvent::NonSmoothEvent(const NonSmoothEvent& newNonSmoothEvent): Event(newNonSmoothEvent)
{}

NonSmoothEvent::NonSmoothEvent(const unsigned long int& time): Event(time, "NonSmoothEvent")
{}

// NonSmoothEvent(NonSmoothEventXML*, const std::string& ):timeOfNonSmoothEvent(0), type("undefined from xml")
//{}

NonSmoothEvent::~NonSmoothEvent()
{}

void NonSmoothEvent::process(Simulation* simulation)
{
  if (!(simulation->getOneStepNSProblems().empty()))
  {
    EventDriven * eventDriven = static_cast<EventDriven*>(simulation);

    // Compute y[0], y[1] and update index sets.
    simulation->updateOutput(0, 1);

    simulation->updateIndexSets();

    VectorOfSetOfUnitaryRelations indexSets = eventDriven->getIndexSets();

    // ---> solve impact LCP if IndexSet[1]\IndexSet[2] is not empty.
    if (!(indexSets[1] - indexSets[2]).isEmpty())
    {
      simulation->nextStep();  // To save pre-impact values
      // solve the LCP-impact => y[1],lambda[1]
      eventDriven->computeOneStepNSProblem("impact"); // solveLCPImpact();

      // compute p[1], post-impact velocity, y[1] and indexSet[2]
      eventDriven->update(1);

      //    for(itOsns=allOSNS.begin();itOsns!=allOSNS.end();++itOsns)
      //      (itOsns->second)->updateOutput(0);

      //    //  update indexSet that depends on y[0]
      //    eventDriven->updateIndexSet(1);

      // check that IndexSet[1]-IndexSet[2] is now empty

      indexSets = eventDriven->getIndexSets();

      //    if( !(indexSets[1]-indexSets[2]).isEmpty())
      //RuntimeException::selfThrow("NonSmoothEvent::process, error after impact-LCP solving.");
    }

    if (!((indexSets[2]).isEmpty()))
    {
      cout << "SOLVE LCP ACCELERATION " << endl;

      // Update the state of the DS
      OSIIterator itOSI;
      for (itOSI = simulation->getOneStepIntegrators().begin(); itOSI != simulation->getOneStepIntegrators().end() ; ++itOSI)
        (*itOSI)->updateState(2);

      // solve LCP-acceleration
      cout << "LCP acc solving ..." << endl;
      eventDriven->computeOneStepNSProblem("acceleration"); //solveLCPAcceleration();

      // for all index in IndexSets[2], update the index set according to y[2] and/or lambda[2] sign.
      eventDriven->updateIndexSetsWithDoubleCondition();
    }

    simulation->nextStep();
  }
}
