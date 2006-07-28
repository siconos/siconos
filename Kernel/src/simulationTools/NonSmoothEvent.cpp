/* Siconos-Kernel version 1.2.0, Copyright INRIA 2005-2006.
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
  //  cout << "Non Smooth Event ... " << endl;
  if (!(simulation->getOneStepNSProblems().empty()))
  {
    EventDriven * eventDriven = static_cast<EventDriven*>(simulation);

    // Compute y[0], y[1] and update index sets.

    OSNSIterator itOsns;
    OneStepNSProblems allOSNS = simulation->getOneStepNSProblems();
    for (itOsns = allOSNS.begin(); itOsns != allOSNS.end(); ++itOsns)
      //      for(itOsns=simulation->getOneStepNSProblems().begin();itOsns!=simulation->getOneStepNSProblems().end();++itOsns)
      (itOsns->second)->updateOutput(0, 1);

    simulation->updateIndexSets();

    VectorOfSetOfUnitaryRelations indexSets = eventDriven->getIndexSets();
    cout << " NS EVENT PROCESS " << indexSets[1].size() << " " << indexSets[2].size() << endl;
    // ---> solve impact LCP if IndexSet[1]\IndexSet[2] is not empty.
    if (!(indexSets[1] - indexSets[2]).isEmpty())
    {
      cout << "SOLVE LCP IMPACT " << endl;
      // solve the LCP-impact => y[1],lambda[1]
      eventDriven->computeOneStepNSProblem("impact"); // solveLCPImpact();

      // compute p[1] and post-impact velocity
      eventDriven->updateImpactState();

      //  update indexSet that depends on y[1]
      eventDriven->updateIndexSet(2);

      // check that IndexSet[1]-IndexSet[2] is now empty
      if (!(indexSets[1] - indexSets[2]).isEmpty())
        RuntimeException::selfThrow("EventDriven advanceToEvent, error after impact-LCP solving.");
    }

    if (!((indexSets[2]).isEmpty()))
    {
      cout << "SOLVE LCP ACCELERATION " << endl;

      // Update the state of the DS
      OSIIterator itOSI;
      double time = simulation->getModelPtr()->getCurrentT();
      for (itOSI = simulation->getOneStepIntegrators().begin(); itOSI != simulation->getOneStepIntegrators().end() ; ++itOSI)
        (*itOSI)->updateState(time, 2);

      // solve LCP-acceleration
      eventDriven->computeOneStepNSProblem("acceleration"); //solveLCPAcceleration();

      // for all index in IndexSets[2], update the index set according to y[2] and/or lambda[2] sign.
      eventDriven->updateIndexSetsWithDoubleCondition();
    }

    simulation->nextStep();

  }
  cout << "End of Non Smooth Event ... " << endl;

}
