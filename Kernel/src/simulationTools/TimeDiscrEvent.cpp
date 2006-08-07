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
#include "TimeDiscrEvent.h"
#include "TimeDiscretisation.h"
using namespace std;

// Default constructor
TimeDiscrEvent::TimeDiscrEvent(): Event(DEFAULT_EVENT_TIME, "TimeDiscretisationEvent")
{}

// copy constructor
TimeDiscrEvent::TimeDiscrEvent(const TimeDiscrEvent& newTimeDiscrEvent): Event(newTimeDiscrEvent)
{}

TimeDiscrEvent::TimeDiscrEvent(const unsigned long int& time): Event(time, "TimeDiscretisationEvent")
{}

TimeDiscrEvent::~TimeDiscrEvent()
{}

void TimeDiscrEvent::process(Simulation* simulation)
{
  OSNSIterator itOsns;
  OneStepNSProblems allOSNS = simulation->getOneStepNSProblems();
  //  for(itOsns=(simulation->getOneStepNSProblems()).begin();itOsns!=(simulation->getOneStepNSProblems()).end();++itOsns)
  for (itOsns = allOSNS.begin(); itOsns != allOSNS.end(); ++itOsns)
    (itOsns->second)->updateOutput(0, 1);

  TimeDiscretisation * td =  simulation->getTimeDiscretisationPtr();
  td->increment();
  Model * model =  simulation->getModelPtr();
  model->setCurrentT(model->getCurrentT() + td->getH());
  simulation->nextStep();

  //todo
  // cout << " Time Discr Event processing: nothing implemented" << endl;
}
