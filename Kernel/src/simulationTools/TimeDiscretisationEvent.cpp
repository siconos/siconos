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
#include "TimeDiscretisationEvent.h"
#include "EventFactory.h"
#include "Simulation.h"
using namespace std;
using namespace EventFactory;

// Default constructor
TimeDiscretisationEvent::TimeDiscretisationEvent(): Event(0.0, "TimeDiscretisationEvent")
{}

TimeDiscretisationEvent::TimeDiscretisationEvent(double time, const std::string& name): Event(time, "TimeDiscretisationEvent")
{}

TimeDiscretisationEvent::~TimeDiscretisationEvent()
{}

void TimeDiscretisationEvent::process(Simulation* simulation)
{
  // Update y[i] values in Interactions with new DS states.
  //simulation->updateOutput(0, 1);
  // Save state(s) in Memories (DS and Interactions, through OSI and OSNS).
  simulation->saveInMemory();
}

AUTO_REGISTER_EVENT("TimeDiscretisationEvent", TimeDiscretisationEvent);
