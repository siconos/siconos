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
#include "TimeDiscretisationEvent.hpp"
#include "EventFactory.hpp"
#include "Simulation.hpp"
using namespace std;
using namespace EventFactory;

// Default constructor
TimeDiscretisationEvent::TimeDiscretisationEvent(): Event(0.0, TD_EVENT)
{}

TimeDiscretisationEvent::TimeDiscretisationEvent(double time, int notUsed): Event(time, TD_EVENT)
{}

TimeDiscretisationEvent::~TimeDiscretisationEvent()
{}

void TimeDiscretisationEvent::process(SP::Simulation simulation)
{
  // Update y[i] values in Interactions with new DS states.
  //simulation->updateOutput(0, 1);
  // Save state(s) in Memories (DS and Interactions, through OSI and OSNS).
  simulation->saveInMemory();
}

AUTO_REGISTER_EVENT(TD_EVENT, TimeDiscretisationEvent);
