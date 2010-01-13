/* Siconos-Kernel, Copyright INRIA 2005-2010.
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

#include "Sensor.hpp"
#include "SensorEvent.hpp"
#include "Model.hpp"
#include "TimeDiscretisation.hpp"
#include "EventFactory.hpp"
#include "Simulation.hpp"
#include <iostream>
using namespace std;

Sensor::Sensor(): _type(0), _id("none")
{}

Sensor::Sensor(int name, SP::TimeDiscretisation t): _type(name), _id("none"), _timeDiscretisation(t)
{}

Sensor::~Sensor()
{}

void Sensor::initialize()
{
  // == Create an event linked to the present Actuator. ==
  // Uses the events factory to insert the new event.
  EventFactory::Registry& regEvent(EventFactory::Registry::get());
  _eSensor = regEvent.instantiate(_timeDiscretisation->currentTime(), 4);
  boost::static_pointer_cast<SensorEvent>(_eSensor)->setSensorPtr(shared_from_this());
}

// Add the present sensor into the Simulation process
// i.e. add eSensor into the EventsManager of the simulation
void Sensor::recordInSimulation()
{
  model()->simulation()->eventsManager()->insertEvent(_eSensor);
}

void Sensor::display() const
{
  cout << "=====> Sensor of type " << _type << ", named " << _id ;
  if (_model)
    cout << " and linked to model named " << model()->title() << "." << endl;
  else
    cout << " and not linked to a model." << endl;
  cout << "======" << endl ;
  cout << endl;
}
