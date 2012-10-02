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

#include "Sensor.hpp"
#include "SensorEvent.hpp"
#include "DynamicalSystem.hpp"
#include "TimeDiscretisation.hpp"
#include "EventFactory.hpp"
#include "Simulation.hpp"
#include <iostream>
using namespace std;

Sensor::Sensor(): _type(0), _id("none")
{}

Sensor::Sensor(int name, SP::TimeDiscretisation t, SP::DynamicalSystem ds): _type(name), _id("none"), _DS(ds), _timeDiscretisation(t)
{
  _nDim = _DS->getN();
}

Sensor::~Sensor()
{}

void Sensor::initialize(SP::Model m)
{
  _model = m;
  // == Create an event linked to the present Sensor. ==
  // Uses the events factory to insert the new event.
  EventFactory::Registry& regEvent(EventFactory::Registry::get());
  _eSensor = regEvent.instantiate(_timeDiscretisation->currentTime(), SENSOR_EVENT);
  std11::static_pointer_cast<SensorEvent>(_eSensor)->setSensorPtr(shared_from_this());
}

// Add the present sensor into the Simulation process
// i.e. add eSensor into the EventsManager of the simulation
void Sensor::recordInSimulation()
{
  _model->simulation()->eventsManager()->insertEvent(_eSensor);
}

void Sensor::display() const
{
  cout << "=====> Sensor of type " << _type << ", named " << _id ;
  if (_DS)
    cout << " and linked to the DynamicalSystem number " << _DS->number() << "." << endl;
  else
    cout << " and not linked to a DynamicalSystem." << endl;
  cout << "======" << endl ;
  cout << endl;
}
