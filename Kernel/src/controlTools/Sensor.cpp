/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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

#include "Sensor.h"
#include "SensorEvent.h"
#include "Model.h"
#include "TimeDiscretisation.h"
#include "EventFactory.h"
#include "Simulation.h"
#include <iostream>
using namespace std;

Sensor::Sensor(): type(0), id("none"), model(NULL), timeDiscretisation(NULL), eSensor(NULL)
{}

Sensor::Sensor(int name, TimeDiscretisation* t): type(name), id("none"), model(t->getModelPtr()), timeDiscretisation(t), eSensor(NULL)
{}

Sensor::~Sensor()
{}

void Sensor::initialize()
{
  // == Create an event linked to the present Actuator. ==
  // Uses the events factory to insert the new event.
  EventFactory::Registry& regEvent(EventFactory::Registry::get());
  eSensor = regEvent.instantiate(timeDiscretisation->getCurrentTime(), 4);
  static_cast<SensorEvent*>(eSensor)->setSensorPtr(this);
}

// Add the present sensor into the Simulation process
// i.e. add eSensor into the EventsManager of the simulation
void Sensor::recordInSimulation()
{
  model->getSimulationPtr()->getEventsManagerPtr()->insertEvent(eSensor);
}

void Sensor::display() const
{
  cout << "=====> Sensor of type " << type << ", named " << id ;
  if (model != NULL)
    cout << " and linked to model named " << model->getTitle() << "." << endl;
  else
    cout << " and not linked to a model." << endl;
  cout << "======" << endl ;
  cout << endl;
}
