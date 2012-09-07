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

#include "Actuator.hpp"
#include "ActuatorEvent.hpp"
#include "Sensor.hpp"
#include "Model.hpp"
#include "TimeDiscretisation.hpp"
#include "EventFactory.hpp"
#include "DynamicalSystem.hpp"
#include "Simulation.hpp"
#include <iostream>
using namespace std;

Actuator::Actuator(): _type(0), _id("none")
{
  _allDS.reset(new DynamicalSystemsSet());
  _allSensors.reset(new Sensors());
}

Actuator::Actuator(int name, SP::TimeDiscretisation t, SP::DynamicalSystem ds): _type(name), _id("none"), _DS(ds), _timeDiscretisation(t)
{
  _allDS.reset(new DynamicalSystemsSet());
  _allSensors.reset(new Sensors());
  addDynamicalSystemPtr(ds);
  _nDim = ds->getN();
}

Actuator::Actuator(int name, SP::TimeDiscretisation t, SP::DynamicalSystem ds, const Sensors& sensorList): _type(name), _id("none"), _DS(ds), _timeDiscretisation(t)
{
  _allDS.reset(new DynamicalSystemsSet());
  _allSensors.reset(new Sensors());
  addDynamicalSystemPtr(ds);
  _nDim = ds->getN();
  addSensors(sensorList);
}

Actuator::~Actuator()
{
  _allDS->clear();
  _allSensors->clear();
}

void Actuator::addSensors(const Sensors& newSensors)
{
  // Add all the Sensor of newSensors into allSensors.
  // => allSensors is not cleared and so all existing Sensors remain.
  // => no copy of Sensors but copy of the pointers
  for (SensorsIterator itS = newSensors.begin(); itS != newSensors.end(); ++itS)
    _allSensors->insert(*itS);

}

void Actuator::addSensorPtr(SP::Sensor newSensor)
{
  // Add a Sensor into allSensors set: no copy, pointer link.
  _allSensors->insert(newSensor);
}

void Actuator::addDynamicalSystems(const DynamicalSystemsSet& newDSs)
{
  // Add all the DS of newDSs into allDS.
  // => allDS is not cleared and so all existing DSs remain.
  // => no copy of DS but copy of the pointers
  for (ConstDSIterator itDS = newDSs.begin(); itDS != newDSs.end(); ++itDS)
    _allDS->insert(*itDS);
}

void Actuator::addDynamicalSystemPtr(SP::DynamicalSystem newDS)
{
  // Add a DS into allDS set: no copy, pointer link.
  _allDS->insert(newDS);
}

void Actuator::initialize(SP::Model m)
{
  _model = m;
  // == Create an event linked to the present Actuator. ==
  // Uses the events factory to insert the new event.
  EventFactory::Registry& regEvent(EventFactory::Registry::get());
  _eActuator = regEvent.instantiate(_timeDiscretisation->currentTime(), ACTUATOR_EVENT);
  cpp11ns::static_pointer_cast<ActuatorEvent>(_eActuator)->setActuatorPtr(shared_from_this());

  // Warning: no Sensors initialization. They are supposed to be up to date when added in the Actuator.
}

// Add the present actuator into the Simulation process
// i.e. add eActuator into the EventsManager of the simulation
void Actuator::recordInSimulation()
{
  _model->simulation()->eventsManager()->insertEvent(_eActuator);
}

void Actuator::display() const
{
  cout << "=====> Actuator of type " << _type << ", named " << _id ;
  cout << "The associated Sensors are: " << endl;
  for (SensorsIterator itS = _allSensors->begin(); itS != _allSensors->end(); ++itS)
    (*itS)->display();

  cout << "The associated DynamicalSystems are: " << endl;
  for (DSIterator itDS = _allDS->begin(); itDS != _allDS->end(); ++itDS)
    cout << " - Numbers : " << (*itDS)->number()  << endl;
  cout << "======" << endl;
  cout << endl;
}
