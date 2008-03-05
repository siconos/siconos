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

#include "Actuator.h"
#include "ActuatorEvent.h"
#include "Sensor.h"
#include "Model.h"
#include "TimeDiscretisation.h"
#include "EventFactory.h"
#include "DynamicalSystem.h"
#include <iostream>
using namespace std;

Actuator::Actuator(): type("generic"), id("none"), allSensors(NULL), allDS(NULL),  model(NULL), timeDiscretisation(NULL), eventsSet()
{
  allDS = new DynamicalSystemsSet();
  allSensors = new Sensors();
}

Actuator::Actuator(const std::string& name, TimeDiscretisation* t): type(name), id("none"), allSensors(NULL), allDS(NULL), model(t->getModelPtr()), timeDiscretisation(t), eventsSet()
{
  allDS = new DynamicalSystemsSet();
  allSensors = new Sensors();
}

Actuator::Actuator(const std::string& name, TimeDiscretisation* t, const Sensors& sensorList): type(name), id("none"), allSensors(NULL), allDS(NULL), model(t->getModelPtr()), timeDiscretisation(t), eventsSet()
{
  allDS = new DynamicalSystemsSet();
  allSensors = new Sensors();
}

Actuator::~Actuator()
{
  allDS->clear();
  allSensors->clear();
  delete allDS;
  delete allSensors;
}

void Actuator::addSensors(const Sensors& newSensors)
{
  // Add all the Sensor of newSensors into allSensors.
  // => allSensors is not cleared and so all existing Sensors remain.
  // => no copy of Sensors but copy of the pointers

  SensorsIterator itS;
  for (itS = newSensors.begin(); itS != newSensors.end(); ++itS)
    allSensors->insert(*itS);

}

void Actuator::addSensorPtr(Sensor * newSensor)
{
  // Add a Sensor into allSensors set: no copy, pointer link.
  allSensors->insert(newSensor);
}

void Actuator::addDynamicalSystems(const DynamicalSystemsSet& newDSs)
{
  // Add all the DS of newDSs into allDS.
  // => allDS is not cleared and so all existing DSs remain.
  // => no copy of DS but copy of the pointers

  DSIterator itDS;
  for (itDS = newDSs.begin(); itDS != newDSs.end(); ++itDS)
    allDS->insert(*itDS);
}

void Actuator::addDynamicalSystemPtr(DynamicalSystem * newDS)
{
  // Add a DS into allDS set: no copy, pointer link.
  allDS->insert(newDS);
}

void Actuator::initialize()
{
  // == Init. time discretisation data. ==
  timeDiscretisation->initialize();
  SiconosVector * tk = timeDiscretisation->getTkPtr();
  unsigned int sizeTk = tk->size();

  // == Create Events linked to the present Actuator. ==

  EventsContainerIterator it; // to check if insertion succeed or not.
  string type = "ActuatorEvent";

  // Uses the events factory to insert the new event.
  EventFactory::Registry& regEvent(EventFactory::Registry::get()) ;
  for (unsigned int i = 0; i < sizeTk; ++i)
    it = eventsSet.insert(regEvent.instantiate((*tk)(i), type));

  // == Set Actuator object of all Events to this ==
  for (it = eventsSet.begin(); it != eventsSet.end(); ++it)
    static_cast<ActuatorEvent*>((*it))->setActuatorPtr(this);

  // Warning: no Sensors initialization. They are supposed to be up to date when added in the Actuator.
}

void Actuator::display() const
{
  cout << "=====> Actuator of type " << type << ", named " << id ;
  if (model != NULL)
    cout << " and linked to model named " << model->getTitle() << "." << endl;
  else
    cout << " and not linked to a model." << endl;
  cout << "The associated Sensors are: " << endl;
  SensorsIterator itS;
  for (itS = allSensors->begin(); itS != allSensors->end(); ++itS)
    (*itS)->display();
  cout << "The associated DynamicalSystems are: " << endl;
  DSIterator itDS;
  for (itDS = allDS->begin(); itDS != allDS->end(); ++itDS)
    cout << " - Number and Id: " << (*itDS)->getNumber() << ", " << (*itDS)->getId() << endl;
  cout << "======" << endl;
  cout << endl;
}
