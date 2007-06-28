/* Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
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

#include "ControlManager.h"
#include "EventsManager.h"
#include "Model.h"
#include "Sensor.h"
#include "SensorFactory.h"
#include "ActuatorFactory.h"
#include "Simulation.h"
#include "TimeDiscretisation.h"

using namespace std;

ControlManager::ControlManager(): model(NULL)
{}

ControlManager::ControlManager(Model* m): model(m)
{
  if (model == NULL)
    RuntimeException::selfThrow("ControlManager::constructor failed. The given Model is a NULL pointer.");
}

ControlManager::~ControlManager()
{
  // Delete Sensors
  SensorsIterator itS;
  for (itS = allSensors.begin(); itS != allSensors.end(); ++itS)
    if ((*itS) != NULL) delete(*itS);
  allSensors.clear();
  // Delete Actuators
  ActuatorsIterator itA;
  for (itA = allActuators.begin(); itA != allActuators.end(); ++itA)
    if ((*itA) != NULL) delete(*itA);
  allActuators.clear();
}

void ControlManager::initialize()
{
  EventsManager * eventsManager = model->getSimulationPtr()->getEventsManagerPtr();
  // Initialize all the Sensors and insert their events into the EventsManager of the Simulation.
  SensorsIterator itS;
  for (itS = allSensors.begin(); itS != allSensors.end(); ++itS)
  {
    (*itS)->initialize();
    eventsManager->insertEvents((*itS)->getEvents());
  }
  // Initialize all the Actuators and insert their events into the EventsManager of the Simulation.
  ActuatorsIterator itA;
  for (itA = allActuators.begin(); itA != allActuators.end(); ++itA)
  {
    (*itA)->initialize();
    eventsManager->insertEvents((*itA)->getEvents());
  }
}

Sensor* ControlManager::addSensor(const string& type, TimeDiscretisation* t)
{
  if (t->getModelPtr() != model)
    RuntimeException::selfThrow("ControlManager::addSensor(...,timeDiscretisation) failed. The Model linked to the controlManager is different from the one of the TimeDiscretisation.");

  SensorFactory::Registry& regSensor(SensorFactory::Registry::get()) ;
  return (* (allSensors.insert(regSensor.instantiate(type, t))).first);
}

Actuator* ControlManager::addActuator(const string& type, TimeDiscretisation* t)
{
  if (t->getModelPtr() != model)
    RuntimeException::selfThrow("ControlManager::addActuator(...,timeDiscretisation) failed. The Model linked to the controlManager is different from the one of the TimeDiscretisation.");

  ActuatorFactory::Registry& regActuator(ActuatorFactory::Registry::get()) ;
  return (* (allActuators.insert(regActuator.instantiate(type, t))).first);
}

void ControlManager::display() const
{
  cout << "=========> ControlManager " ;
  if (model != NULL)
    cout << "linked to model named: " << model->getTitle() << "." << endl;
  else
    cout << "not linked to a model." << endl;
  cout << "It handles the following objects: " << endl;
  SensorsIterator itS;
  for (itS = allSensors.begin(); itS != allSensors.end(); ++itS)
    (*itS)->display();
  // Initialize all the Actuators.
  ActuatorsIterator itA;
  for (itA = allActuators.begin(); itA != allActuators.end(); ++itA)
    (*itA)->display();
  cout << "==========" << endl;
  cout << endl;
}
