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

#include "ControlManager.h"
#include "EventsManager.h"
#include "Model.h"
#include "Sensor.h"
#include "SensorFactory.h"
#include "ActuatorFactory.h"
#include "Simulation.h"
#include "TimeDiscretisation.h"

using namespace std;

ControlManager::ControlManager(SP::Model m): model(m)
{
  if (!model)
    RuntimeException::selfThrow("ControlManager::constructor failed. The given Model is a NULL pointer.");
}

ControlManager::~ControlManager()
{}

void ControlManager::initialize()
{
  // Initialize all the Sensors and insert their events into the EventsManager of the Simulation.
  for (SensorsIterator itS = allSensors.begin(); itS != allSensors.end(); ++itS)
  {
    (*itS)->initialize();
    (*itS)->recordInSimulation();
  }
  // Initialize all the Actuators and insert their events into the EventsManager of the Simulation.
  for (ActuatorsIterator itA = allActuators.begin(); itA != allActuators.end(); ++itA)
  {
    (*itA)->initialize();
    (*itA)->recordInSimulation();
  }
}

SP::Sensor ControlManager::addSensor(int type, SP::TimeDiscretisation t)
{
  SensorFactory::Registry& regSensor(SensorFactory::Registry::get()) ;
  return (* (allSensors.insert(regSensor.instantiate(type, t))).first);
}

SP::Sensor ControlManager::addAndRecordSensor(int type, SP::TimeDiscretisation t)
{
  double currentTime = model->getSimulationPtr()->getNextTime();
  while (t->getCurrentTime() < currentTime)
    t->increment();
  SensorFactory::Registry& regSensor(SensorFactory::Registry::get()) ;
  SP::Sensor tmp = *(allSensors.insert(regSensor.instantiate(type, t))).first;
  tmp->initialize();
  tmp->recordInSimulation();

  return tmp;
}

SP::Actuator ControlManager::addActuator(int type, SP::TimeDiscretisation t)
{
  ActuatorFactory::Registry& regActuator(ActuatorFactory::Registry::get()) ;
  return (* (allActuators.insert(regActuator.instantiate(type, t))).first);
}

SP::Actuator ControlManager::addAndRecordActuator(int type, SP::TimeDiscretisation t)
{
  ActuatorFactory::Registry& regActuator(ActuatorFactory::Registry::get()) ;
  SP::Actuator tmp = *(allActuators.insert(regActuator.instantiate(type, t))).first;
  tmp->initialize();
  tmp->recordInSimulation();
  return tmp;
}

void ControlManager::display() const
{
  cout << "=========> ControlManager " ;
  if (model)
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
