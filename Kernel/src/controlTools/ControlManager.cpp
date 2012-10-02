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

#include "ControlManager.hpp"
#include "EventsManager.hpp"
#include "Model.hpp"
#include "Sensor.hpp"
#include "SensorFactory.hpp"
#include "ActuatorFactory.hpp"
#include "Simulation.hpp"
#include "TimeDiscretisation.hpp"

using namespace std;

ControlManager::ControlManager(SP::Model m): _model(m)
{
  if (!_model)
    RuntimeException::selfThrow("ControlManager::constructor failed. The given Model is a NULL pointer.");
}

ControlManager::~ControlManager()
{}

void ControlManager::initialize()
{
  // Initialize all the Sensors and insert their events into the
  // EventsManager of the Simulation.
  for (SensorsIterator itS = _allSensors.begin();
       itS != _allSensors.end(); ++itS)
  {
    (*itS)->initialize(_model);
    (*itS)->recordInSimulation();
  }
  // Initialize all the Actuators and insert their events into the
  // EventsManager of the Simulation.
  for (ActuatorsIterator itA = _allActuators.begin();
       itA != _allActuators.end(); ++itA)
  {
    (*itA)->initialize(_model);
    (*itA)->recordInSimulation();
  }
}

SP::Sensor ControlManager::addSensor(int type, SP::TimeDiscretisation t, unsigned int number)
{
  SensorFactory::Registry& regSensor(SensorFactory::Registry::get()) ;
  return (* (_allSensors.insert(regSensor.instantiate(type, t, getDSFromModel(number)))).first);
}

SP::Sensor ControlManager::addAndRecordSensor(int type, SP::TimeDiscretisation t, unsigned int number)
{
  double currentTime = _model->simulation()->nextTime();
  while (t->currentTime() < currentTime)
    t->increment();
  SensorFactory::Registry& regSensor(SensorFactory::Registry::get()) ;
  SP::Sensor tmp = *(_allSensors.insert(regSensor.instantiate(type, t, getDSFromModel(number)))).first;
  tmp->initialize(_model);
  tmp->recordInSimulation();
  return tmp;
}

SP::Actuator ControlManager::addActuator(int type, SP::TimeDiscretisation t, unsigned int number)
{
  ActuatorFactory::Registry& regActuator(ActuatorFactory::Registry::get()) ;
  return (* (_allActuators.insert(regActuator.instantiate(type, t, getDSFromModel(number)))).first);
}

SP::Actuator ControlManager::addAndRecordActuator(int type, SP::TimeDiscretisation t, unsigned int number)
{
  ActuatorFactory::Registry& regActuator(ActuatorFactory::Registry::get()) ;
  SP::Actuator tmp = *(_allActuators.insert(regActuator.instantiate(type, t, getDSFromModel(number)))).first;
  tmp->initialize(_model);
  tmp->recordInSimulation();
  return tmp;
}

void ControlManager::addSensorPtr(SP::Sensor s)
{
  _allSensors.insert(s);
}

void ControlManager::addAndRecordSensorPtr(SP::Sensor s)
{
  _allSensors.insert(s);
  s->initialize(_model);
  s->recordInSimulation();
}

void ControlManager::addActuatorPtr(SP::Actuator act)
{
  _allActuators.insert(act);
}

void ControlManager::addAndRecordActuatorPtr(SP::Actuator act)
{
  _allActuators.insert(act);
  act->initialize(_model);
  act->recordInSimulation();
}

void ControlManager::display() const
{
  cout << "=========> ControlManager " ;
  if (model())
    cout << "linked to model named: " << model()->title() << "." << endl;
  else
    cout << "not linked to a model." << endl;
  cout << "It handles the following objects: " << endl;
  SensorsIterator itS;
  for (itS = _allSensors.begin(); itS != _allSensors.end(); ++itS)
    (*itS)->display();
  // Initialize all the Actuators.
  ActuatorsIterator itA;
  for (itA = _allActuators.begin(); itA != _allActuators.end(); ++itA)
    (*itA)->display();
  cout << "==========" << endl;
  cout << endl;
}
