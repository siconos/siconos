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
#include "Observer.hpp"
#include "ObserverFactory.hpp"
#include "Simulation.hpp"
#include "TimeDiscretisation.hpp"



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
  }
  // Initialize all the Actuators and insert their events into the
  // EventsManager of the Simulation.
  for (ActuatorsIterator itA = _allActuators.begin();
       itA != _allActuators.end(); ++itA)
  {
    (*itA)->initialize(_model);
  }

  // Initialize all the Actuators and insert their events into the
  // EventsManager of the Simulation.
  for (ObserversIterator itO = _allObservers.begin();
       itO != _allObservers.end(); ++itO)
  {
    (*itO)->initialize(_model);
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
  return tmp;
}

SP::Observer ControlManager::addObserver(int type, SP::TimeDiscretisation t, unsigned int number)
{
  ObserverFactory::Registry& regObserver(ObserverFactory::Registry::get()) ;
  return (* (_allObservers.insert(regObserver.instantiate(type, t, getDSFromModel(number)))).first);
}

SP::Observer ControlManager::addAndRecordObserver(int type, SP::TimeDiscretisation t, unsigned int number)
{
  ObserverFactory::Registry& regObserver(ObserverFactory::Registry::get()) ;
  SP::Observer tmp = *(_allObservers.insert(regObserver.instantiate(type, t, getDSFromModel(number)))).first;
  tmp->initialize(_model);
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
}

void ControlManager::addActuatorPtr(SP::Actuator act)
{
  _allActuators.insert(act);
}

void ControlManager::addAndRecordActuatorPtr(SP::Actuator act)
{
  _allActuators.insert(act);
  act->initialize(_model);
}

void ControlManager::addObserverPtr(SP::Observer obs)
{
  _allObservers.insert(obs);
}

void ControlManager::addAndRecordObserverPtr(SP::Observer obs)
{
  _allObservers.insert(obs);
  obs->initialize(_model);
}

void ControlManager::display() const
{
  std::cout << "=========> ControlManager " ;
  if (model())
    std::cout << "linked to model named: " << model()->title() << "." <<std::endl;
  else
    std::cout << "not linked to a model." <<std::endl;
  std::cout << "It handles the following objects: " <<std::endl;
  SensorsIterator itS;
  for (itS = _allSensors.begin(); itS != _allSensors.end(); ++itS)
    (*itS)->display();
  ActuatorsIterator itA;
  for (itA = _allActuators.begin(); itA != _allActuators.end(); ++itA)
    (*itA)->display();
  ObserversIterator itO;
  for (itO = _allObservers.begin(); itO != _allObservers.end(); ++itO)
    (*itO)->display();
  std::cout << "==========" << std::endl;
  std::cout << std::endl;
}
