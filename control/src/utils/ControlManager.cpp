/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

#include "EventsManager.hpp"
#include "Simulation.hpp"
#include "Simulation.hpp"
#include "TimeDiscretisation.hpp"
#include "OneStepIntegrator.hpp"
#include "NonSmoothDynamicalSystem.hpp"

#include "ControlManager.hpp"
#include "Sensor.hpp"
#include "SensorFactory.hpp"
#include "ActuatorFactory.hpp"
#include "Observer.hpp"
#include "ObserverFactory.hpp"
#include "SensorEvent.hpp"
#include "ActuatorEvent.hpp"
#include "ObserverEvent.hpp"
#include "ExtraAdditionalTerms.hpp"
//#define DEBUG_BEGIN_END_ONLY
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

ControlManager::ControlManager(SP::Simulation sim): _sim(sim)
{
  if(!_sim)
    THROW_EXCEPTION("ControlManager::constructor failed. The given Simulation is a NULL pointer.");
}

ControlManager::~ControlManager()
{}

void ControlManager::initialize(const NonSmoothDynamicalSystem& nsds)
{
  DEBUG_BEGIN("ControlManager::initialize(const NonSmoothDynamicalSystem& nsds)\n")
  // Initialize all the Sensors and insert their events into the
  // EventsManager of the Simulation.
  for(SensorsIterator itS = _allSensors.begin();
      itS != _allSensors.end(); ++itS)
  {
    (*itS)->initialize(nsds);
  }
  // Initialize all the Actuators and insert their events into the
  // EventsManager of the Simulation.
  for(ActuatorsIterator itA = _allActuators.begin();
      itA != _allActuators.end(); ++itA)
  {
    (*itA)->initialize(nsds,*_sim);
  }

  // Initialize all the Observer and insert their events into the
  // EventsManager of the Simulation.
  for(ObserversIterator itO = _allObservers.begin();
      itO != _allObservers.end(); ++itO)
  {
    (*itO)->initialize(nsds,*_sim);
  }

  // init the control terms, if any
  //OSISet& allOSI = *m.simulation()->oneStepIntegrators();
  OSISet& allOSI = *_sim->oneStepIntegrators();
  DynamicalSystemsGraph& DSG0 = *nsds.topology()->dSG(0);
  for(OSIIterator itosi = allOSI.begin(); itosi != allOSI.end(); ++itosi)
  {
    if((*itosi)->extraAdditionalTerms())
    {
      // would be nice to check is those are for Control
      (*itosi)->extraAdditionalTerms()->init(DSG0, nsds, _sim->eventsManager()->timeDiscretisation());
    }
  }
  DEBUG_END("ControlManager::initialize(const NonSmoothDynamicalSystem& nsds)\n");
}

SP::Sensor ControlManager::addSensor(int type, SP::TimeDiscretisation td, SP::DynamicalSystem ds)
{
  SensorFactory::Registry& regSensor(SensorFactory::Registry::get()) ;
  SP::Sensor s = (* (_allSensors.insert(regSensor.instantiate(type, ds))).first);
  linkSensorSimulation(s, td);
  return s;
}

SP::Sensor ControlManager::addAndRecordSensor(int type, SP::TimeDiscretisation td, SP::DynamicalSystem ds, const NonSmoothDynamicalSystem& nsds)
{
  SensorFactory::Registry& regSensor(SensorFactory::Registry::get()) ;
  SP::Sensor s = *(_allSensors.insert(regSensor.instantiate(type, ds))).first;
  linkSensorSimulation(s, td);
  s->initialize(nsds);
  return s;
}

SP::Actuator ControlManager::addActuator(int type, SP::TimeDiscretisation td, SP::ControlSensor sensor)
{
  if(!sensor)
    THROW_EXCEPTION("ControlManager::addActuator - sensor is not valid !");
  ActuatorFactory::Registry& regActuator(ActuatorFactory::Registry::get()) ;
  SP::Actuator act = (* (_allActuators.insert(regActuator.instantiate(type, sensor))).first);
  linkActuatorSimulation(act, td);
  return act;
}

SP::Actuator ControlManager::addAndRecordActuator(int type, SP::TimeDiscretisation td, SP::ControlSensor sensor, const NonSmoothDynamicalSystem& nsds)
{
  if(!sensor)
    THROW_EXCEPTION("ControlManager::addActuator - sensor is not valid !");
  ActuatorFactory::Registry& regActuator(ActuatorFactory::Registry::get()) ;
  SP::Actuator act = *(_allActuators.insert(regActuator.instantiate(type, sensor))).first;
  linkActuatorSimulation(act, td);
  act->initialize(nsds,*_sim);
  return act;
}

SP::Observer ControlManager::addObserver(int type, SP::TimeDiscretisation td, SP::ControlSensor sensor, const SiconosVector& xHat0)
{
  if(!sensor)
    THROW_EXCEPTION("ControlManager::addActuator - sensor is not valid !");
  ObserverFactory::Registry& regObserver(ObserverFactory::Registry::get()) ;
  SP::Observer obs = (* (_allObservers.insert(regObserver.instantiate(type, sensor, xHat0))).first);
  linkObserverSimulation(obs, td);
  return obs;
}

SP::Observer ControlManager::addAndRecordObserver(int type, SP::TimeDiscretisation td, SP::ControlSensor sensor, const SiconosVector& xHat0, const NonSmoothDynamicalSystem& nsds)
{
  ObserverFactory::Registry& regObserver(ObserverFactory::Registry::get()) ;
  SP::Observer obs = *(_allObservers.insert(regObserver.instantiate(type, sensor, xHat0))).first;
  linkObserverSimulation(obs, td);
  obs->initialize(nsds,*_sim);
  return obs;
}

void ControlManager::addSensorPtr(SP::Sensor s, SP::TimeDiscretisation td)
{
  if(!s)
    THROW_EXCEPTION("ControlManager::addActuator - sensor is not valid !");
  _allSensors.insert(s);
  linkSensorSimulation(s, td);
}

void ControlManager::addAndRecordSensorPtr(SP::Sensor s, SP::TimeDiscretisation td, const NonSmoothDynamicalSystem& nsds)
{
  if(!s)
    THROW_EXCEPTION("ControlManager::addActuator - sensor is not valid !");
  _allSensors.insert(s);
  linkSensorSimulation(s, td);
  s->initialize(nsds);
}

void ControlManager::addActuatorPtr(SP::Actuator act, SP::TimeDiscretisation td)
{
  if(!act)
    THROW_EXCEPTION("ControlManager::addActuator - actuator is not valid !");
  _allActuators.insert(act);
  linkActuatorSimulation(act, td);
}

void ControlManager::addAndRecordActuatorPtr(SP::Actuator act, SP::TimeDiscretisation td, const NonSmoothDynamicalSystem& nsds)
{
  if(!act)
    THROW_EXCEPTION("ControlManager::addActuator - actuator is not valid !");
  _allActuators.insert(act);
  linkActuatorSimulation(act, td);
  act->initialize(nsds,*_sim);
}

void ControlManager::addObserverPtr(SP::Observer obs, SP::TimeDiscretisation td)
{
  if(!obs)
    THROW_EXCEPTION("ControlManager::addActuator - observer is not valid !");
  _allObservers.insert(obs);
  linkObserverSimulation(obs, td);
}

void ControlManager::addAndRecordObserverPtr(SP::Observer obs, SP::TimeDiscretisation td, const NonSmoothDynamicalSystem& nsds)
{
  if(!obs)
    THROW_EXCEPTION("ControlManager::addActuator - observer is not valid !");
  _allObservers.insert(obs);
  linkObserverSimulation(obs, td);
  obs->initialize(nsds,*_sim);
}

void ControlManager::linkSensorSimulation(SP::Sensor s, SP::TimeDiscretisation td)
{
  Event& ev = _sim->eventsManager()->insertEvent(SENSOR_EVENT, td);
  static_cast<SensorEvent&>(ev).setSensorPtr(s);
  s->setTimeDiscretisation(*td);
}

void ControlManager::linkActuatorSimulation(SP::Actuator act, SP::TimeDiscretisation td)
{
  Event& ev = _sim->eventsManager()->insertEvent(ACTUATOR_EVENT, td);
  static_cast<ActuatorEvent&>(ev).setActuatorPtr(act);
  act->setTimeDiscretisation(*td);
}

void ControlManager::linkObserverSimulation(SP::Observer obs, SP::TimeDiscretisation td)
{
  Event& ev = _sim->eventsManager()->insertEvent(OBSERVER_EVENT, td);
  static_cast<ObserverEvent&>(ev).setObserverPtr(obs);
  obs->setTimeDiscretisation(*td);
}

void ControlManager::display() const
{
  std::cout << "=========> ControlManager " ;
  std::cout << "It handles the following objects: " <<std::endl;
  SensorsIterator itS;
  for(itS = _allSensors.begin(); itS != _allSensors.end(); ++itS)
    (*itS)->display();
  ActuatorsIterator itA;
  for(itA = _allActuators.begin(); itA != _allActuators.end(); ++itA)
    (*itA)->display();
  ObserversIterator itO;
  for(itO = _allObservers.begin(); itO != _allObservers.end(); ++itO)
    (*itO)->display();
  std::cout << "==========" << std::endl;
  std::cout << std::endl;
}
