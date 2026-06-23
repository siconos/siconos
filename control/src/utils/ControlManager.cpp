/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

#include "ControlManager.hpp"

#include "Actuator.hpp"
#include "ActuatorEvent.hpp"
#include "EventsManager.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "Observer.hpp"
#include "ObserverEvent.hpp"
#include "OneStepIntegrator.hpp"
#include "Sensor.hpp"
#include "SensorEvent.hpp"
#include "SiconosException.hpp"
#include "Simulation.hpp"
#include "Topology.hpp"
// #define DEBUG_BEGIN_END_ONLY
//  #define DEBUG_NOCOLOR
//  #define DEBUG_STDOUT
//  #define DEBUG_MESSAGES
#include "siconos_debug.h"

siconos::control::ControlManager::ControlManager(
    std::shared_ptr<siconos::simulation::Simulation> sim)
    : _sim(sim) {
  if (!_sim)
    THROW_EXCEPTION(
        "siconos::control::ControlManager::constructor failed. The given Simulation is a NULL "
        "pointer.");
}

void siconos::control::ControlManager::initialize(
    const siconos::modeling::NonSmoothDynamicalSystem& nsds) {
  DEBUG_BEGIN(
      "siconos::control::ControlManager::initialize(const "
      "siconos::modeling::NonSmoothDynamicalSystem& nsds)\n")
  // Initialize all the Sensors and insert their events into the
  // EventsManager of the Simulation.
  for (auto itS : _allSensors) {
    itS->initialize(nsds);
  }
  // Initialize all the Actuators and insert their events into the
  // EventsManager of the Simulation.
  for (auto itA : _allActuators) {
    itA->initialize(nsds, *_sim);
  }

  // Initialize all the Observer and insert their events into the
  // EventsManager of the Simulation.
  for (auto itO : _allObservers) {
    itO->initialize(nsds, *_sim);
  }

  // init the control terms, if any
  // OSISet& allOSI = *m.simulation()->oneStepIntegrators();
  auto& allOSI = *_sim->oneStepIntegrators();
  auto& DSG0 = *nsds.topology()->dSG(0);
  for (auto itosi : allOSI) {
    if (itosi->extraAdditionalTerms()) {
      // would be nice to check is those are for Control
      itosi->extraAdditionalTerms()->init(DSG0, nsds,
                                          _sim->eventsManager()->timeDiscretisation());
    }
  }
  DEBUG_END(
      "siconos::control::ControlManager::initialize(const "
      "siconos::modeling::NonSmoothDynamicalSystem& nsds)\n");
}

std::shared_ptr<siconos::control::Sensor> siconos::control::ControlManager::addSensor(
    SensorType type, std::shared_ptr<siconos::simulation::TimeDiscretisation> td,
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds) {
  auto s = *(_allSensors.insert(SensorFactory::instance()->create(ds, type))).first;
  linkSensorSimulation(s, td);
  return s;
}

std::shared_ptr<siconos::control::Sensor> siconos::control::ControlManager::addAndRecordSensor(
    SensorType type, std::shared_ptr<siconos::simulation::TimeDiscretisation> td,
    std::shared_ptr<siconos::modeling::DynamicalSystem> ds,
    const siconos::modeling::NonSmoothDynamicalSystem& nsds) {
  auto s = *(_allSensors.insert(SensorFactory::instance()->create(ds, type))).first;
  linkSensorSimulation(s, td);
  s->initialize(nsds);
  return s;
}

std::shared_ptr<siconos::control::Actuator> siconos::control::ControlManager::addActuator(
    ActuatorType type, std::shared_ptr<siconos::simulation::TimeDiscretisation> td,
    std::shared_ptr<ControlSensor> sensor) {
  if (!sensor)
    THROW_EXCEPTION("siconos::control::ControlManager::addActuator - sensor is not valid !");
  auto act = *(_allActuators.insert(ActuatorFactory::instance()->create(sensor, type))).first;
  linkActuatorSimulation(act, td);
  return act;
}

std::shared_ptr<siconos::control::Actuator>
siconos::control::ControlManager::addAndRecordActuator(
    ActuatorType type, std::shared_ptr<siconos::simulation::TimeDiscretisation> td,
    std::shared_ptr<ControlSensor> sensor,
    const siconos::modeling::NonSmoothDynamicalSystem& nsds) {
  if (!sensor)
    THROW_EXCEPTION("siconos::control::ControlManager::addActuator - sensor is not valid !");
  auto act = *(_allActuators.insert(ActuatorFactory::instance()->create(sensor, type))).first;
  linkActuatorSimulation(act, td);
  act->initialize(nsds, *_sim);
  return act;
}

std::shared_ptr<siconos::control::Observer> siconos::control::ControlManager::addObserver(
    ObserverType type, std::shared_ptr<siconos::simulation::TimeDiscretisation> td,
    std::shared_ptr<ControlSensor> sensor, const siconos::algebra::SiconosVector& xHat0) {
  if (!sensor)
    THROW_EXCEPTION("siconos::control::ControlManager::addActuator - sensor is not valid !");
  auto obs =
      *(_allObservers.insert(ObserverFactory::instance()->create(sensor, xHat0, type))).first;
  linkObserverSimulation(obs, td);
  return obs;
}

std::shared_ptr<siconos::control::Observer>
siconos::control::ControlManager::addAndRecordObserver(
    ObserverType type, std::shared_ptr<siconos::simulation::TimeDiscretisation> td,
    std::shared_ptr<ControlSensor> sensor, const siconos::algebra::SiconosVector& xHat0,
    const siconos::modeling::NonSmoothDynamicalSystem& nsds) {
  auto obs =
      *(_allObservers.insert(ObserverFactory::instance()->create(sensor, xHat0, type))).first;
  linkObserverSimulation(obs, td);
  obs->initialize(nsds, *_sim);
  return obs;
}

void siconos::control::ControlManager::addSensorPtr(
    std::shared_ptr<Sensor> s, std::shared_ptr<siconos::simulation::TimeDiscretisation> td) {
  if (!s)
    THROW_EXCEPTION("siconos::control::ControlManager::addActuator - sensor is not valid !");
  _allSensors.insert(s);
  linkSensorSimulation(s, td);
}

void siconos::control::ControlManager::addAndRecordSensorPtr(
    std::shared_ptr<Sensor> s, std::shared_ptr<siconos::simulation::TimeDiscretisation> td,
    const siconos::modeling::NonSmoothDynamicalSystem& nsds) {
  if (!s)
    THROW_EXCEPTION("siconos::control::ControlManager::addActuator - sensor is not valid !");
  _allSensors.insert(s);
  linkSensorSimulation(s, td);
  s->initialize(nsds);
}

void siconos::control::ControlManager::addActuatorPtr(
    std::shared_ptr<Actuator> act,
    std::shared_ptr<siconos::simulation::TimeDiscretisation> td) {
  if (!act)
    THROW_EXCEPTION("siconos::control::ControlManager::addActuator - actuator is not valid !");
  _allActuators.insert(act);
  linkActuatorSimulation(act, td);
}

void siconos::control::ControlManager::addAndRecordActuatorPtr(
    std::shared_ptr<Actuator> act, std::shared_ptr<siconos::simulation::TimeDiscretisation> td,
    const siconos::modeling::NonSmoothDynamicalSystem& nsds) {
  if (!act)
    THROW_EXCEPTION("siconos::control::ControlManager::addActuator - actuator is not valid !");
  _allActuators.insert(act);
  linkActuatorSimulation(act, td);
  act->initialize(nsds, *_sim);
}

void siconos::control::ControlManager::addObserverPtr(
    std::shared_ptr<Observer> obs,
    std::shared_ptr<siconos::simulation::TimeDiscretisation> td) {
  if (!obs)
    THROW_EXCEPTION("siconos::control::ControlManager::addActuator - observer is not valid !");
  _allObservers.insert(obs);
  linkObserverSimulation(obs, td);
}

void siconos::control::ControlManager::addAndRecordObserverPtr(
    std::shared_ptr<Observer> obs, std::shared_ptr<siconos::simulation::TimeDiscretisation> td,
    const siconos::modeling::NonSmoothDynamicalSystem& nsds) {
  if (!obs)
    THROW_EXCEPTION("siconos::control::ControlManager::addActuator - observer is not valid !");
  _allObservers.insert(obs);
  linkObserverSimulation(obs, td);
  obs->initialize(nsds, *_sim);
}

void siconos::control::ControlManager::linkSensorSimulation(
    std::shared_ptr<Sensor> s, std::shared_ptr<siconos::simulation::TimeDiscretisation> td) {
  auto& ev = _sim->eventsManager()->insertEvent(siconos::simulation::EventType::Sensor, td);
  static_cast<SensorEvent&>(ev).setSensorPtr(s);
  s->setTimeDiscretisation(*td);
}

void siconos::control::ControlManager::linkActuatorSimulation(
    std::shared_ptr<Actuator> act,
    std::shared_ptr<siconos::simulation::TimeDiscretisation> td) {
  auto& ev = _sim->eventsManager()->insertEvent(siconos::simulation::EventType::Actuator, td);
  static_cast<ActuatorEvent&>(ev).setActuatorPtr(act);
  act->setTimeDiscretisation(*td);
}

void siconos::control::ControlManager::linkObserverSimulation(
    std::shared_ptr<Observer> obs,
    std::shared_ptr<siconos::simulation::TimeDiscretisation> td) {
  auto& ev = _sim->eventsManager()->insertEvent(siconos::simulation::EventType::Observer, td);
  static_cast<ObserverEvent&>(ev).setObserverPtr(obs);
  obs->setTimeDiscretisation(*td);
}

void siconos::control::ControlManager::display() const {
  std::cout << "=========> ControlManager ";
  std::cout << "It handles the following objects: " << std::endl;
  for (auto itS : _allSensors) itS->display();

  for (auto itA : _allActuators) itA->display();

  for (auto itO : _allObservers) itO->display();
  std::cout << "==========" << std::endl;
  std::cout << std::endl;
}
