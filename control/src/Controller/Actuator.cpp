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

#include "Actuator.hpp"
#include "SimpleMatrix.hpp"
#include "Topology.hpp"
#include "SiconosVector.hpp"
#include "ControlSensor.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "Tools.hpp"
#include "PluggedObject.hpp"
#include "SiconosException.hpp"

siconos::control::Actuator::Actuator(ActuatorType type, std::shared_ptr<ControlSensor> sensor)
    : _type(type), _sensor(sensor)
{
}

siconos::control::Actuator::Actuator(ActuatorType type, std::shared_ptr<ControlSensor> sensor,
                                     std::shared_ptr<siconos::algebra::SimpleMatrix> B)
    : _type(type), _B(B), _sensor(sensor)
{
  if (B) {
    _u = std::make_shared<siconos::algebra::SiconosVector>(B->size(1));
    _u->setZero();
  }
}

void siconos::control::Actuator::addSensorPtr(std::shared_ptr<ControlSensor> newSensor)
{
  _sensor = newSensor;
}

void siconos::control::Actuator::initialize(
    const siconos::modeling::NonSmoothDynamicalSystem& nsds,
    const siconos::simulation::Simulation& s)
{
  if (!_sensor) {
    THROW_EXCEPTION(
        "siconos::control::Actuator::initialize - No Sensor given to the Actuator");
  }

  // Init the control variable and add the necessary properties
  auto& DSG0 = *nsds.topology()->dSG(0);
  auto dsgVD = DSG0.descriptor(_sensor->getDS());
  if (_B) {
    DSG0.B[dsgVD] = _B;
  }
  else if (!_plugingName.empty()) {
    DSG0.pluginU[dsgVD] = std::make_shared<siconos::plugins::PluggedObject>(_plugingName);
    if (!_pluginJacgxName.empty()) {
      DSG0.pluginJacgx[dsgVD] =
          std::make_shared<siconos::plugins::PluggedObject>(_plugingName);
    }
    if (!_u) {
      THROW_EXCEPTION(
          "siconos::control::Actuator::initialize - u should have already been initialized");
    }
  }
  else {
    THROW_EXCEPTION(
        "siconos::control::Actuator::initialize - neither the matrix B or the plugin g are "
        "not initialized");
  }

  DSG0.u[dsgVD] = _u;
}

void siconos::control::Actuator::setSizeu(unsigned size)
{
  if (_B && size != _B->size(1)) {
  }
  _u= std::make_shared<siconos::algebra::SiconosVector>(size);
  _u->setZero();
}

std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> siconos::control::Actuator::getInternalNSDS() const
{
  return std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem>();
}

void siconos::control::Actuator::display() const
{
  std::cout << "=====> Actuator of type " << siconos::tools::enum_to_string(_type) << ", named " << _id << std::endl;
  ;
  std::cout << "The associated Sensor is: " << std::endl;
  if (_sensor) _sensor->display();
  std::cout << "======" << std::endl;
  std::cout << "The value of the control is: " << std::endl;
  _u->display();
  std::cout << std::endl;
}
