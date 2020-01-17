/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

#include "RuntimeException.hpp"
#include "Actuator.hpp"
#include "ActuatorEvent.hpp"
#include "ControlSensor.hpp"
#include "TimeDiscretisation.hpp"
#include "EventFactory.hpp"
#include "Simulation.hpp"
#include <iostream>
#include "NonSmoothDynamicalSystem.hpp"

Actuator::Actuator(): _type(0), _id("none")
{
}

Actuator::Actuator(unsigned int type, SP::ControlSensor sensor): _type(type), _id("none"), _sensor(sensor)
{
}

Actuator::Actuator(unsigned int type, SP::ControlSensor sensor, SP::SimpleMatrix B): _type(type), _id("none"), _B(B), _sensor(sensor)
{
  if(B)
  {
    _u = std11::make_shared<SiconosVector>(B->size(1), 0);
  }
}

Actuator::~Actuator()
{
}

void Actuator::addSensorPtr(SP::ControlSensor newSensor)
{
  _sensor = newSensor;
}

void Actuator::initialize(const NonSmoothDynamicalSystem& nsds, const Simulation& s)
{
  if(!_sensor)
  {
    RuntimeException::selfThrow("Actuator::initialize - No Sensor given to the Actuator");
  }

  // Init the control variable and add the necessary properties
  DynamicalSystemsGraph& DSG0 = *nsds.topology()->dSG(0);
  DynamicalSystemsGraph::VDescriptor dsgVD = DSG0.descriptor(_sensor->getDS());
  if(_B)
  {
    DSG0.B[dsgVD] = _B;
  }
  else if(!_plugingName.empty())
  {
    DSG0.pluginU[dsgVD].reset(new PluggedObject(_plugingName));
    if(!_pluginJacgxName.empty())
    {
      DSG0.pluginJacgx[dsgVD].reset(new PluggedObject(_plugingName));
    }
    if(!_u)
    {
      RuntimeException::selfThrow("Actuator::initialize - u should have already been initialized");
    }
  }
  else
  {
    RuntimeException::selfThrow("Actuator::initialize - neither the matrix B or the plugin g are not initialized");
  }

  DSG0.u[dsgVD] = _u;
}

void Actuator::setSizeu(unsigned size)
{
  if(_B && size != _B->size(1))
  {

  }
  _u.reset(new SiconosVector(size, 0));
}

SP::NonSmoothDynamicalSystem Actuator::getInternalNSDS() const
{
  return std11::shared_ptr<NonSmoothDynamicalSystem>();
}

void Actuator::display() const
{
  std::cout << "=====> Actuator of type " << _type << ", named " << _id << std::endl;;
  std::cout << "The associated Sensor is: " << std::endl;
  if(_sensor)
    _sensor->display();
  std::cout << "======" <<std::endl;
  std::cout << "The value of the control is: " << std::endl;
  _u->display();
  std::cout <<std::endl;
}
