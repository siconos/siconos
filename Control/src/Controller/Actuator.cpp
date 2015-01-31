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

#include "RuntimeException.hpp"
#include "Actuator.hpp"
#include "ActuatorEvent.hpp"
#include "ControlSensor.hpp"
#include "Model.hpp"
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

Actuator::~Actuator()
{
}

void Actuator::addSensorPtr(SP::ControlSensor newSensor)
{
  _sensor = newSensor;
}

void Actuator::initialize(const Model& m)
{
  if (!_sensor)
  {
    RuntimeException::selfThrow("Actuator::initialize - No Sensor given to the Actuator");
  }

  // Init the control variable and add the necessary properties
  DynamicalSystemsGraph& DSG0 = *m.nonSmoothDynamicalSystem()->topology()->dSG(0);
  DynamicalSystemsGraph::VDescriptor dsgVD = DSG0.descriptor(_sensor->getDS());
  if (_B)
  {
    DSG0.B[dsgVD] = _B;
  }
  else if (!_plugingName.empty())
  {
    DSG0.pluginU[dsgVD].reset(new PluggedObject(_plugingName));
    if (!_pluginJacgxName.empty())
    {
      DSG0.pluginJacgx[dsgVD].reset(new PluggedObject(_plugingName));
    }
    if (!_u)
    {
      RuntimeException::selfThrow("Actuator::initialize - u should have already been initialized");
    }
  }
  else
  {
    RuntimeException::selfThrow("Actuator::initialize - niether the matrix B or the plugin g are not initialized");
  }

  DSG0.u[dsgVD] = _u;
}

void Actuator::setSizeu(unsigned size)
{
  _u.reset(new SiconosVector(size, 0));
}

SP::Model Actuator::getInternalModel() const
{
  return std11::shared_ptr<Model>();
}

void Actuator::display() const
{
  std::cout << "=====> Actuator of type " << _type << ", named " << _id << std::endl;;
  std::cout << "The associated Sensor is: " << std::endl;
  if (_sensor)
    _sensor->display();
  std::cout << "======" <<std::endl;
  std::cout << "The value of the control is: " << std::endl;
  _u->display();
  std::cout <<std::endl;
}
