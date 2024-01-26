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

#include "Sensor.hpp"

#include <iostream>

#include "DynamicalSystem.hpp"
#include "Tools.hpp"  // enum_to_string

siconos::control::Sensor::Sensor(SensorType type,
                                 std::shared_ptr<siconos::modeling::DynamicalSystem> ds)
    : _type(type), _id("none"), _DS(ds)
{
  _DSx = _DS->x();
}

void siconos::control::Sensor::display() const
{
  std::cout << "=====> Sensor of type " << siconos::tools::enum_to_string(_type) << ", named "
            << _id;
  if (_DS)
    std::cout << " and linked to the DynamicalSystem number " << _DS->number() << "."
              << std::endl;
  else
    std::cout << " and not linked to a DynamicalSystem." << std::endl;
  std::cout << "======" << std::endl;
  std::cout << std::endl;
}
