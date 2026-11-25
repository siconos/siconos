// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "Sensor.hpp"

#include <iostream>

#include "DynamicalSystem.hpp"
#include "Tools.hpp"  // enum_to_string

siconos::control::Sensor::Sensor(SensorType type,
                                 std::shared_ptr<siconos::modeling::DynamicalSystem> ds)
    : _type(type), _id("none"), _DS(ds) {
  _DSx = _DS->x();
}

void siconos::control::Sensor::display() const {
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
