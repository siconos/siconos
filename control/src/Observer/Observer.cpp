// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "Observer.hpp"

#include "SiconosException.hpp"
#include "SiconosVector.hpp"
#include "TimeDiscretisation.hpp"
#include "Tools.hpp"  // enum_to_string

siconos::control::Observer::Observer(ObserverType type, std::shared_ptr<ControlSensor> sensor,
                                     const siconos::algebra::SiconosVector& xHat0,
                                     const std::string& newId)
    : _type(type), _sensor(sensor), _id(newId) {
  _xHat = std::make_shared<siconos::algebra::SiconosVector>(xHat0);
}

siconos::control::Observer::Observer(ObserverType type, std::shared_ptr<ControlSensor> sensor,
                                     const siconos::algebra::SiconosVector& xHat0,
                                     std::shared_ptr<siconos::modeling::DynamicalSystem> ds,
                                     const std::string& newId)
    : _type(type), _DS(ds), _sensor(sensor), _id(newId) {
  _xHat = std::make_shared<siconos::algebra::SiconosVector>(xHat0);
}

void siconos::control::Observer::initialize(
    const siconos::modeling::NonSmoothDynamicalSystem& nsds,
    const siconos::simulation::Simulation& s) {
  // Get the dimension of the output
  // XXX What if there is more than one sensor ...
  if (!_sensor) {
    THROW_EXCEPTION("siconos::control::Observer::initialize - the no ControlSensor was given");
  }
}

void siconos::control::Observer::display() const {
  std::cout << "=====> Observer of type " << siconos::tools::enum_to_string(_type)
            << ", named " << _id << "\n";
}

void siconos::control::Observer::setTimeDiscretisation(
    const siconos::simulation::TimeDiscretisation& td) {
  _td = std::make_shared<siconos::simulation::TimeDiscretisation>(td);
}
