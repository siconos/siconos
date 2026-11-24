// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "SensorEvent.hpp"

#include "Sensor.hpp"

void siconos::control::SensorEvent::process(siconos::simulation::Simulation& sim) {
  _sensor->capture();
}
