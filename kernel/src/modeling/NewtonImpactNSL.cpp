// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "NewtonImpactNSL.hpp"

#include <iostream>

#include "SiconosException.hpp"

bool siconos::modeling::NewtonImpactNSL::isVerified() const {
  bool res = false;
  // to do
  THROW_EXCEPTION("NewtonImpactFrictionNSL:: isVerified, not yet implemented!");
  return res;
}

void siconos::modeling::NewtonImpactNSL::display() const {
  std::cout
      << "==============================================================================="
      << std::endl;
  std::cout << "=== Newton impact (frictionless) non-smooth law coefficient of restitution: "
            << _e << std::endl;
  std::cout
      << "==============================================================================="
      << std::endl;
}
