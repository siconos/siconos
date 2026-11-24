// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "NewtonImpactRollingFrictionNSL.hpp"

#include <iostream>

#include "SiconosException.hpp"

bool siconos::modeling::NewtonImpactRollingFrictionNSL::isVerified() const {
  bool res = false;
  // to do
  THROW_EXCEPTION("NewtonImpactRollingFrictionNSL:: isVerified, not yet implemented!");
  return res;
}

void siconos::modeling::NewtonImpactRollingFrictionNSL::display() const {
  std::cout << "=== Newton impact-friction non-smooth law data display ===" << std::endl;
  std::cout << " Normal Newton coefficient of restitution: " << _en << std::endl;
  std::cout << " Tangential Newton coefficient of restitution: " << _et << std::endl;
  std::cout << "Friction coefficient: " << _mu << std::endl;
  std::cout << "Rolling friction coefficient: " << _muR << std::endl;
  std::cout << "==========================================================" << std::endl;
}
