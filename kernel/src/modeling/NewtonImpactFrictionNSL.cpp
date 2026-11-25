/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#include "NewtonImpactFrictionNSL.hpp"

#include <iostream>

#include "SiconosException.hpp"

bool siconos::modeling::NewtonImpactFrictionNSL::isVerified() const {
  bool res = false;
  // to do
  THROW_EXCEPTION("NewtonImpactFrictionNSL:: isVerified, not yet implemented!");
  return res;
}

void siconos::modeling::NewtonImpactFrictionNSL::display() const {
  std::cout << "=== Newton impact-friction non-smooth law data display ===" << std::endl;
  std::cout << " Normal Newton coefficient of restitution: " << _en << std::endl;
  std::cout << " Tangential Newton coefficient of restitution: " << _et << std::endl;
  std::cout << "Friction coefficient: " << _mu << std::endl;
  std::cout << "==========================================================" << std::endl;
}
