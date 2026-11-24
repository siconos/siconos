/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include "FremondImpactFrictionNSL.hpp"

#include <iostream>

void siconos::modeling::FremondImpactFrictionNSL::display() const {
  std::cout << "=== Fremond impact-friction non-smooth law data display ===\n";
  std::cout << " Normal Fremond coefficient of restitution: " << _en << "\n";
  std::cout << " Tangential Fremond coefficient of restitution: " << _et << "\n";
  std::cout << "Friction coefficient: " << _mu << "\n";
  std::cout << "==========================================================\n";
}
