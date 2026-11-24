/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#include <iostream>

#include "MohrCoulombPlasticityNSL.hpp"
#include "SiconosException.hpp"

bool siconos::modeling::MohrCoulombPlasticityNSL::isVerified() const {
  bool res = false;
  // to do
  THROW_EXCEPTION("MohrCoulombPlasticityNSL:: isVerified, not yet implemented!");
  return res;
}

void siconos::modeling::MohrCoulombPlasticityNSL::display() const {
  std::cout << "=== Mohr Coulomb Plasticity non-smooth law data display ===" << std::endl;
  std::cout << " Cohesion: " << _c << std::endl;
  std::cout << " Friction Angle: " << _phi << std::endl;
  std::cout << "==========================================================" << std::endl;
}
