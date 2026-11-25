// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "VonMisesPlasticityNSL.hpp"

#include <iostream>

#include "SiconosException.hpp"

bool siconos::modeling::VonMisesPlasticityNSL::isVerified() const {
  THROW_EXCEPTION("VonMisesPlasticityNSL::isVerified, not yet implemented!");
  return false;
}

void siconos::modeling::VonMisesPlasticityNSL::display() const {
  std::cout << "=== Von Mises plasticity non-smooth law ===" << std::endl;
  std::cout << " yield stress: " << _sigma_y << std::endl;
  std::cout << " hypothesis: "
            << (_hypothesis == Hypothesis::plane_stress ? "plane stress" : "plane strain")
            << std::endl;
  std::cout << "===========================================" << std::endl;
}
