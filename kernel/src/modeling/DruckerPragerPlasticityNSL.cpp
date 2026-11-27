// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "DruckerPragerPlasticityNSL.hpp"

#include <cmath>
#include <iostream>
#include <memory>

#include "SiconosException.hpp"

std::shared_ptr<siconos::modeling::DruckerPragerPlasticityNSL>
siconos::modeling::DruckerPragerPlasticityNSL::mohrCoulombPlaneStrain(double c, double phi,
                                                                      double psi) {
  return std::make_shared<DruckerPragerPlasticityNSL>(
      std::sin(phi), std::sin(phi) - std::sin(psi), std::sqrt(2.0) * c * std::cos(phi));
}

bool siconos::modeling::DruckerPragerPlasticityNSL::isVerified() const {
  THROW_EXCEPTION("DruckerPragerPlasticityNSL::isVerified, not yet implemented!");
  return false;
}

void siconos::modeling::DruckerPragerPlasticityNSL::display() const {
  std::cout << "=== Drucker-Prager plasticity non-smooth law ===" << std::endl;
  std::cout << " eta (slope): " << _eta << std::endl;
  std::cout << " theta (dilatancy parameter): " << _theta << std::endl;
  std::cout << " cohesion: " << _cohesion << std::endl;
  std::cout << "================================================" << std::endl;
}
