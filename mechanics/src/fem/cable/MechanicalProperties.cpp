// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "MechanicalProperties.h"

#include <iostream>

void siconos::fem::cable::MechanicalProperties::set_T(double a_T) {
  if (a_T > 0.0) initialTension_ = a_T;
}

void siconos::fem::cable::MechanicalProperties::set_rho(double a_rho) {
  if (a_rho > 0.0) linearDensity_ = a_rho;
}

void siconos::fem::cable::MechanicalProperties::display() const {
  std::cout << "--- Mechanical properties: \n";
  std::cout << "- Rigidity: " << crossSectionRigidity_;
  std::cout << " - Linear density: " << linearDensity_;
  std::cout << " - Initial tension: " << initialTension_;
  std::cout << "\n --------------\n\n";
}
