// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "RelayNSL.hpp"

#include <iostream>

siconos::modeling::RelayNSL::RelayNSL(siconos::algebra::Index size, double lb, double ub)
    : NonSmoothLaw(size), _lb(lb), _ub(ub) {
  if (ub < lb)
    THROW_EXCEPTION("RelayNSL::RelayNSL). Upper bound ub must be greater than lower bound lb");
}

bool siconos::modeling::RelayNSL::RelayNSL::isVerified(void) const {
  bool res = false;
  // to do
  return res;
}

void siconos::modeling::RelayNSL::RelayNSL::display() const {
  std::cout << "------------------------------------" << std::endl;
  std::cout << "____ data of the RelayNSL" << std::endl;
  std::cout << "| nSLawSize : " << _size << std::endl;
  std::cout << "| lb : " << _lb << std::endl;
  std::cout << "| ub : " << _ub << std::endl;
  std::cout << "____________________________" << std::endl;
  std::cout << "------------------------------------" << std::endl;
}
