// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "NormalConeNSL.hpp"

#include <cassert>
#include <iostream>

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

siconos::modeling::NormalConeNSL::NormalConeNSL(
    siconos::algebra::Index size, std::shared_ptr<siconos::algebra::SiconosMatrix> H,
    std::shared_ptr<siconos::algebra::SiconosVector> K)
    : NonSmoothLaw(size), _H(H), _K(K) {
  assert(H->cols() == size &&
         "NormalConeNSL::NormalConeNSL - the number of columns in H and the declared size are "
         "not equal, check your code !");
}

bool siconos::modeling::NormalConeNSL::isVerified(void) const {
  bool res = false;
  // to do
  return res;
}

void siconos::modeling::NormalConeNSL::display() const {
  std::cout << "------------------------------------" << std::endl;
  std::cout << "____ data of the NormalConeNSL" << std::endl;
  std::cout << "| nSLawSize : " << _size << std::endl;
  std::cout << "| H : " << std::endl;
  siconos::algebra::print(*_H);
  std::cout << "| K : " << std::endl;
  siconos::algebra::print(*_K);
  std::cout << "____________________________" << std::endl;
  std::cout << "------------------------------------" << std::endl;
}
