// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "SecondOrderDS.hpp"

#include <iostream>

#include "BoundaryCondition.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

void siconos::modeling::SecondOrderDS::setBoundaryConditions(
    std::shared_ptr<siconos::modeling::BoundaryCondition> newbd) {
  if (boundaryConditions_) {
    std::cout << "Warning : SecondOrderDS::setBoundaryConditions. old boundary conditions "
                 "were pre-existing. \n";
  }
  boundaryConditions_ = newbd;
  reactionToBoundaryConditions_ =
      std::make_shared<siconos::algebra::SiconosVector>(boundaryConditions_->size());
};
