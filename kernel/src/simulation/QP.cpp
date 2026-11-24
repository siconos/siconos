// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "QP.hpp"

#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

int siconos::nonsmooth_formulations::QP::compute(double) {
  THROW_EXCEPTION("siconos::nonsmooth_formulations::QP::compute not yet implemented");
  return 1;
}

void siconos::nonsmooth_formulations::QP::display() const {
  THROW_EXCEPTION("siconos::nonsmooth_formulations::QP::compute not yet implemented");
}

void siconos::nonsmooth_formulations::QP::setQ(
    const siconos::algebra::SiconosMatrix& newValue) {
  *_Q = newValue;
}
void siconos::nonsmooth_formulations::QP::setP(
    const siconos::algebra::SiconosVector& newValue) {
  *_p = newValue;
}
