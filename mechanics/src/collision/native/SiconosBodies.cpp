/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include "SiconosBodies.hpp"

#include "SiconosException.hpp"
#include "Simulation.hpp"

void siconos::collision::native::SiconosBodies::compute() {
  try {
    _sim->advanceToEvent();
    _sim->processEvents();
  }

  catch (...) {
    siconos::exception::process();
  }
}
