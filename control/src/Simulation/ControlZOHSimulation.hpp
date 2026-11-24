// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*!\file ControlZOHSimulation.hpp
 * \brief Simulation with Control using ZOH integration scheme
 */

#ifndef ControlZOHSimulation_hpp
#define ControlZOHSimulation_hpp

#include "ControlSimulation.hpp"

namespace siconos::control {
class ControlZOHSimulation : public ControlSimulation {
 public:
  ControlZOHSimulation(double t0, double T, double h);

  virtual void run();
};
}  // namespace siconos::control

#endif
