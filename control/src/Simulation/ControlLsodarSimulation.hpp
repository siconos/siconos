// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*!\file ControlLsodarSimulation.hpp
 * \brief Simulation with Control using Lsodar
 */

#ifndef ControlLsodarSimulation_hpp
#define ControlLsodarSimulation_hpp

#include "ControlSimulation.hpp"

namespace siconos::control {
class ControlLsodarSimulation : public ControlSimulation {
 public:
  ControlLsodarSimulation(double t0, double T, double h);

  virtual void run();
};
}  // namespace siconos::control
#endif
