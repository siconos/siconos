// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifdef WITH_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>

#include "RegisterSimulation.hpp"
#include "SiconosFull.hpp"

void RegisterSimulationObin(std::ofstream& ofs,
                            std::shared_ptr<siconos::simulation::Simulation>& sim) {
  boost::archive::binary_oarchive ar(ofs);
  siconos_io_register_Numerics(ar);
  siconos_io_register_Kernel(ar);
  siconos_io_register_Mechanics(ar);
  siconos_io_register_Control(ar);
  ar << NVP(sim);
}
#endif
