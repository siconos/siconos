/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#ifndef RegisterModel_hpp
#define RegisterModel_hpp

#ifdef WITH_SERIALIZATION

#include <SiconosFwd.hpp>
#include <fstream>

void RegisterSimulationOxml(std::ofstream& ofs,
                            std::shared_ptr<siconos::simulation::Simulation>&);
void RegisterSimulationObin(std::ofstream& ofs,
                            std::shared_ptr<siconos::simulation::Simulation>&);
void RegisterSimulationIxml(std::ifstream& ifs,
                            std::shared_ptr<siconos::simulation::Simulation>&);
void RegisterSimulationIbin(std::ifstream& ifs,
                            std::shared_ptr<siconos::simulation::Simulation>&);

#endif

#endif
