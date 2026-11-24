// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file SiconosRestart.hpp
  \brief provides pre-compiled functions for a full Siconos Model
  serialization through libSiconosIO library */

#ifndef SICONOSRESTART_HPP
#define SICONOSRESTART_HPP

#include <SiconosFwd.hpp>
#include <string>

/** SICONOS
 */
namespace siconos {

/** save a Siconos Simulation with the full simulation state into a
 *  file
 * \param s a Simulation
 * \param filename with extension : .xml, .bin (binary archive)
 */
void save(std::shared_ptr<siconos::simulation::Simulation> s, const std::string& filename);

/** load a Siconos Model with the full simulation state from file
 * \param filename
 * \return a Simulation
 */
std::shared_ptr<siconos::simulation::Simulation> load(const std::string& filename);

}  // namespace siconos

#endif
