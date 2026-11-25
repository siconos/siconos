// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifdef WITH_SERIALIZATION
#include "SiconosRestart.hpp"

#include <boost/filesystem.hpp>

#include "RegisterSimulation.hpp"
#include "Simulation.hpp"

namespace siconos {

void save(std::shared_ptr<siconos::simulation::Simulation> s, const std::string& filename) {
  boost::filesystem::path tempf = boost::filesystem::path(filename + ".tmp");

  boost::filesystem::path destf = boost::filesystem::path(filename);

  std::ofstream ofs(tempf.c_str());
  {
    if (destf.extension() == ".xml") {
      RegisterSimulationOxml(ofs, s);
    }

    else if (destf.extension() == ".bin") {
      RegisterSimulationObin(ofs, s);
    }
  }

  // atomic
  boost::filesystem::rename(tempf, destf);
}

/** load Siconos model from file
 * \param filename
 */
std::shared_ptr<siconos::simulation::Simulation> load(const std::string& filename) {
  std::shared_ptr<siconos::simulation::Simulation> s;

  std::ifstream ifs(filename.c_str());
  {
    if (boost::filesystem::path(filename).extension() == ".xml") {
      RegisterSimulationIxml(ifs, s);
    } else if (boost::filesystem::path(filename).extension() == ".bin") {
      RegisterSimulationIbin(ifs, s);
    }
  }
  return s;
}
}  // namespace siconos
#else
// #include "SiconosRestart.hpp"
#include "SiconosException.hpp"
#include "Simulation.hpp"
namespace siconos {

void save(std::shared_ptr<siconos::simulation::Simulation> s, const std::string& filename) {
  THROW_EXCEPTION("Siconos/IO must be compiled with serialization support for this service.");
}

std::shared_ptr<siconos::simulation::Simulation> load(const std::string& filename) {
  THROW_EXCEPTION("Siconos/IO must be compiled with serialization support for this service.");
  /* Dummy return to make every compiler happy  */
  return std::shared_ptr<siconos::simulation::Simulation>();
}
}  // namespace siconos
#endif
