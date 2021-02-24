/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

#include "SiconosConfig.h"
#ifdef WITH_SERIALIZATION
#include "SiconosRestart.hpp"

#include <boost/filesystem.hpp>

#include "Simulation.hpp"
#include "RegisterSimulation.hpp"


namespace Siconos
{

void save(SP::Simulation s, const std::string& filename)
{
  boost::filesystem::path tempf =
    boost::filesystem::path(filename + ".tmp");

  boost::filesystem::path destf =
    boost::filesystem::path(filename);

  std::ofstream ofs(tempf.c_str());
  {
    if(destf.extension() == ".xml")
    {
      RegisterSimulationOxml(ofs, s);
    }

    else if(destf.extension() == ".bin")
    {
      RegisterSimulationObin(ofs, s);
    }
  }

  // atomic
  boost::filesystem::rename(tempf, destf);
}

/** load Siconos model from file
 * \param filename
 */
SP::Simulation load(const std::string& filename)
{
  //SP::Simulation s(new Simulation());
  SP::Simulation s;

  std::ifstream ifs(filename.c_str());
  {
    if(boost::filesystem::path(filename).extension() == ".xml")
    {
      RegisterSimulationIxml(ifs, s);
    }
    else if(boost::filesystem::path(filename).extension() == ".bin")
    {
      RegisterSimulationIbin(ifs, s);
    }
  }
  return s;
}
}
#else
#include "SiconosRestart.hpp"
#include <SiconosException.hpp>
namespace Siconos
{

void save(SP::Simulation s, const std::string& filename)
{
  THROW_EXCEPTION("Siconos/IO must be compiled with serialization support for this service.");
}

SP::Simulation load(const std::string& filename)
{
  THROW_EXCEPTION("Siconos/IO must be compiled with serialization support for this service.");
  /* Dummy return to make every compiler happy  */
  return std::shared_ptr<Simulation>();
}
}
#endif
