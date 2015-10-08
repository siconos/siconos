/* Siconos-IO, Copyright INRIA 2005-2011.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

#include "SiconosConfig.h"
#ifdef WITH_SERIALIZATION
#include "SiconosRestart.hpp"

#include <boost/filesystem.hpp>

#include "RegisterModel.hpp"

#include <Model.hpp>

namespace Siconos
{

  void save(SP::Model model, std::string filename)
  {
    boost::filesystem::path tempf =
      boost::filesystem::path(filename + ".tmp");

    boost::filesystem::path destf =
      boost::filesystem::path(filename);

    std::ofstream ofs(tempf.c_str());
    {
      if (destf.extension() == ".xml")
      {
        RegisterModelOxml(ofs, model);
      }

      else if (destf.extension() == ".bin")
      {
        RegisterModelObin(ofs, model);
      }
    }

    // atomic
    boost::filesystem::rename(tempf, destf);
  }

/** load Siconos model from file
 * \param filename
 */
  SP::Model load(std::string filename)
  {
    SP::Model model(new Model());

    std::ifstream ifs(filename.c_str());
    {
      if (boost::filesystem::path(filename).extension() == ".xml")
      {
        RegisterModelIxml(ifs, model);
      }
      else if (boost::filesystem::path(filename).extension() == ".bin")
      {
        RegisterModelIbin(ifs, model);
      }
    }
    return model;
  }
}
#else
#include "SiconosRestart.hpp"
#include <RuntimeException.hpp>
namespace Siconos
{

  void save(SP::Model model, std::string filename)
  {
    RuntimeException::selfThrow("Siconos/IO must be compiled with serialization support for this service.");
  }

  SP::Model load(std::string filename)
  {
    RuntimeException::selfThrow("Siconos/IO must be compiled with serialization support for this service.");
  }
}
#endif
