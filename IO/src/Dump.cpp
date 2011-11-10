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

#include "Dump.hpp"

#include <boost/filesystem.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>


namespace Siconos
{
namespace IO
{

/** save Siconos model into a file
 *  \param model
 *  \param filename with extension : .xml, .dat
 */
void save(SP::Model model, std::string filename)
{
  std::ofstream ofs(filename.c_str());
  {
    if (boost::filesystem::path(filename).extension() == ".xml")
    {
      boost::archive::xml_oarchive oa(ofs);
      siconos_io_register(oa);
      oa << NVP(model);
    }

    else if (boost::filesystem::path(filename).extension() == ".dat")
    {
      boost::archive::binary_oarchive oa(ofs);
      siconos_io_register(oa);
      oa << NVP(model);  // remove NVP
    }
  }
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
      boost::archive::xml_iarchive ia(ifs);
      siconos_io_register(ia);
      ia >> NVP(model);
      return model;
    }
    else if (boost::filesystem::path(filename).extension() == ".dat")
    {
      boost::archive::binary_iarchive ia(ifs);
      siconos_io_register(ia);
      ia >> NVP(model);
      return model;

    }
  }
}
}
}
