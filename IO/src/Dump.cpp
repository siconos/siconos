

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
