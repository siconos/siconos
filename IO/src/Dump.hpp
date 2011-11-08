
#include "SiconosFull.hpp"

namespace Siconos
{
namespace IO
{

/** save Siconos model into a file
 *  \param model
 *  \param filename with extension : .xml, .dat
 */
void save(SP::Model model, std::string filename);

/** load Siconos model from file
 * \param filename
 * \return a SP::Model
 */
SP::Model load(std::string filename);
}
}
