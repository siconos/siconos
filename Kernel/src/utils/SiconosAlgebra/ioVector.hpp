/* Siconos-Kernel, Copyright INRIA 2005-2012.
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

/*! \file ioVector.hpp
   input/output for SiconosVector

*/

#ifndef __ioVector__
#define __ioVector__

#include <string>

class SiconosVector;


namespace ioVector
{
/** Read a SiconosMatrix from a file
 * \param fileName the file containing the matrix
 * \param[in,out] m the SiconosVector containing the matrix
 * \return bool true if read ok, else false ...
 */
bool read(const std::string& fileName, const std::string& Mode, SiconosVector& m);

/** Write a SiconosVector to a file
    \param[in] SiconosVector the vector to be read
    \param[in] string type of output:
    Type of Output for write function:
    - "boost": boost way: \n
    [row] (a0, a1,..)
    - "python"(default): \n
    row \n
    a0 a1 a2 ... \n
    - "noDim": \n
    a0 a1 a2 ... \n
    Reading input format is the one corresponding to "python".
    \return bool true if read ok, else false ...
*/
bool write(const std::string& fileName, const std::string& Mode, const SiconosVector& m, const std::string& outputType = "python");
}
#endif
