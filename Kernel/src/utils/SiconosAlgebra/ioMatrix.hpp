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

/*! \file ioMatrix.hpp
    \brief input/output for SiconosMatrix

*/

#ifndef __ioMatrix__
#define __ioMatrix__

#include <string>
#include "SiconosFwd.hpp"

/** io object specialization */

namespace ioMatrix
{
/** Specialization to read a SiconosMatrix
    \param[in] fileName the name of the file to read
    \param[in] mode the storage type used in the file (either ascii or binary)
    \param[in,out] m the SiconosMatrix to be filled
    \return true if read ok, else false ...
*/
bool read(const std::string& fileName, const std::string& mode, SiconosMatrix& m);

/** Specialization to write a SiconosMatrix
    \param[in] fileName the name of the file to write in
    \param[in] mode the storage type used in the file (either ascii or binary)
    \param[in] m the SiconosMatrix to write
    \param[in] outputType type of output:
    - "python"(default):
    row col
    a00 a01 a02 ...
    a10 ...
    - "noDim":
    a00 a01 a02 ...
    a10 ...
    Reading input format is the one corresponding to "python".
    \return true if read ok, else false ...
*/
bool write(const std::string& fileName, const std::string& mode, const SiconosMatrix& m, const std::string& outputType = "python");

}

#endif
