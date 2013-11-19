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
#include<iostream>

class SiconosVector;


namespace ioVector
{

  typedef std::ios_base::openmode openmode;
  /** Format to read binary data */ 
  const openmode BINARY_IN = std::ios::in|std::ios::binary;
  /** Format to write binary data */ 
  const openmode BINARY_OUT = std::ios::out|std::ios::binary; //|std::ios::scientific;
  // Note FP : ios::scientific result in file.good() always false, even if writing process goes well ...
  // Don't know why.
  /** Format to read ascii data */ 
  const openmode ASCII_IN = std::ios::in;
  /** Format to write ascii data */ 
  const openmode ASCII_OUT = std::ios::out|std::ios::scientific;
    

/** Read a SiconosVector from a file
    \param[in] fileName the file containing the vector
    \param[in] ios_base::openmode, mode for reading (like  ios::in|ios:binary ...)
    default = ascii
    \param[in,out] m the SiconosVector to be filled
    \param[in] inputType (see outputType in write function)
    \return bool true if read ok, else false ...
 */
  bool read(const std::string& fileName, SiconosVector& m, const openmode& = ASCII_IN, const std::string& inputType = "python");
  
/** Write a SiconosVector to a file
    \param[in] output file name
    \param[in] ios_base::openmode, mode for writing (like  ios::out|ios:binary ...)
    default = ascii
    \param[in,out] m the SiconosVector to be written
    \param[in] std::string type of output:
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
  bool write(const std::string& fileName, const SiconosVector& m, const openmode& = ASCII_OUT, const std::string& outputType = "python");
}
#endif
