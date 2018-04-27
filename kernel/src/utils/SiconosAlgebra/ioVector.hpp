/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

/*! \file ioVector.hpp
   input/output for SiconosVector

*/

#ifndef __ioVector__
#define __ioVector__

#include <string>
#include <iostream>
#include "SiconosFwd.hpp"


namespace ioVector
{

  typedef std::ios_base::openmode openmode;
  /** Format to read binary data */ 
  const openmode BINARY_IN = std::ios::in|std::ios::binary;
  /** Format to write binary data */ 
  const openmode BINARY_OUT = std::ios::out|std::ios::binary; //|std::ios::scientific;
  // Note FP : ios::scientific result in file.good() always false, even if writing process goes well ...
  // Don't know why.
  // Note MB : ios::scientific -> fmtflags, std::ios:out -> openmode 
  /* stdc++ : 
   mode
    Flags describing the requested i/o mode for the file.
    This is an object of the bitmask member type openmode that consists of a combination of the following member constants:
    member constant	stands for	access
    in *	input	File open for reading: the internal stream buffer supports input operations.
    out	output	File open for writing: the internal stream buffer supports output operations.
    binary	binary	Operations are performed in binary mode rather than text.
    ate	at end	The output position starts at the end of the file.
    app	append	All output operations happen at the end of the file, appending to its existing contents.
    trunc	truncate	Any contents that existed in the file before it is open are discarded.
  */

  /** Format to read ascii data */ 
  const openmode ASCII_IN = std::ios::in;
  /** Format to write ascii data */ 
  const openmode ASCII_OUT = std::ios::out; //|std::ios::scientific;

/** Read a SiconosVector from a file
 *  \param[in] fileName the file containing the vector
 *  \param[in,out] m the SiconosVector to be filled
 *  \param[in] mode ios_base::openmode, mode for reading (like  ios::in|ios:binary ...)
 *       default = ascii
 *  \param[in] precision value for float output. Default = 15.
 *  \param[in] inputType (see outputType in write function)
 *  \param[in] flags for reading
 *  \return bool true if read ok, else false ...
 */
  bool read(const std::string& fileName, 
            SiconosVector& m, 
            const openmode&  mode = ASCII_IN,
            int precision =15,
            const std::string& inputType = "python",
            const std::ios::fmtflags& flags = std::cin.flags());
  
/** Write a SiconosVector to a file
    \param[in] fileName output file name
    \param[in] mode ios_base::openmode, mode for writing (like  ios::out|ios:binary ...)
    default = ascii
    \param[in] flags
    \param[in,out] m the SiconosVector to be written
    \param[in] precision value for float output. Default = 15.
    \param[in] outputType std::string type of output:
        Type of Output for write function:
        - "boost": boost way: \n
        [row] (a0, a1,..)
        - "python"(default): \n
        row \n
        a0 a1 a2 ... \n
        - "noDim": \n
        a0 a1 a2 ... \n
        Reading input format is the one corresponding to "python".
    \param[in] flags
    \return bool true if read ok, else false ...
*/
  bool write(const std::string& fileName,
             const SiconosVector& m, 
             const openmode& mode= ASCII_OUT,
             int precision =15, 
             const std::string& outputType = "python",
             const std::ios_base::fmtflags & flags  = std::cout.flags());
  
}
#endif
