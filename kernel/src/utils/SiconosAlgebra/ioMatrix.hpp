/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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

/*! \file ioMatrix.hpp
    \brief input/output for SiconosMatrix

*/

#ifndef __ioMatrix__
#define __ioMatrix__

#include <string>
#include "SiconosFwd.hpp"
#include "SiconosAlgebraTypeDef.hpp"

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

/** Function to load data from a file and compare it with the provided
 * data.  Returns true if the file was loaded and the comparison was
 * performed.  Caller needs to check diff <= epsilon to verify the
 * result.
 * \param data The data to compare against the file.
 * \param filename The name of the file to load and compare.
 * \param epsilon The comparison threshold.
 * \param error An optional pointer to double to receive scalar error.
 * \param index An optional list of column indexes, size==0 indicates all columns.
 * \param ref If provided, loaded matrix is returned in this pointer.
 * \param mode Mode string to pass to ioMatrix::read.
 * \param verbose True to print verbose output.
 * \return True if the file was loaded and the comparison was performed.
 */
bool compareRefFile(const SimpleMatrix& data, std::string filename, double epsilon,
                    double *error=0, Index index=Index(), SP::SimpleMatrix *ref=0,
                    std::string mode="ascii", bool verbose=true);
}

#endif
