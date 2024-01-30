/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

/*! \file io.hpp
   namespace to handle file input/output for vectors and matrices.

*/

#ifndef __ioVector__
#define __ioVector__

#include <iostream>
#include <nlohmann/json.hpp>  // json in/out
#include <vector>
// #include <string>
#include "SiconosVector.hpp"
#include "SiconosMatrix.hpp"
#include "SimpleMatrix.hpp"

/** utilities to handle file input/output for vectors and matrices */
namespace siconos::algebra {

// class SiconosMatrix;
// class SiconosVector;
// class SimpleMatrix;
namespace io {
/** Format to read binary data */
constexpr std::ios_base::openmode BINARY_IN = std::ios_base::in | std::ios_base::binary;
/** Format to write binary data */
constexpr std::ios_base::openmode BINARY_OUT = std::ios_base::out | std::ios_base::binary;
/** Format to read ascii data */
constexpr std::ios_base::openmode ASCII_IN = std::ios_base::in;
/** Format to write ascii data */
constexpr std::ios_base::openmode ASCII_OUT = std::ios_base::out;

/** Select which data are written in the file
    - python (default):
    row
    a0 a1 a2 ...
    - nodim:
    a0 a1 a2 ...

 */
enum class WriteType {
  python, /**< default: dimensions on the first line and then the content of the object */
  nodim,  /**< only matrix/vector values, no dimensions at the beginning of the file */
};

/** Read a SiconosVector from a file
 *  \param[in] fileName the file containing the vector
 *  \param[in,out] m the SiconosVector to be filled
 *  \param[in] mode ios_base::openmode, mode for reading (like  ios::in|ios:binary ...)
 *       default = ascii
 *  \param[in] precision value for float output. Default = 15.
 *  \param[in] with or without dimensions (see WriteType enum)
 *  \param[in] flags for reading
 *  \return bool true if read ok, else false ...
 */
bool read(const std::string &fileName, SiconosVector &m,
          const std::ios_base::openmode &mode = ASCII_IN, int precision = 15,
          const WriteType = WriteType::python,
          const std::ios::fmtflags &flags = std::cin.flags());

/** Write a SiconosVector to a file
    \param[in] fileName output file name
    \param[in] mode ios_base::openmode, mode for writing (like  ios::out|ios:binary ...)
    default = ascii
    \param[in] flags
    \param[in,out] m the SiconosVector to be written
    \param[in] precision value for float output. Default = 15.
    \param[in] outputType output format, choose between
        - python (default):
        row
        a0 a1 a2 ...
        - noDim:
        a0 a1 a2 ...
        Reading input format is the one corresponding to "python".
    \param[in] flags
    \return bool true if read ok, else false ...
*/
bool write(const std::string &fileName, const SiconosVector &m,
           const std::ios_base::openmode &mode = ASCII_OUT, int precision = 15,
           const WriteType = WriteType::python,
           const std::ios_base::fmtflags &flags = std::cout.flags());

/**
   Read a SiconosMatrix

   \param[in] fileName the name of the file to read
   \param[in,out] m the SiconosMatrix to be filled
   \param[in] mode the storage type used in the file (either ASCII_IN or BINARY_IN)
   \return true if read ok, else false ...
*/
bool read(const std::string &fileName, SiconosMatrix &m,
          const std::ios_base::openmode &mode = ASCII_IN);

/**
   Write a SiconosMatrix

   \param[in] fileName the name of the file to write in
   \param[in] m the SiconosMatrix to write
   \param[in] mode the storage type used in the file (either ASCII_OUT or BINARY_OUT)
   \param[in] outputType type of output:
   - python (default):
   row col
   a00 a01 a02 ...
   a10 ...
   - noDim:
   a00 a01 a02 ...
   a10 ...
   Reading input format is the one corresponding to "python".
   \return true if read ok, else false ...
*/
bool write(const std::string &fileName, const SiconosMatrix &m,
           const std::ios_base::openmode &mode = ASCII_OUT,
           const WriteType = WriteType::python);

/** Function to load data from a file and compare it with the provided
 *  data.  Returns the measured difference between files if the file
 *  was loaded and the comparison was performed, which must be >= 0.0,
 *  otherwise -1.0 is returned.  Caller needs to check diff <= epsilon
 *  to verify the result.
 *
 *  \param data The data to compare against the file.
 *  \param filename The name of the file to load and compare.
 *  \param epsilon The comparison threshold.
 *  \param index An optional list of column indexes, size==0 indicates all columns.
 *  \param[in] mode the storage type used in the file (either ASCII_IN or BINARY_IN)
 *  \param verbose True to print verbose output.
 *  \return Positive or 0.0 if the file was loaded and the comparison was performed,
 *  otherwise -1.
 */
double compareRefFile(const SimpleMatrix &data, std::string filename, double epsilon,
                      std::vector<int> index = {},
                      const std::ios_base::openmode mode = std::ios_base::in,
                      bool verbose = true);

/** \returns a pointer to a SiconosVector, built from json input
    \param jin json input
*/
std::shared_ptr<SiconosVector> readVectorFromJson(const nlohmann::json &jin);
}  // namespace io
}  // namespace siconos::algebra
#endif
