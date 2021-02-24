/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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

/*! \file SiconosMemory.hpp
  \brief class SiconosMemory

*/


#ifndef SICONOSMEMORY_H
#define SICONOSMEMORY_H

#include "SiconosPointers.hpp"
#include <deque>
#include "SiconosAlgebraTypeDef.hpp"
#include "SiconosSerialization.hpp" // For ACCEPT_SERIALIZATION

/** Container used to save vectors in SiconosMemory */
typedef std::vector<SiconosVector> MemoryContainer;

/** a set of memory vectors */
typedef std::vector<SiconosMemory> VectorOfMemories;

TYPEDEF_SPTR(MemoryContainer)

/** Interface to stl container of SiconosVector.

    This class is used as a backup during simulation, to save vectors (e.g. state) computed during previous time steps.

    - The size of the container is fixed, with a first-in first-out mechanism
    used through swap method.
    - All saved vectors must have the same dimension.

    This class must be reviewed and backup should probably be moved to graph rather than in this object.

*/
class SiconosMemory : public MemoryContainer
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(SiconosMemory);

  /** the real number of SiconosVectors saved in the Memory (ie the ones for which memory has been allocated) */
  MemoryContainer::size_type _nbVectorsInMemory = 0;

  /** index to avoid removal and creation of vectors.
   this[_indx] is to the oldest element in the set */
  MemoryContainer::size_type _indx = 0;


  /* Forbid  assignment */
  //SiconosMemory(const SiconosMemory& Mem) = delete;
  //void operator=(const SiconosMemory&) = delete;
  SiconosMemory(const MemoryContainer&) = delete;
  void operator=(const MemoryContainer& V);

public:

  /** creates an empty SiconosMemory. */
  SiconosMemory() = default;

  /** creates a SiconosMemory
   * \param size number of elements in the container
   * \param vectorSize size of each vector in the container
   */
  SiconosMemory(const unsigned int size, const unsigned int vectorSize);

  /** creates a SiconosMemory, copy constructor
   * Required because of resize call in DS initMemory function.
   * \param mem a SiconosMemory
   */
  SiconosMemory(const SiconosMemory& mem);


  /** destructor */
  ~SiconosMemory(){};

  /** Assignment
   */
  void operator=(const SiconosMemory&);

  /** To get SiconosVector number i of the memory
   * \param int i: the position in the memory of the wanted SiconosVector
   * \return a SP::SiconosVector
   */
  const SiconosVector& getSiconosVector(const unsigned int) const;

  /** To get SiconosVector number i of the memory as mutable reference.
   * Use should be avoided whenever possible.
   * (Used in LinearSMC::actuate)
   * \param int i: the position in the memory of the wanted SiconosVector
   * \return a SP::SiconosVector
   */
  SiconosVector& getSiconosVectorMutable(const unsigned int);

  /** set size of the SiconosMemory (number of vectors and size of vector)
   * \param steps the max size for this SiconosMemory, size of the container
   * \param vectorSize size of each vector of the container
   */
  void setMemorySize(const unsigned int steps,
                     const unsigned int vectorSize);

  /** gives the numbers of SiconosVectors currently stored in the memory
   * \return int >= 0
   */
  inline unsigned int nbVectorsInMemory() const
  {
    return _nbVectorsInMemory;
  };

  /** puts a SiconosVector into the memory
   * \param v the SiconosVector we want to put in memory
   */
  void swap(const SiconosVector& v);

  /** puts a SiconosVector into the memory
   * \param v the SiconosVector we want to put in memory, or do nothing if v is null
   */
  void swap(SP::SiconosVector v);

  /** displays the data of the memory object
   */
  void display() const;

};


#endif

