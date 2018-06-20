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

/*! \file SiconosMemory.hpp
    \brief class SiconosMemory

*/


#ifndef SICONOSMEMORY_H
#define SICONOSMEMORY_H

#include "SiconosMemoryException.hpp"
#include "SiconosPointers.hpp"
#include <deque>
#include "SiconosAlgebraTypeDef.hpp"

/** Container used to save vectors in SiconosMemory */
typedef std::vector<SiconosVector> MemoryContainer;
TYPEDEF_SPTR(MemoryContainer)

/** This class is a backup for vectors of previous time step

    There is a max number of saved vector (memorySize) and all the vector (simple or block)
    should have the same size.

*/
class SiconosMemory : public MemoryContainer
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(SiconosMemory);

  /** the real number of SiconosVectors saved in the Memory (ie the ones for which memory has been allocated) */
  MemoryContainer::size_type _nbVectorsInMemory;

  /** index to avoid removal and creation of vectors */
  MemoryContainer::size_type _indx;

public:

  /** default constructor. */
  SiconosMemory() : _nbVectorsInMemory(0), _indx(0) {};

  /** constructor with size parameter.
   * \param size size of the MemoryContainer
   * \param vectorSize the size of the SiconosVector to store
   */
  SiconosMemory(const unsigned int size, const unsigned int vectorSize);

  /** constructor with deque parameter.
   * \param deque MemoryContainer, the deque of siconosVector which must be stored
   * _size is set to the size of the deque given in parameters
   */
  SiconosMemory(const MemoryContainer& deque);

  /** constructor with size and deque parameter.
   * \param size int , the size of the memory
   * \param deque MemoryContainer, the deque of siconosVector which must be stored
   * this constructor is useful if the deque given in parameters has a size lower than the normal size of the memory
   */
  SiconosMemory(const unsigned int size, const MemoryContainer& deque);

  /** Copy constructor
   * \param Mem a SiconosMemory
   */
  SiconosMemory(const SiconosMemory& Mem);

  /** destructor */
  ~SiconosMemory();

  /** Assignment
   */
  void operator=(const SiconosMemory&);

  /** Assignment from container
   * Assumes all entries are valid SiconosVectors
   */
  void operator=(const MemoryContainer& V);

  /*************************************************************************/

  /** fill the memory with a vector of siconosVector
   * \param v MemoryContainer
   * \param size of the input container
   */
  void setVectorMemory(const MemoryContainer& v, MemoryContainer::size_type size);

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

  /** gives the size of the memory
   * \return int >= 0
   */
  inline unsigned int getMemorySize() const
  {
    return size();
  };

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

typedef std::vector<SiconosMemory> VectorOfMemories;

#endif

