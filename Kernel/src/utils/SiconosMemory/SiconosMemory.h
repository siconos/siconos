/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/

/*! \file SiconosMemory.h
    \brief class SiconosMemory

*/


#ifndef SICONOSMEMORY_H
#define SICONOSMEMORY_H

#include "SiconosMemoryXML.h"
#include "SiconosMemoryException.h"
#include <deque>

/** This class is used to save vectors for several steps
 *
 *  There is a max number of saved vector (memorySize) and all the vector (simple or block)
 * should have the same size.
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 3.0.0.
 *   \date (Creation) 07/06/2004
 */
class SiconosMemory
{

public:

  /**
   * constructor with size parameter.
   * \param int : the size of the memory
   * memorySize is set with the parameter, and the memory is allocated for this number of SiconosVector
   */
  SiconosMemory(const unsigned int);

  /**
   * xml constructor + (optional) size of memory; if not set, size is read in xml file.
   * \param SiconosMemoryXML * : the XML object which contains the data of the memory
   * \param unsigned int : the size of the memory (OPTIONAL)
   */
  SiconosMemory(SiconosMemoryXML *, const unsigned int = 1);

  /**
   * constructor with deque parameter.
   * \param deque<SiconosVector*> : the deque of siconosVector which must be stored
   * memorySize is set to the size of the deque given in parameters
   */
  SiconosMemory(const std::deque<SiconosVector*>&);

  /**
   * constructor with size and deque parameter.
   * \param int : the size of the memory
   * \param deque<SiconosVector*> : the deque of siconosVector which must be stored
   * this constructor is useful if the deque given in parameters has a size lower than the normal size of the memory
   */
  SiconosMemory(const unsigned int, const std::deque<SiconosVector*>&);

  /**
   * constructor by copy
   */
  SiconosMemory(const SiconosMemory&);

  /**
   * destructor
   * delete the SiconosVectors allocated in vectorMemory if memorySize > 0
   */
  ~SiconosMemory();

  /*************************************************************************/

  /**
   * fill the memory with a vector of siconosVector
   * \param deque<SiconosVector*>
   * memorySize is set to the size of the deque given in parameters
   */
  void setVectorMemory(const std::deque<SiconosVector*>&);

  /**
   * To get SiconosVector number i of the memory
   * \param int i: the position in the memory of the wanted SiconosVector
   * \return a SiconosVector*
   */
  SiconosVector* getSiconosVector(const unsigned int) const;

  /**
   * gives the size of the memory
   * \return int >= 0
   */
  inline const unsigned int getMemorySize() const
  {
    return memorySize;
  };

  /**
   * set the max size of the SiconosMemory
   * \param int : the max size for this SiconosMemory
   */
  inline void setSiconosMemorySize(const unsigned int max)
  {
    memorySize = max;
  };

  /**
   * gives the numbers of SiconosVectors currently stored in the memory
   * \return int >= 0
   */
  inline const unsigned int getNbVectorsInMemory() const
  {
    return nbVectorsInMemory;
  };

  /**
   * gives the vector of SiconosVectors of the memory
   * \return stl deque od siconosVector
   */
  inline std::deque<SiconosVector*> getVectorMemory() const
  {
    return vectorMemory;
  };

  /** Allows to get the SiconosMemoryXML of the SiconosMemory
   *  \return SiconosMemoryXML* : the object SiconosMemoryXML of the SiconosMemory
   */
  inline SiconosMemoryXML* getSiconosMemoryXML()
  {
    return memoryXML;
  }

  /**
   * puts a SiconosVector into the memory
   * \param SiconosVector* : the SiconosVector we want to put in memory
   */
  void swap(SiconosVector *);

  /**
   * displays the data of the memory object
   */
  void display() const;

  /** Copy the max size of the SiconosMemory to the XML DOM tree
   *  \exception SiconosMemoryException
   */
  inline void saveMemorySizeToXML()
  {
    if (memoryXML != NULL) memoryXML->setSiconosMemorySize(memorySize);
    else SiconosMemoryException::selfThrow("SiconosMemory::saveMemorySizeToXML() - memoryXML object == NULL");
  }

  SiconosMemory& operator = (const SiconosMemory&);


private:
  /**
  * default constructor.
  */
  SiconosMemory();

  /** the maximum size of the memory (i.e the max numbers of SiconosVectors it can store) */
  unsigned int memorySize;

  /** the real number of SiconosVectors saved in the Memory (ie the ones for which memory has been allocated) */
  unsigned int nbVectorsInMemory;

  /** the stl deque which contains the SiconosVectors kept in memory */
  /** we use deque rather than vector to access to push_front function */
  std::deque<SiconosVector*> vectorMemory;

  /** stl vector to provides info about memory allocation */
  std::deque<bool> isVectorMemoryAllocated;

  /** link to the XML for SiconosMemory objects */
  SiconosMemoryXML * memoryXML;

};

#endif

