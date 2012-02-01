/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
/*! \file SiconosMemoryXML.hpp

*/
#ifndef SICONOSMEMORYXML_H
#define SICONOSMEMORYXML_H

#include "SiconosPointers.hpp"
#include "SiconosDOMTreeTools.hpp"
#include "SiconosMemory.hpp"

/** Maximum size for the memory */
const std::string SM_MEMORYSIZE = "sizeMax";

const std::string SM_MEMORY = "Memory";

/** XML management for SiconosMemory
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 3.0.0.
 *   \date 07/13/2004
 *
 *
 */
class SiconosMemoryXML
{
  xmlNodePtr _memoryNode;
  xmlNodePtr _parentNode;

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(SiconosMemoryXML);


  /** Return a vector of SiconosVector computed from a memory node
  *   \return A  deque of SiconosVector
  */
  SP::MemoryContainer getVectorMemoryValue();

  /** Change values of a _memoryNode from a deque<SiconosVector>
  *   \param memory : the memory you want to copy the value in the _memoryNode
  *   \exception XMLException
  */
  void setVectorMemoryValue(const MemoryContainer& memory);

  /** Default constructor */
  SiconosMemoryXML();

public:

  /** Basic constructor
   * \param newMemoryNode a xmlNodePtr
   * \param newParentNode a xmlNodePtr (optional)
   * \param name a string (optional)
   */
  SiconosMemoryXML(xmlNodePtr newMemoryNode, xmlNodePtr newParentNode = NULL, const std::string& name = "default");

  /* Copy constructor
   * \param MemXML a SiconosMemoryXML to copy
   */
  SiconosMemoryXML(const SiconosMemoryXML & MemXML);

  /** Destructor */
  ~SiconosMemoryXML() {};

  /** allows to get the deque of SiconosVector from a SiconosMemory in the XML
  *  \return deque<SiconosVector*>
  */
  inline SP::MemoryContainer getSiconosMemoryVector()
  {
    return getVectorMemoryValue();
  }

  /** allows to get the size max of the SiconosMemory
  *  \return int : the max size of the SiconosMemory
  */
  inline int getSiconosMemorySize()
  {
    return SiconosDOMTreeTools::getAttributeValue<int>(_memoryNode, NUMBER_ATTRIBUTE);
  }

  /** allows to set the vector of SiconosVector of a SiconosMemory in the XML
  *  \param MemoryContainer to set
  */
  inline void setSiconosMemoryVector(const MemoryContainer& v)
  {
    setVectorMemoryValue(v);
  }

  /** allows to set the value of the max size of the SiconosMemory
  *  \param int : the value to set
  */
  inline void setSiconosMemorySize(const unsigned int& s)
  {
    SiconosDOMTreeTools::setIntegerAttributeValue(_memoryNode, SM_MEMORYSIZE, s);
  }

  /** determines if the SiconosMemoryXML contains memory objects
  *  \return bool : true if there's memories defined, false otherwise
  */
  inline bool hasMemory()
  {
    return SiconosDOMTreeTools::findNodeChild(_memoryNode, SM_MEMORY);
  }

  /** returns the xmlNode of the SiconosMemoryXML
  *  \return _memoryNode
  */
  inline xmlNodePtr getSiconosMemoryXMLNode() const
  {
    return _memoryNode;
  }

  /** returns the xmlNode of the SiconosMemoryXML
  *  \return _parentNode
  */
  inline xmlNodePtr getSiconosParentXMLNode() const
  {
    return _parentNode;
  }

  /** deletes the nodes which won't be used (when there's more than maxSize = nbGoodNode SiconosVector in the SiconosMemory)
   */
  void deleteUnusedMemoryNodes(const int&);

};
#endif // SICONOSMEMORYXML_H
