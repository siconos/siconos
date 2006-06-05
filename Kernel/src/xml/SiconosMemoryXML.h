/* Siconos-Kernel version 1.2.0, Copyright INRIA 2005-2006.
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
#ifndef SICONOSMEMORYXML_H
#define SICONOSMEMORYXML_H

#include "SiconosDOMTreeTools.h"
#include<deque>

/** \class SiconosMemoryXML
 *   \brief This class manages SiconosMemory data part
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.2.0.
 *   \date 07/13/2004
 *
 *
 */

const std::string SM_MEMORYSIZE = "sizeMax";
const std::string SM_MEMORY = "Memory";

class SiconosMemoryXML
{
public:

  SiconosMemoryXML();
  SiconosMemoryXML(xmlNode* memoryNode, xmlNode* parentNode = NULL, const std::string& name = "default");
  ~SiconosMemoryXML();

  /** \fn deque<SiconosVector*> getSiconosMemoryVector()
   *  \brief allows to get the deque of SiconosVector from a SiconosMemory in the XML
   *  \return deque<SiconosVector*>
   */
  inline std::deque<SiconosVector*> getSiconosMemoryVector()
  {
    return getVectorMemoryValue();
  }

  /** \fn int getSiconosMemorySize()
   *  \brief allows to get the size max of the SiconosMemory
   *  \return int : the max size of the SiconosMemory
   */
  inline int getSiconosMemorySize()
  {
    return SiconosDOMTreeTools::getAttributeValue<int>(memoryNode, NUMBER_ATTRIBUTE);
  }

  /** \fn void setSiconosMemoryVector(deque<SiconosVector*> v)
   *  \brief allows to set the vector of SiconosVector of a SiconosMemory in the XML
   *  \param deque<SiconosVector*> to set
   */
  inline void setSiconosMemoryVector(const std::deque<SiconosVector*>& v)
  {
    setVectorMemoryValue(v);
  }

  /** \fn void setSiconosMemorySize(int s)
   *  \brief allows to set the value of the max size of the SiconosMemory
   *  \param int : the value to set
   */
  inline void setSiconosMemorySize(const unsigned int& s)
  {
    SiconosDOMTreeTools::setIntegerAttributeValue(memoryNode, SM_MEMORYSIZE, s);
  }

  /** \fn bool hasMemory()
   *  \brief determines if the SiconosMemoryXML contains memory objects
   *  \return bool : true if there's memories defined, false otherwise
   */
  inline bool hasMemory()
  {
    bool res = false;
    if (SiconosDOMTreeTools::findNodeChild(memoryNode, SM_MEMORY) != NULL) res = true;
    return res;
  }

  /** \fn xmlNode* getSiconosMemoryXMLNode()
   *  \brief returns the xmlNode of the SiconosMemoryXML
   *  \return xmlNode* : the value of the xmlNode of the SiconosMemoryXML
   */
  inline xmlNode* getSiconosMemoryXMLNode()
  {
    return memoryNode;
  }


  /** \fn void deleteUnusedMemoryNodes( int nbGoodNode )
   *  \brief deletes the nodes which won't be used (when there's more than maxSize = nbGoodNode SiconosVector in the SiconosMemory)
   */
  void deleteUnusedMemoryNodes(const int&);

private:

  /** \fn deque<SiconosVector*> getVectorMemoryValue()
   *   \brief Return a vector of SiconosVector computed from a memory node
   *   \param memoryNode : the memory node you want to get in a vector of SiconosVector type
   *   \return A  deque of SiconosVector
   */
  std::deque<SiconosVector*> getVectorMemoryValue();

  /** \fn void setVectorMemoryValue(deque<SiconosVector*>& memory)
   *   \brief Change values of a memoryNode from a deque<SiconosVector>
   *   \param memoryNode : the memory node you want to set
   *   \param memory : the memory you want to copy the value in the memoryNode
   *   \exception XMLException
   */
  void setVectorMemoryValue(const std::deque<SiconosVector*>& memory);

  xmlNode * memoryNode;
  xmlNode * parentNode;
};

#endif // SICONOSMEMORYXML_H
