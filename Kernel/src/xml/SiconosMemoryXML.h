
#ifndef SICONOSMEMORYXML_H
#define SICONOSMEMORYXML_H

//#include "SiconosMemory.h"
#include <string>
#include <libxml/tree.h>
#include "SiconosDOMTreeTools.h"

#include "XMLException.h"


//using namespace std;

/** \class SiconosMemoryXML
*   \brief This class manages SiconosMemory data part
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date 07/13/2004
*
*
*/


const string SM_MEMORYSIZE = "sizeMax";
const string SM_MEMORY = "Memory";

class SiconosMemoryXML
{
public:

  SiconosMemoryXML();
  SiconosMemoryXML(xmlNode* memoryNode, xmlNode* parentNode = NULL, string name = "default");
  ~SiconosMemoryXML();

  /** \fn vector<SiconosVector*> getSiconosMemoryVector()
   *  \brief allows to get the vector of SiconosVector from a SiconosMemory in the XML
   *  \return vector<SiconosVector*>
   */
  inline vector<SiconosVector*> getSiconosMemoryVector()
  {
    return getVectorMemoryValue();
  }

  /** \fn int getSiconosMemorySize()
   *  \brief allows to get the size max of the SiconosMemory
   *  \return int : the max size of the SiconosMemory
   */
  inline int getSiconosMemorySize()
  {
    return SiconosDOMTreeTools::getIntegerAttributeValue(this->memoryNode, SM_MEMORYSIZE);
  }

  /** \fn void setSiconosMemoryVector(vector<SiconosVector*> v)
   *  \brief allows to set the vector of SiconosVector of a SiconosMemory in the XML
   *  \param vector<SiconosVector*> to set
   */
  inline void setSiconosMemoryVector(vector<SiconosVector*> v)
  {
    this->setVectorMemoryValue(v);
  }

  /** \fn void setSiconosMemorySize(int s)
   *  \brief allows to set the value of the max size of the SiconosMemory
   *  \param int : the value to set
   */
  inline void setSiconosMemorySize(int s)
  {
    SiconosDOMTreeTools::setIntegerAttributeValue(this->memoryNode, SM_MEMORYSIZE, s);
  }

  /** \fn bool hasMemory()
   *  \brief determines if the SiconosMemoryXML contains memory objects
   *  \return bool : true if there's memories defined, false otherwise
   */
  inline bool hasMemory()
  {
    bool res = false;
    if (SiconosDOMTreeTools::findNodeChild(this->memoryNode, SM_MEMORY) != NULL) res = true;
    return res;
  }

  /** \fn xmlNode* getSiconosMemoryXMLNode()
   *  \brief returns the xmlNode of the SiconosMemoryXML
   *  \return xmlNode* : the value of the xmlNode of the SiconosMemoryXML
   */
  inline xmlNode* getSiconosMemoryXMLNode()
  {
    return this->memoryNode;
  }


  /** \fn void deleteUnusedMemoryNodes( int nbGoodNode )
   *  \brief deletes the nodes which won't be used (when there's more than maxSize = nbGoodNode SiconosVector in the SiconosMemory)
   */
  void deleteUnusedMemoryNodes(int nbGoodNode);

private:

  /** \fn vector<SiconosVector*> getVectorMemoryValue()
  *   \brief Return a vector of SiconosVector computed from a memory node
  *   \param memoryNode : the memory node you want to get in a vector of SiconosVector type
  *   \return A  vector of SiconosVector
  */
  vector<SiconosVector*> getVectorMemoryValue();

  /** \fn void setVectorMemoryValue(vector<SiconosVector*> memory)
  *   \brief Change values of a memoryNode from a vector<SiconosVector>
  *   \param memoryNode : the memory node you want to set
  *   \param memory : the memory you want to copy the value in the memoryNode
  *   \exception XMLException
  */
  void setVectorMemoryValue(const vector<SiconosVector*> memory);

  xmlNode * memoryNode;
  xmlNode * parentNode;
};

#endif // SICONOSMEMORYXML_H
