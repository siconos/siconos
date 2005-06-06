#include "SiconosMemoryXML.h"
using namespace std;

SiconosMemoryXML::SiconosMemoryXML()
{}

SiconosMemoryXML::SiconosMemoryXML(xmlNode* memoryNode, xmlNode* parentNode, string name)
{
  /*
   * if parent == NULL, the parent node is not defined because the memoryNode exists
   * otherwise, we must create the SiconosMemoryNode with the parentNode
   */
  if (parentNode == NULL)
  {
    if (memoryNode != NULL)
    {
      this->memoryNode = memoryNode;
    }
    else
    {
      //XMLException::selfThrow("SiconosMemoryXML - constructor : element not found.");
      XMLException::selfThrow("SiconosMemoryXML - constructor : element '" + name + "' not found, memoryNode == NULL and parentNode == NULL.");
    }
    this->parentNode = NULL;
  }
  else if (name != "default")
  {
    /*
     * we create the node for the memory with no attributes
     */
    xmlNode * node;
    this->parentNode = parentNode;

    node = xmlNewChild(parentNode, NULL, BAD_CAST name.c_str(), NULL);
    xmlNewProp(node, (xmlChar*)(SM_MEMORYSIZE.c_str()), (xmlChar*)"0");
    this->memoryNode = node;
  }
  else XMLException::selfThrow("SiconosMemoryXML - constructor : illegal tag name.");
}

SiconosMemoryXML::~SiconosMemoryXML()
{}



vector<SiconosVector*> SiconosMemoryXML::getVectorMemoryValue()
{
  vector<SiconosVector*> v;
  xmlNode *node = SiconosDOMTreeTools::findNodeChild(this->memoryNode, SM_MEMORY);

  int cpt = 0;
  while (node != NULL)
  {
    SiconosVector *sv = new /*SiconosVector*/SimpleVector();

    *sv = SiconosDOMTreeTools::getSiconosVectorValue(node);

    v.push_back(sv);
    node = SiconosDOMTreeTools::findFollowNode(node, SM_MEMORY);
    cpt ++;
  }

  return v;
}

//SiconosMemory SiconosMemoryXML::getMemoryValue()
//{
//  vector<SiconosVector*> v = getVectorMemoryValue();
//  int size = SiconosDOMTreeTools::getIntegerAttributeValue( this->memoryNode, SM_MEMORY );
//  SiconosMemory sm( size, v );
//  return sm;
//}


void SiconosMemoryXML::setVectorMemoryValue(const vector<SiconosVector*> memory)
{
  IN("SiconosMemoryXML::setVectorMemoryValue\n");

  xmlNode *oldNode = SiconosDOMTreeTools::findNodeChild(this->memoryNode, SM_MEMORY);
  xmlNode *node; /* oldNode is the node before node */
  string stringValue;
  stringstream sstr;
  node = oldNode;

  int i = 0;

  /*
   * if oldNode == NULL now, that's because the memory is declared in the xml file with no vector
   * only the max size of the memory is given
   */
  if (oldNode != NULL)
  {
    while ((node != NULL) && (i < memory.size()))
    {
      SiconosDOMTreeTools::setSiconosVectorNodeValue(node, *(memory[i]));
      oldNode = node;
      node = SiconosDOMTreeTools::findFollowNode(node, SM_MEMORY);
      i++;
    }
  }

  while ((i < memory.size()) && (memory[i]->size() > 0)) //not enought nodes in the DOM tree to save memory
  {
    if (oldNode == NULL) node = xmlNewChild(this->memoryNode, NULL, BAD_CAST SM_MEMORY.c_str(), NULL);
    else node = xmlNewNode(NULL, BAD_CAST SM_MEMORY.c_str());

    stringValue = memory[i]->toString();
    xmlNodeSetContent(node, (const xmlChar *)stringValue.c_str());
    if (oldNode != NULL) oldNode->next = node;
    oldNode = node;
    i++;
  }

  OUT("SiconosMemoryXML::setVectorMemoryValue\n");
}

void SiconosMemoryXML::deleteUnusedMemoryNodes(int nbGoodNode)
{
  IN("SiconosMemoryXML::deleteUnusedMemoryNodes( int nbGoodNode )\n");

  xmlNode * node = SiconosDOMTreeTools::findNodeChild(this->memoryNode, SM_MEMORY);
  xmlNode * tmp;

  for (int i = 0; i < nbGoodNode; i++)
  {
    node = node->next;
  }

  /*
   * now, these are the nodes to delete
   */
  while (node != NULL)
  {
    tmp = node->next;
    xmlUnlinkNode(node);
    xmlFreeNode(node);

    node = tmp;
  }

  OUT("SiconosMemoryXML::deleteUnusedMemoryNodes( int nbGoodNode )\n");
}

//void SiconosMemoryXML::setMemoryValue( const SiconosMemory & memory)
//{
//  int i=0;
//  xmlNode *oldNode;
//  xmlNode *node, *parent;
//  string stringValue;
//
//  if( memoryNode != NULL )
//  {
//    oldNode = SiconosDOMTreeTools::findNodeChild(this->memoryNode, SM_MEMORY);
//    node=oldNode;
//
//    while ( (node!=NULL)&&(i<memory.getNbVectorsInMemory()) )
//    {
//      SiconosDOMTreeTools::setSiconosVectorValue( node, memory.getSiconosVector(i) );
//      oldNode = node;
//      node = SiconosDOMTreeTools::findFollowNode(node, SM_MEMORY);
//      i++;
//    }
//
//    while ( (i<memory.getMemorySize()) && (i<memory.getNbVectorsInMemory()))  //not enought nodes in the DOM tree to save memory
//    {
//      node = xmlNewNode(NULL, BAD_CAST SM_MEMORY.c_str());
//      stringValue = memory.getSiconosVector(i)->toString();
//      xmlNodeSetContent(node, (const xmlChar *)stringValue.c_str());
//      xmlAddChild( this->memoryNode, node );
////      oldNode->next = node;
////      oldNode = node;
//      i++;
//    }
//    SiconosDOMTreeTools::setIntegerAttributeValue( this->memoryNode, SM_MEMORYSIZE, memory.getMemorySize() );
//  }
//  else
//  {
//    XMLException::selfThrow("SiconosMemoryXML - setMemoryValue : memoryNode == NULL");
//  }
//
//
//}

