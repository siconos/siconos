
#include "SiconosMemoryXML.h"

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
      SiconosDOMTreeTools::setSiconosVectorValue(node, *(memory[i]));
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

//$Log: SiconosMemoryXML.cpp,v $
//Revision 1.7  2004/09/10 11:26:29  charlety
//
//_ Integration of the new version of the SiconosVector in the platform. the class simpleVector is used mostly to replace old SiconosVector. When a vector can be composite or simple, like the state of a dynamical system, a pointer on SiconosVector is used, and the vector is initialized simple or composite when the system is initialized.
//
//_ All the tests which worked with the previous version of the vector are OK with the new version.
//
//_ Example SICONOS and bouncingBall are OK
//
//_ some comments have still to be adapted to NewSiconosVector .
//
//_ This version of NewSiconosVector could be called 0.9. some details have to be fixed, it will be done before the end of September.
//
//Revision 1.6  2004/08/04 11:03:23  jbarbier
//- about the SiconosMemory : when a SiconosMemory has a maxSize greater than the
//number of steps in memory required by an integrator, the oldest SiconosVector
//are deleted
//
//- the way to initialize the SiconosMemory by the integrator has been updated to
//match with these changes
//
//Revision 1.5  2004/08/03 12:07:12  jbarbier
//- all test on th eModel are successfull
//
//- new tests on the Model with the opening of XML file
//
//- link TimeDiscretisation -> Strategy
//
//- attribute T of the Model is now optional
//
//Revision 1.4  2004/08/02 09:26:26  jbarbier
//- xml save for SiconosMemory corrected
//- temporary operation in Moreau::integrate because of the current version of
//SiconosVector
//
//Revision 1.3  2004/07/30 14:37:16  jbarbier
//- saving methods for DynamicalSystemXML and LagrangianNLDSXML
//
//Revision 1.2  2004/07/29 14:25:45  jbarbier
