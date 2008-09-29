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
#include "SiconosMemoryXML.h"
#include "SimpleVector.h"

using namespace std;

SiconosMemoryXML::SiconosMemoryXML(): memoryNode(NULL), parentNode(NULL)
{}

SiconosMemoryXML::SiconosMemoryXML(xmlNode* newMemoryNode, xmlNode* newParentNode, const string& name):
  memoryNode(newMemoryNode), parentNode(newParentNode)
{
  /*
   * if parent == NULL, the parent node is not defined because the memoryNode exists
   * otherwise, we must create the SiconosMemoryNode with the parentNode
   */
  if (parentNode == NULL)
  {
    if (memoryNode == NULL)
      XMLException::selfThrow("SiconosMemoryXML - constructor : element '" + name + "' not found, memoryNode == NULL and parentNode == NULL.");
  }
  else if (name != "default")
  {
    // we create the node for the memory with no attributes
    xmlNode * node;
    node = xmlNewChild(parentNode, NULL, BAD_CAST name.c_str(), NULL);
    xmlNewProp(node, (xmlChar*)(SM_MEMORYSIZE.c_str()), (xmlChar*)"0");
    memoryNode = node;
  }
  else XMLException::selfThrow("SiconosMemoryXML - constructor : illegal tag name.");
}

SiconosMemoryXML::~SiconosMemoryXML()
{}

deque<SP::SiconosVector> SiconosMemoryXML::getVectorMemoryValue()
{
  deque<SP::SiconosVector> v;
  xmlNode *node = SiconosDOMTreeTools::findNodeChild(memoryNode, SM_MEMORY);

  int cpt = 0;
  while (node != NULL)
  {
    v.push_back(SP::SimpleVector(new SimpleVector(SiconosDOMTreeTools::getSiconosVectorValue(node))));
    node = SiconosDOMTreeTools::findFollowNode(node, SM_MEMORY);
    cpt ++;
  }
  return v;
}

void SiconosMemoryXML::setVectorMemoryValue(const deque<SP::SiconosVector>& memory)
{
  xmlNode *oldNode = SiconosDOMTreeTools::findNodeChild(memoryNode, SM_MEMORY);
  xmlNode *node; /* oldNode is the node before node */
  string stringValue;
  stringstream sstr;
  node = oldNode;

  unsigned int i = 0;

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
    if (oldNode == NULL) node = xmlNewChild(memoryNode, NULL, BAD_CAST SM_MEMORY.c_str(), NULL);
    else node = xmlNewNode(NULL, BAD_CAST SM_MEMORY.c_str());

    stringValue = memory[i]->toString();
    xmlNodeSetContent(node, (const xmlChar *)stringValue.c_str());
    if (oldNode != NULL) oldNode->next = node;
    oldNode = node;
    i++;
  }
}

void SiconosMemoryXML::deleteUnusedMemoryNodes(const int& nbGoodNode)
{
  xmlNode * node = SiconosDOMTreeTools::findNodeChild(memoryNode, SM_MEMORY);
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
}
