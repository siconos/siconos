/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
#include "SiconosMemoryXML.hpp"
#include "SiconosVector.hpp"

using namespace std;

SiconosMemoryXML::SiconosMemoryXML(xmlNodePtr newMemoryNode, xmlNodePtr newParentNode, const string& name):
  _memoryNode(newMemoryNode), _parentNode(newParentNode)
{
  /*
   * if parent == NULL, the parent node is not defined because the _memoryNode exists
   * otherwise, we must create the SiconosMemoryNode with the _parentNode
   */
  if (!_parentNode)
  {
    if (!_memoryNode)
      XMLException::selfThrow("SiconosMemoryXML - constructor : element '" + name + "' not found, memoryNode == NULL and parentNode == NULL.");
  }
  else if (name != "default")
  {
    // we create the node for the memory with no attributes
    xmlNodePtr node;
    node = xmlNewChild(_parentNode, NULL, BAD_CAST name.c_str(), NULL);
    xmlNewProp(node, (xmlChar*)(SM_MEMORYSIZE.c_str()), (xmlChar*)"0");
    _memoryNode = node;
  }
  else XMLException::selfThrow("SiconosMemoryXML - constructor : illegal tag name.");
}

// Copy constructor
SiconosMemoryXML::SiconosMemoryXML(const SiconosMemoryXML & MemXML):
  _memoryNode(MemXML.getSiconosMemoryXMLNode()), _parentNode(MemXML.getSiconosParentXMLNode())
{
}


SP::MemoryContainer SiconosMemoryXML::getVectorMemoryValue()
{
  SP::MemoryContainer v;
  xmlNodePtr node = SiconosDOMTreeTools::findNodeChild(_memoryNode, SM_MEMORY);

  int cpt = 0;
  while (node)
  {
    v->push_back(SP::SiconosVector(new SiconosVector(SiconosDOMTreeTools::getSiconosVectorValue(node))));
    node = SiconosDOMTreeTools::findFollowNode(node, SM_MEMORY);
    cpt ++;
  }
  return v;
}

void SiconosMemoryXML::setVectorMemoryValue(const MemoryContainer& memory)
{
  xmlNodePtr oldNode = SiconosDOMTreeTools::findNodeChild(_memoryNode, SM_MEMORY);
  xmlNodePtr node; /* oldNode is the node before node */
  string stringValue;
  stringstream sstr;
  node = oldNode;

  unsigned int i = 0;

  /*
   * if oldNode == NULL now, that's because the memory is declared in the xml file with no vector
   * only the max size of the memory is given
   */
  if (oldNode)
  {
    while ((node) && (i < memory.size()))
    {
      SiconosDOMTreeTools::setSiconosVectorNodeValue(node, *(memory[i]));
      oldNode = node;
      node = SiconosDOMTreeTools::findFollowNode(node, SM_MEMORY);
      i++;
    }
  }

  while ((i < memory.size()) && (memory[i]->size() > 0)) //not enought nodes in the DOM tree to save memory
  {
    if (! oldNode) node = xmlNewChild(_memoryNode, NULL, BAD_CAST SM_MEMORY.c_str(), NULL);
    else node = xmlNewNode(NULL, BAD_CAST SM_MEMORY.c_str());

    stringValue = memory[i]->toString();
    xmlNodeSetContent(node, (const xmlChar *)stringValue.c_str());
    if (oldNode) oldNode->next = node;
    oldNode = node;
    i++;
  }
}

void SiconosMemoryXML::deleteUnusedMemoryNodes(const int& nbGoodNode)
{
  xmlNodePtr  node = SiconosDOMTreeTools::findNodeChild(_memoryNode, SM_MEMORY);
  xmlNodePtr  tmp;

  for (int i = 0; i < nbGoodNode; i++)
  {
    node = node->next;
  }

  /*
   * now, these are the nodes to delete
   */
  while (node)
  {
    tmp = node->next;
    xmlUnlinkNode(node);
    xmlFreeNode(node);

    node = tmp;
  }
}
