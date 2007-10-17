/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
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
#include "FirstOrderRXML.h"
using namespace std;

FirstOrderRXML::FirstOrderRXML(): RelationXML(), gNode(NULL), hNode(NULL)
{
  jacobianGNode.resize(1, NULL);
  jacobianHNode.resize(2, NULL);
}

FirstOrderRXML::FirstOrderRXML(xmlNodePtr relationNode): RelationXML(relationNode), gNode(NULL), hNode(NULL)
{
  xmlNodePtr node;
  // g function
  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, "g")) != NULL)
    gNode = node;

  // h function
  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, "h")) != NULL)
    hNode = node;

  // Gradients ...
  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, "jacobianG")) != NULL)
  {
    // get number of functions given
    unsigned int size = SiconosDOMTreeTools::getAttributeValue<unsigned int>(node, "number");
    if (size == 0)
      XMLException::selfThrow("FirstOrderRXML:: constructor failed. Some Gradients plug-in names are required (jacobianG).");
    jacobianGNode.resize(size, NULL);
    // get corresponding nodes
    jacobianGNode[0] = SiconosDOMTreeTools::findNodeChild(node, "matrix")  ;
    if (jacobianGNode[0] == NULL)
      XMLException::selfThrow("FirstOrderRXML:: constructor, jacobianG0 is missing");

    for (unsigned int i = 1; i < size; i++)
    {
      jacobianGNode[i] =  SiconosDOMTreeTools::findFollowNode(jacobianGNode[i - 1]);
      if (jacobianGNode[i] == NULL)
        XMLException::selfThrow("FirstOrderRXML:: constructor, another gradient of G is required");
    }
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, "jacobianH")) != NULL)
  {
    // get number of functions given
    unsigned int size = SiconosDOMTreeTools::getAttributeValue<unsigned int>(node, "number");
    if (size == 0)
      XMLException::selfThrow("FirstOrderRXML:: constructor failed. Some Gradients plug-in names are required (jacobianH).");
    jacobianHNode.resize(size, NULL);
    // get corresponding nodes
    jacobianHNode[0] = SiconosDOMTreeTools::findNodeChild(node, "matrix")  ;
    if (jacobianHNode[0] == NULL)
      XMLException::selfThrow("FirstOrderRXML:: constructor, jacobianH0 is missing");

    for (unsigned int i = 1; i < size; i++)
    {
      jacobianHNode[i] =  SiconosDOMTreeTools::findFollowNode(jacobianHNode[i - 1]);
      if (jacobianHNode[i] == NULL)
        XMLException::selfThrow("FirstOrderRXML:: constructor, another gradient of H is required");
    }

  }
}

FirstOrderRXML::~FirstOrderRXML()
{}

// ================== g ==================

void FirstOrderRXML::setGPlugin(const string&  plugin)
{
  if (gNode == NULL)
  {
    gNode = SiconosDOMTreeTools::createSingleNode(rootNode, "computeG");
    xmlNewProp(gNode, (xmlChar*)"plugin", (xmlChar*)plugin.c_str());
  }
  else SiconosDOMTreeTools::setStringAttributeValue(gNode, "plugin", plugin);
}

string FirstOrderRXML::getGPlugin() const
{
  if (!isGPlugin())
    XMLException::selfThrow("FirstOrderRXML - getComputeGPlugin : g is not calculated from a plugin");
  return  SiconosDOMTreeTools::getStringAttributeValue(gNode, "plugin");
}

// ================== h ==================

void FirstOrderRXML::setHPlugin(const string&  plugin)
{
  if (hNode == NULL)
  {
    hNode = SiconosDOMTreeTools::createSingleNode(rootNode, "computeH");
    xmlNewProp(hNode, (xmlChar*)"plugin", (xmlChar*)plugin.c_str());
  }
  else SiconosDOMTreeTools::setStringAttributeValue(hNode, "plugin", plugin);
}

string FirstOrderRXML::getHPlugin() const
{
  if (!isHPlugin())
    XMLException::selfThrow("FirstOrderRXML - getComputeHPlugin : h is not calculated from a plugin");
  return  SiconosDOMTreeTools::getStringAttributeValue(hNode, "plugin");
}

// ================== jacobianG ==================

void FirstOrderRXML::setJacobianGPlugin(const std::string&, unsigned int)
{
  XMLException::selfThrow("FirstOrderRXML -  setJacobianGPlugin: not yet implemented.");
}

string FirstOrderRXML::getJacobianGPlugin(unsigned int index) const
{
  if (index >= jacobianGNode.size())
    XMLException::selfThrow("FirstOrderRXML - getJacobianGPlugin(index), index out of range");

  if (!isJacobianGPlugin(index))
    XMLException::selfThrow("FirstOrderRXML -  getJacobianGPlugin: not computed with a plugin");

  return SiconosDOMTreeTools::getStringAttributeValue(jacobianGNode[index], "matrixPlugin");
}

SimpleMatrix FirstOrderRXML::getJacobianGMatrix(unsigned int index) const
{
  if (index >= jacobianGNode.size())
    XMLException::selfThrow("FirstOrderRXML - getJacobianGMatrix(index), index out of range");
  if (isJacobianGPlugin(index))
    XMLException::selfThrow("FirstOrderRXML - getJacobianGMatrix: JacobianG is computed using a plug-in");
  return  SiconosDOMTreeTools::getSiconosMatrixValue(jacobianGNode[index]);
}

// ================== jacobianH ==================

void FirstOrderRXML::setJacobianHPlugin(const std::string&, unsigned int)
{
  XMLException::selfThrow("FirstOrderRXML -  setJacobianHPlugin: not yet implemented.");
}

string FirstOrderRXML::getJacobianHPlugin(unsigned int index) const
{
  if (index >= jacobianHNode.size())
    XMLException::selfThrow("FirstOrderRXML - getJacobianHPlugin(index), index out of range");

  if (!isJacobianHPlugin(index))
    XMLException::selfThrow("FirstOrderRXML -  getJacobianHPlugin: not computed with a plugin");

  return SiconosDOMTreeTools::getStringAttributeValue(jacobianHNode[index], "matrixPlugin");
}

SimpleMatrix FirstOrderRXML::getJacobianHMatrix(unsigned int index) const
{
  if (index >= jacobianHNode.size())
    XMLException::selfThrow("FirstOrderRXML - getJacobianHMatrix(index), index out of range");
  if (isJacobianHPlugin(index))
    XMLException::selfThrow("FirstOrderRXML - getJacobianHMatrix: JacobianH is computed using a plug-in");
  return  SiconosDOMTreeTools::getSiconosMatrixValue(jacobianHNode[index]);
}


