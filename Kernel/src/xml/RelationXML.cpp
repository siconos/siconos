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
#include "RelationXML.hpp"
#include "SimpleMatrix.hpp"
using namespace std;
using namespace RELATION;

RelationXML::RelationXML(xmlNodePtr relationNode): rootNode(relationNode), hNode(NULL), gNode(NULL), hDotNode(NULL)
{
  xmlNodePtr node;
  // g function
  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, "g")))
    gNode = node;

  // h function
  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, "h")))
    hNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, "hDot")))
    hDotNode = node;

  // Jacobian ...
  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, "jacobianH")))
  {
    // get number of functions given
    unsigned int size = SiconosDOMTreeTools::getAttributeValue<unsigned int>(node, "number");
    if (size == 0)
      XMLException::selfThrow("RelationXML:: constructor failed. Some jacobian plug-in names are required (jacobianH).");
    jacobianHNode.resize(size, NULL);
    // get corresponding nodes
    jacobianHNode[0] = SiconosDOMTreeTools::findNodeChild(node, "matrix")  ;
    if (!jacobianHNode[0])
      XMLException::selfThrow("RelationXML:: constructor, jacobianH0 is missing");

    for (unsigned int i = 1; i < size; i++)
    {
      jacobianHNode[i] =  SiconosDOMTreeTools::findFollowNode(jacobianHNode[i - 1]);
      if (!jacobianHNode[i])
        XMLException::selfThrow("RelationXML:: constructor, another gradient of H is required");
    }

  }

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, "jacobianG")))
  {
    // get number of functions given
    unsigned int size = SiconosDOMTreeTools::getAttributeValue<unsigned int>(node, "number");
    if (size == 0)
      XMLException::selfThrow("RelationXML:: constructor failed. Some jacobian plug-in names are required (jacobianG).");
    jacobianGNode.resize(size, NULL);
    // get corresponding nodes
    jacobianGNode[0] = SiconosDOMTreeTools::findNodeChild(node, "matrix")  ;
    if (!jacobianGNode[0])
      XMLException::selfThrow("RelationXML:: constructor, jacobianG0 is missing");

    for (unsigned int i = 1; i < size; i++)
    {
      jacobianGNode[i] =  SiconosDOMTreeTools::findFollowNode(jacobianGNode[i - 1]);
      if (!jacobianGNode[i])
        XMLException::selfThrow("RelationXML:: constructor, another gradient of G is required");
    }
  }
}

const RELATION::TYPES RelationXML::getType() const
{
  std::string type((char*)rootNode->name);
  if (type == "LagrangianRelation")
    return Lagrangian;
  else if (type == "FirstOrderRelation")
    return FirstOrder;
  else
  {
    XMLException::selfThrow("RelationXML - getType: unknown type of Relation.");
    return Lagrangian;
  }
}

const RELATION::SUBTYPES RelationXML::getSubType() const
{
  std::string res = SiconosDOMTreeTools::getStringAttributeValue(rootNode, "type");
  if (res == "NonLinear")
    return NonLinearR;
  else if (res == "Linear")
    return LinearR;
  else if (res == "Type1")
    return Type1R;
  else if (res == "LinearTI")
    return LinearTIR;
  else if (res == "Scleronomous")
    return ScleronomousR;
  else if (res == "Rheonomous")
    return RheonomousR;
  else if (res == "Compliant")
    return CompliantR;
  else
  {
    XMLException::selfThrow("RelationXML - getType: unknown type of Relation.");
    return NonLinearR;
  }
}

unsigned int RelationXML::getNumberOfJacobians(unsigned int index) const
{
  string input;
  if (index == 0)
    input = "jacobianH";
  else
    input = "jacobianG";

  xmlNodePtr node = SiconosDOMTreeTools::findNodeChild(rootNode, input);
  if (!node)
    return 0;
  else
    return SiconosDOMTreeTools::getAttributeValue<unsigned int>(node, "number");
}

void RelationXML::setHPlugin(const string&  plugin)
{
  if (!hNode)
  {
    hNode = SiconosDOMTreeTools::createSingleNode(rootNode, "computeh");
    xmlNewProp(hNode, (xmlChar*)"plugin", (xmlChar*)plugin.c_str());
  }
  else SiconosDOMTreeTools::setStringAttributeValue(hNode, "plugin", plugin);
}

string RelationXML::getHPlugin() const
{
  if (!isHPlugin())
    XMLException::selfThrow("RelationXML - getComputeHPlugin : h is not calculated from a plugin");
  return  SiconosDOMTreeTools::getStringAttributeValue(hNode, "plugin");
}

void RelationXML::setGPlugin(const string&  plugin)
{
  if (!gNode)
  {
    gNode = SiconosDOMTreeTools::createSingleNode(rootNode, "computeG");
    xmlNewProp(gNode, (xmlChar*)"plugin", (xmlChar*)plugin.c_str());
  }
  else SiconosDOMTreeTools::setStringAttributeValue(gNode, "plugin", plugin);
}

string RelationXML::getGPlugin() const
{
  if (!isGPlugin())
    XMLException::selfThrow("RelationXML - getComputeGPlugin : g is not calculated from a plugin");
  return  SiconosDOMTreeTools::getStringAttributeValue(gNode, "plugin");
}

void RelationXML::setHDotVector(const SiconosVector&v)
{
  if (!hasHDot())
    hDotNode = SiconosDOMTreeTools::createVectorNode(rootNode, "hDot", v);
  else
    SiconosDOMTreeTools::setSiconosVectorNodeValue(hDotNode, v);
}

void RelationXML::setHDotPlugin(const string& plugin)
{
  if (! hDotNode)
  {
    hDotNode = SiconosDOMTreeTools::createSingleNode(rootNode, "hDot");
    xmlNewProp(hDotNode, (xmlChar*)"vectorPlugin", (xmlChar*)plugin.c_str());
  }
  else
    SiconosDOMTreeTools::setStringAttributeValue(hDotNode, "plugin", plugin);
}

void RelationXML::setJacobianHPlugin(const std::string&, unsigned int)
{
  XMLException::selfThrow("RelationXML -  setJacobianHPlugin: not yet implemented.");
}

string RelationXML::getJacobianHPlugin(unsigned int index) const
{
  if (index >= jacobianHNode.size())
    XMLException::selfThrow("RelationXML - getJacobianHPlugin(index), index out of range");

  if (!isJacobianHPlugin(index))
    XMLException::selfThrow("RelationXML -  getJacobianHPlugin: not computed with a plugin");

  return SiconosDOMTreeTools::getStringAttributeValue(jacobianHNode[index], "matrixPlugin");
}

SimpleMatrix RelationXML::getJacobianHMatrix(unsigned int index) const
{
  if (index >= jacobianHNode.size())
    XMLException::selfThrow("RelationXML - getJacobianHMatrix(index), index out of range");
  if (isJacobianHPlugin(index))
    XMLException::selfThrow("RelationXML - getJacobianHMatrix: JacobianH is computed using a plug-in");
  return  SiconosDOMTreeTools::getSiconosMatrixValue(jacobianHNode[index]);
}

void RelationXML::setJacobianGPlugin(const std::string&, unsigned int)
{
  XMLException::selfThrow("RelationXML -  setJacobianGPlugin: not yet implemented.");
}

string RelationXML::getJacobianGPlugin(unsigned int index) const
{
  if (index >= jacobianGNode.size())
    XMLException::selfThrow("RelationXML - getJacobianGPlugin(index), index out of range");

  if (!isJacobianGPlugin(index))
    XMLException::selfThrow("RelationXML -  getJacobianGPlugin: not computed with a plugin");

  return SiconosDOMTreeTools::getStringAttributeValue(jacobianGNode[index], "matrixPlugin");
}

SimpleMatrix RelationXML::getJacobianGMatrix(unsigned int index) const
{
  if (index >= jacobianGNode.size())
    XMLException::selfThrow("RelationXML - getJacobianGMatrix(index), index out of range");
  if (isJacobianGPlugin(index))
    XMLException::selfThrow("RelationXML - getJacobianGMatrix: JacobianG is computed using a plug-in");
  return  SiconosDOMTreeTools::getSiconosMatrixValue(jacobianGNode[index]);
}
