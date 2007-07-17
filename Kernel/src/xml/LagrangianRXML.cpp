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
#include "LagrangianRXML.h"
#include "SimpleMatrix.h"

using namespace std;

LagrangianRXML::LagrangianRXML(): RelationXML(), hNode(NULL), hDotNode(NULL)
{}

LagrangianRXML::LagrangianRXML(xmlNodePtr LNLRelationNode): RelationXML(LNLRelationNode), hNode(NULL), hDotNode(NULL)
{
  xmlNodePtr node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, "h")) != NULL)
    hNode = node;
  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, "hDot")) != NULL)
    hDotNode = node;
  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, "G")) != NULL)
  {
    // get number of G functions given
    unsigned int size = SiconosDOMTreeTools::getAttributeValue<unsigned int>(node, "number");
    if (size == 0)
      XMLException::selfThrow("LagrangianRXML:: constructor failed. Some Gradients plug-in names are required (G).");
    GNode.resize(size, NULL);
    // get corresponding nodes
    GNode[0] = SiconosDOMTreeTools::findNodeChild(node, "matrix")  ;
    if (GNode[0] == NULL)
      XMLException::selfThrow("LagrangianRXML:: constructor, G0 is missing");
    for (unsigned int i = 1; i < size; i++)
    {
      GNode[i] =  SiconosDOMTreeTools::findFollowNode(GNode[i - 1]);
      if (GNode[i] == NULL)
        XMLException::selfThrow("LagrangianRXML:: constructor, another G is required");
    }
  }
}

LagrangianRXML::~LagrangianRXML()
{}

unsigned int LagrangianRXML::getNumberOfGradients() const
{
  xmlNodePtr node = SiconosDOMTreeTools::findNodeChild(rootNode, "G");
  if (node == NULL)
    return 0;
  else
    return SiconosDOMTreeTools::getAttributeValue<unsigned int>(node, "number");
}

void LagrangianRXML::setHDotVector(const SiconosVector&v)
{
  if (!hasHDot())
    hDotNode = SiconosDOMTreeTools::createVectorNode(rootNode, "hDot", v);
  else
    SiconosDOMTreeTools::setSiconosVectorNodeValue(hDotNode, v);
}

void LagrangianRXML::setHDotPlugin(const string& plugin)
{
  if (hDotNode == NULL)
  {
    hDotNode = SiconosDOMTreeTools::createSingleNode(rootNode, "hDot");
    xmlNewProp(hDotNode, (xmlChar*)"vectorPlugin", (xmlChar*)plugin.c_str());
  }
  else
    SiconosDOMTreeTools::setStringAttributeValue(hDotNode, "plugin", plugin);
}

bool LagrangianRXML::isGPlugin(unsigned int index) const
{
  if (index >= GNode.size())
    XMLException::selfThrow("LagrangianRXML - isGPlugin(index), index out of range");

  return xmlHasProp(GNode[index], (xmlChar *)"matrixPlugin");
}

bool LagrangianRXML::hasG(unsigned int index) const
{
  if (index >= GNode.size())
    XMLException::selfThrow("LagrangianRXML - hasG(index), index out of range");
  return (GNode[index] != NULL);
}

string LagrangianRXML::getGPlugin(unsigned int index) const
{
  if (index >= GNode.size())
    XMLException::selfThrow("LagrangianRXML - getGPlugin(index), index out of range");

  if (!isGPlugin(index))
    XMLException::selfThrow("LagrangianRXML - getGPlugin : G is not computed with a plugin");

  return SiconosDOMTreeTools::getStringAttributeValue(GNode[index], "matrixPlugin");
}

SimpleMatrix LagrangianRXML::getGMatrix(unsigned int index) const
{
  if (index >= GNode.size())
    XMLException::selfThrow("LagrangianRXML - getGMatrix(index), index out of range");
  if (isGPlugin(index))
    XMLException::selfThrow("LagrangianRXML - getGMatrix: G is computed using a plug-in");
  return  SiconosDOMTreeTools::getSiconosMatrixValue(GNode[index]);
}

void LagrangianRXML::setGPlugin(const std::string& plugin, unsigned int index)
{
  // \todo
  XMLException::selfThrow("LagrangianRXML - setGPlugin(...), not yet implemented.");
}

void LagrangianRXML::setGMatrix(SiconosMatrix *newMat, unsigned int index)
{
  // \todo
  XMLException::selfThrow("LagrangianRXML - setGMatrix(...), not yet implemented.");
}



