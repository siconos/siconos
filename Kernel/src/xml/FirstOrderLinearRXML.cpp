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
#include "FirstOrderLinearRXML.h"
using namespace std;

FirstOrderLinearRXML::FirstOrderLinearRXML(): FirstOrderRXML(), CNode(NULL), DNode(NULL), FNode(NULL), eNode(NULL), BNode(NULL)
{}

FirstOrderLinearRXML::FirstOrderLinearRXML(xmlNode * LTIRelationNode):
  FirstOrderRXML(LTIRelationNode), CNode(NULL), DNode(NULL), FNode(NULL), eNode(NULL), BNode(NULL)
{
  xmlNodePtr node;
  if ((node = SiconosDOMTreeTools::findNodeChild(LTIRelationNode, "C")))
    CNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(LTIRelationNode, "D")))
    DNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(LTIRelationNode, "F")))
    FNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(LTIRelationNode, "e")))
    eNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(LTIRelationNode, "B")))
    BNode = node;
}

FirstOrderLinearRXML::~FirstOrderLinearRXML()
{}

void FirstOrderLinearRXML::setC(const SiconosMatrix& matrix)
{
  if (!CNode)
    CNode = SiconosDOMTreeTools::createMatrixNode(rootNode, "C", matrix);
  else SiconosDOMTreeTools::setSiconosMatrixNodeValue(CNode, matrix);
}

void FirstOrderLinearRXML::setD(const SiconosMatrix& matrix)
{
  if (!DNode)
    DNode = SiconosDOMTreeTools::createMatrixNode(rootNode, "D", matrix);
  else SiconosDOMTreeTools::setSiconosMatrixNodeValue(DNode, matrix);
}

void FirstOrderLinearRXML::setF(const SiconosMatrix &matrix)
{
  if (!FNode)
    FNode = SiconosDOMTreeTools::createMatrixNode(rootNode, "F", matrix);
  else SiconosDOMTreeTools::setSiconosMatrixNodeValue(FNode, matrix);
}

void FirstOrderLinearRXML::setE(const SiconosVector& vec)
{
  if (!eNode)
    eNode = SiconosDOMTreeTools::createVectorNode(rootNode, "e", vec);
  else SiconosDOMTreeTools::setSiconosVectorNodeValue(eNode, vec);
}

void FirstOrderLinearRXML::setB(const SiconosMatrix &matrix)
{
  if (! BNode)
    BNode = SiconosDOMTreeTools::createMatrixNode(rootNode, "B", matrix);
  else SiconosDOMTreeTools::setSiconosMatrixNodeValue(BNode, matrix);
}

void FirstOrderLinearRXML::setCPlugin(const std::string& plugin)
{
  if (!CNode)
  {
    CNode = SiconosDOMTreeTools::createSingleNode(rootNode, "C");
    xmlNewProp(CNode, (xmlChar*)("matrixPlugin"), (xmlChar*)plugin.c_str());
  }
  else
    SiconosDOMTreeTools::setStringAttributeValue(CNode, "matrixPlugin", plugin);
}

void FirstOrderLinearRXML::setDPlugin(const std::string& plugin)
{
  if (!DNode)
  {
    DNode = SiconosDOMTreeTools::createSingleNode(rootNode, "D");
    xmlNewProp(DNode, (xmlChar*)("matrixPlugin"), (xmlChar*)plugin.c_str());
  }
  else
    SiconosDOMTreeTools::setStringAttributeValue(DNode, "matrixPlugin", plugin);
}

void FirstOrderLinearRXML::setFPlugin(const std::string& plugin)
{
  if (!FNode)
  {
    FNode = SiconosDOMTreeTools::createSingleNode(rootNode, "F");
    xmlNewProp(FNode, (xmlChar*)("matrixPlugin"), (xmlChar*)plugin.c_str());
  }
  else
    SiconosDOMTreeTools::setStringAttributeValue(FNode, "matrixPlugin", plugin);
}

void FirstOrderLinearRXML::setEPlugin(const std::string& plugin)
{
  if (!eNode)
  {
    eNode = SiconosDOMTreeTools::createSingleNode(rootNode, "e");
    xmlNewProp(eNode, (xmlChar*)("vectorPlugin"), (xmlChar*)plugin.c_str());
  }
  else
    SiconosDOMTreeTools::setStringAttributeValue(eNode, "vectorPlugin", plugin);
}

void FirstOrderLinearRXML::setBPlugin(const std::string& plugin)
{
  if (!BNode)
  {
    BNode = SiconosDOMTreeTools::createSingleNode(rootNode, "B");
    xmlNewProp(BNode, (xmlChar*)("matrixPlugin"), (xmlChar*)plugin.c_str());
  }
  else
    SiconosDOMTreeTools::setStringAttributeValue(BNode, "matrixPlugin", plugin);
}
