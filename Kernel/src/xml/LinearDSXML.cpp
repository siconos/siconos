/* Siconos-Kernel version 1.3.0, Copyright INRIA 2005-2006.
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
#include "LinearDSXML.h"
using namespace std;

LinearDSXML::LinearDSXML() :
  DynamicalSystemXML(), ANode(NULL), MxdotNode(NULL), bNode(NULL)
{}

LinearDSXML::LinearDSXML(xmlNode * LinearDSNode, const bool& isBVP):
  DynamicalSystemXML(LinearDSNode, isBVP), ANode(NULL), MxdotNode(NULL), bNode(NULL)
{
  xmlNode *node;
  // The only required node is A
  node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, LDS_A);
  if (node != NULL)
    ANode = node;
  else
    XMLException::selfThrow("LinearDSXML - loadLinearDSProperties error : tag " + LDS_A + " not found.");
  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, LDS_Mxdot)) != NULL)
    MxdotNode = node;
  if ((node = SiconosDOMTreeTools::findNodeChild(rootDynamicalSystemXMLNode, LDS_B)) != NULL)
    bNode = node;
}

LinearDSXML::~LinearDSXML()
{}

void LinearDSXML::setA(const SiconosMatrix& m)
{
  if (ANode != NULL)
    SiconosDOMTreeTools::setSiconosMatrixNodeValue(ANode, m);
  else
    ANode = SiconosDOMTreeTools::createMatrixNode(rootDynamicalSystemXMLNode, LDS_A, m);
}

void LinearDSXML::setAPlugin(const std::string& plugin)
{
  if (ANode == NULL)
  {
    ANode = SiconosDOMTreeTools::createSingleNode(rootDynamicalSystemXMLNode, "A");
    xmlNewProp(ANode, (xmlChar*)("matrixPlugin"), (xmlChar*)plugin.c_str());
  }
  else
    SiconosDOMTreeTools::setStringAttributeValue(ANode, "matrixPlugin", plugin);
}

void LinearDSXML::setMxdot(const SiconosMatrix& m)
{
  if (MxdotNode != NULL)
    SiconosDOMTreeTools::setSiconosMatrixNodeValue(MxdotNode, m);
  else
    MxdotNode = SiconosDOMTreeTools::createMatrixNode(rootDynamicalSystemXMLNode, LDS_Mxdot, m);
}

void LinearDSXML::setB(const SiconosVector& v)
{
  if (bNode != NULL)
    SiconosDOMTreeTools::setSiconosVectorNodeValue(bNode, v);
  else bNode = SiconosDOMTreeTools::createVectorNode(rootDynamicalSystemXMLNode, LDS_B, v);
}

void LinearDSXML::setBPlugin(const std::string& plugin)
{
  if (bNode == NULL)
  {
    bNode = SiconosDOMTreeTools::createSingleNode(rootDynamicalSystemXMLNode, "b");
    xmlNewProp(bNode, (xmlChar*)("vectorPlugin"), (xmlChar*)plugin.c_str());
  }
  else
    SiconosDOMTreeTools::setStringAttributeValue(bNode, "vectorPlugin", plugin);
}

