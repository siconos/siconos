/* Siconos-Kernel, Copyright INRIA 2005-2010.
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
#include "FirstOrderLinearDSXML.hpp"
using namespace std;

FirstOrderLinearDSXML::FirstOrderLinearDSXML() :
  FirstOrderNonLinearDSXML(), ANode(NULL), bNode(NULL)
{}

FirstOrderLinearDSXML::FirstOrderLinearDSXML(xmlNodePtr nodeDS, bool isBVP):
  FirstOrderNonLinearDSXML(nodeDS, isBVP), ANode(NULL), bNode(NULL)
{
  xmlNodePtr node = SiconosDOMTreeTools::findNodeChild(rootNode, LDS_A);
  if (node)
    ANode = node;
  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, LDS_B)))
    bNode = node;
}

void FirstOrderLinearDSXML::setA(const SiconosMatrix& m)
{
  if (ANode)
    SiconosDOMTreeTools::setSiconosMatrixNodeValue(ANode, m);
  else
    ANode = SiconosDOMTreeTools::createMatrixNode(rootNode, LDS_A, m);
}

void FirstOrderLinearDSXML::setAPlugin(const std::string& plugin)
{
  if (!ANode)
  {
    ANode = SiconosDOMTreeTools::createSingleNode(rootNode, "A");
    xmlNewProp(ANode, (xmlChar*)("matrixPlugin"), (xmlChar*)plugin.c_str());
  }
  else
    SiconosDOMTreeTools::setStringAttributeValue(ANode, "matrixPlugin", plugin);
}

void FirstOrderLinearDSXML::setB(const SiconosVector& v)
{
  if (bNode)
    SiconosDOMTreeTools::setSiconosVectorNodeValue(bNode, v);
  else bNode = SiconosDOMTreeTools::createVectorNode(rootNode, LDS_B, v);
}

void FirstOrderLinearDSXML::setBPlugin(const std::string& plugin)
{
  if (!bNode)
  {
    bNode = SiconosDOMTreeTools::createSingleNode(rootNode, "b");
    xmlNewProp(bNode, (xmlChar*)("vectorPlugin"), (xmlChar*)plugin.c_str());
  }
  else
    SiconosDOMTreeTools::setStringAttributeValue(bNode, "vectorPlugin", plugin);
}

