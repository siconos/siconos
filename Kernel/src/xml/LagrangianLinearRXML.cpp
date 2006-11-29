/* Siconos-Kernel version 2.0.0, Copyright INRIA 2005-2006.
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
#include "LagrangianLinearRXML.h"
using namespace std;

LagrangianLinearRXML::LagrangianLinearRXML():
  LagrangianRXML(), HNode(NULL), bNode(NULL), DNode(NULL)
{}

LagrangianLinearRXML::LagrangianLinearRXML(xmlNode * LLRelationNode)
  : LagrangianRXML(LLRelationNode), HNode(NULL), bNode(NULL), DNode(NULL)
{
  xmlNode *node;

  if ((node = SiconosDOMTreeTools::findNodeChild(LLRelationNode, LLR_H)) != NULL)
    HNode = node;
  else if (hNode == NULL)
    XMLException::selfThrow("LLRelationXML - constructor error : tag " + LLR_H + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(LLRelationNode, LLR_B)) != NULL)
    bNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(LLRelationNode, "D")) != NULL)
    DNode = node;
}

LagrangianLinearRXML::~LagrangianLinearRXML()
{}

void LagrangianLinearRXML::setH(const SiconosMatrix& matrix)
{
  if (HNode == NULL)
  {
    HNode = SiconosDOMTreeTools::createMatrixNode(rootRelationXMLNode, LLR_H, matrix);
  }
  else SiconosDOMTreeTools::setSiconosMatrixNodeValue(HNode, matrix);
}

void LagrangianLinearRXML::setB(const SiconosVector &vec)
{
  if (bNode == NULL)
  {
    bNode = SiconosDOMTreeTools::createVectorNode(rootRelationXMLNode, LLR_B, vec);
  }
  else SiconosDOMTreeTools::setSiconosVectorNodeValue(bNode, vec);
}

void LagrangianLinearRXML::setD(const SiconosMatrix& matrix)
{
  if (DNode == NULL)
  {
    DNode = SiconosDOMTreeTools::createMatrixNode(rootRelationXMLNode, "D", matrix);
  }
  else SiconosDOMTreeTools::setSiconosMatrixNodeValue(DNode, matrix);
}
