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
#include "LagrangianLinearRXML.h"
using namespace std;

LagrangianLinearRXML::LagrangianLinearRXML():
  LagrangianRXML(), HNode(NULL), bNode(NULL), DNode(NULL), FNode(NULL)
{}

LagrangianLinearRXML::LagrangianLinearRXML(xmlNode * LLRelationNode)
  : LagrangianRXML(LLRelationNode), HNode(NULL), bNode(NULL), DNode(NULL), FNode(NULL)
{
  xmlNode *node;

  if ((node = SiconosDOMTreeTools::findNodeChild(LLRelationNode, "H")))
    HNode = node;
  else if (!hNode)
    XMLException::selfThrow("LLRelationXML - constructor error : tag H not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(LLRelationNode, "b")))
    bNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(LLRelationNode, "D")))
    DNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(LLRelationNode, "F")))
    FNode = node;
}

LagrangianLinearRXML::~LagrangianLinearRXML()
{}

void LagrangianLinearRXML::setH(const SiconosMatrix& matrix)
{
  if (!HNode)
  {
    HNode = SiconosDOMTreeTools::createMatrixNode(rootNode, "H", matrix);
  }
  else SiconosDOMTreeTools::setSiconosMatrixNodeValue(HNode, matrix);
}

void LagrangianLinearRXML::setB(const SiconosVector &vec)
{
  if (!bNode)
  {
    bNode = SiconosDOMTreeTools::createVectorNode(rootNode, "b", vec);
  }
  else SiconosDOMTreeTools::setSiconosVectorNodeValue(bNode, vec);
}

void LagrangianLinearRXML::setD(const SiconosMatrix& matrix)
{
  if (!DNode)
  {
    DNode = SiconosDOMTreeTools::createMatrixNode(rootNode, "D", matrix);
  }
  else SiconosDOMTreeTools::setSiconosMatrixNodeValue(DNode, matrix);
}

void LagrangianLinearRXML::setF(const SiconosMatrix& matrix)
{
  if (!FNode)
  {
    FNode = SiconosDOMTreeTools::createMatrixNode(rootNode, "F", matrix);
  }
  else SiconosDOMTreeTools::setSiconosMatrixNodeValue(FNode, matrix);
}
