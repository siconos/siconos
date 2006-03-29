/* Siconos-Kernel version 1.1.4, Copyright INRIA 2005-2006.
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
#include "LinearTIRXML.h"
using namespace std;

LinearTIRXML::LinearTIRXML():
  RelationXML(), CNode(NULL), DNode(NULL), FNode(NULL), eNode(NULL), BNode(NULL)
{}

LinearTIRXML::LinearTIRXML(xmlNode * LTIRelationNode):
  RelationXML(LTIRelationNode), CNode(NULL), DNode(NULL), FNode(NULL), eNode(NULL), BNode(NULL)
{
  xmlNode *node;
  if ((node = SiconosDOMTreeTools::findNodeChild(LTIRelationNode, LTIR_C)) != NULL)
    CNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(LTIRelationNode, LTIR_D)) != NULL)
    DNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(LTIRelationNode, LTIR_F)) != NULL)
    FNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(LTIRelationNode, LTIR_E)) != NULL)
    eNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(LTIRelationNode, LTIR_B)) != NULL)
    BNode = node;
}

LinearTIRXML::~LinearTIRXML()
{}

void LinearTIRXML::setC(const SiconosMatrix& matrix)
{
  if (CNode == NULL)
    CNode = SiconosDOMTreeTools::createMatrixNode(rootRelationXMLNode, LTIR_C, matrix);
  else SiconosDOMTreeTools::setSiconosMatrixNodeValue(CNode, matrix);
}

void LinearTIRXML::setD(const SiconosMatrix& matrix)
{
  if (DNode == NULL)
    DNode = SiconosDOMTreeTools::createMatrixNode(rootRelationXMLNode, LTIR_D, matrix);
  else SiconosDOMTreeTools::setSiconosMatrixNodeValue(DNode, matrix);
}

void LinearTIRXML::setF(const SiconosMatrix &matrix)
{
  if (FNode == NULL)
    FNode = SiconosDOMTreeTools::createMatrixNode(rootRelationXMLNode, LTIR_F, matrix);
  else SiconosDOMTreeTools::setSiconosMatrixNodeValue(FNode, matrix);
}

void LinearTIRXML::setE(const SiconosVector& vec)
{
  if (eNode == NULL)
    eNode = SiconosDOMTreeTools::createVectorNode(rootRelationXMLNode, LTIR_E, vec);
  else SiconosDOMTreeTools::setSiconosVectorNodeValue(eNode, vec);
}

void LinearTIRXML::setB(const SiconosMatrix &matrix)
{
  if (BNode == NULL)
    BNode = SiconosDOMTreeTools::createMatrixNode(rootRelationXMLNode, LTIR_B, matrix);
  else SiconosDOMTreeTools::setSiconosMatrixNodeValue(BNode, matrix);
}

