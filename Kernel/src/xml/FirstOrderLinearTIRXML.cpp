/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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
#include "FirstOrderLinearTIRXML.h"
using namespace std;

FirstOrderLinearTIRXML::FirstOrderLinearTIRXML(): FirstOrderRXML(), CNode(NULL), DNode(NULL), FNode(NULL), eNode(NULL), BNode(NULL)
{}

FirstOrderLinearTIRXML::FirstOrderLinearTIRXML(xmlNode * LTIRelationNode):
  FirstOrderRXML(LTIRelationNode), CNode(NULL), DNode(NULL), FNode(NULL), eNode(NULL), BNode(NULL)
{
  xmlNodePtr node;
  if ((node = SiconosDOMTreeTools::findNodeChild(LTIRelationNode, "C")) != NULL)
    CNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(LTIRelationNode, "D")) != NULL)
    DNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(LTIRelationNode, "F")) != NULL)
    FNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(LTIRelationNode, "e")) != NULL)
    eNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(LTIRelationNode, "B")) != NULL)
    BNode = node;
}

FirstOrderLinearTIRXML::~FirstOrderLinearTIRXML()
{}

void FirstOrderLinearTIRXML::setC(const SiconosMatrix& matrix)
{
  if (CNode == NULL)
    CNode = SiconosDOMTreeTools::createMatrixNode(rootNode, "C", matrix);
  else SiconosDOMTreeTools::setSiconosMatrixNodeValue(CNode, matrix);
}

void FirstOrderLinearTIRXML::setD(const SiconosMatrix& matrix)
{
  if (DNode == NULL)
    DNode = SiconosDOMTreeTools::createMatrixNode(rootNode, "D", matrix);
  else SiconosDOMTreeTools::setSiconosMatrixNodeValue(DNode, matrix);
}

void FirstOrderLinearTIRXML::setF(const SiconosMatrix &matrix)
{
  if (FNode == NULL)
    FNode = SiconosDOMTreeTools::createMatrixNode(rootNode, "F", matrix);
  else SiconosDOMTreeTools::setSiconosMatrixNodeValue(FNode, matrix);
}

void FirstOrderLinearTIRXML::setE(const SiconosVector& vec)
{
  if (eNode == NULL)
    eNode = SiconosDOMTreeTools::createVectorNode(rootNode, "e", vec);
  else SiconosDOMTreeTools::setSiconosVectorNodeValue(eNode, vec);
}

void FirstOrderLinearTIRXML::setB(const SiconosMatrix &matrix)
{
  if (BNode == NULL)
    BNode = SiconosDOMTreeTools::createMatrixNode(rootNode, "B", matrix);
  else SiconosDOMTreeTools::setSiconosMatrixNodeValue(BNode, matrix);
}

