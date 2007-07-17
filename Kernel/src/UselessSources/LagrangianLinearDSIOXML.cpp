/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2006.
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
#include "LagrangianLinearDSIOXML.h"
using namespace std;

LagrangianLinearDSIOXML::LagrangianLinearDSIOXML(): LagrangianDSIOXML()
{
  this->HNode = NULL;
  this->bNode = NULL;
}

LagrangianLinearDSIOXML::LagrangianLinearDSIOXML(xmlNode * dsioNode)
// : DSInputOutputXML(dsioNode)
{
  xmlNode *node;
  this->rootDSIOXMLNode = dsioNode;
  if ((node = SiconosDOMTreeTools::findNodeChild(dsioNode, LLDSIO_H)) != NULL)
    this->HNode = node;
  else
    XMLException::selfThrow("LagrangianLinearDSIOXML - constructor error : tag " + LLDSIO_H + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(dsioNode, LLDSIO_B)) != NULL)
    this->bNode = node;
  else
    XMLException::selfThrow("LagrangianLinearDSIOXML - constructor error : tag " + LLDSIO_B + " not found.");
}

LagrangianLinearDSIOXML::~LagrangianLinearDSIOXML()
{}

void LagrangianLinearDSIOXML::setH(SiconosMatrix *matrix)
{
  if (this->HNode == NULL)
  {
    this->HNode = SiconosDOMTreeTools::createMatrixNode(this->rootDSIOXMLNode, LLDSIO_H, *matrix);
  }
  else SiconosDOMTreeTools::setSiconosMatrixNodeValue(this->HNode, *matrix);
}

void LagrangianLinearDSIOXML::setB(SiconosVector *vector)
{
  if (this->bNode == NULL)
  {
    this->bNode = SiconosDOMTreeTools::createVectorNode(this->rootDSIOXMLNode, LLDSIO_B, *vector);
  }
  else SiconosDOMTreeTools::setSiconosVectorNodeValue(this->bNode, *vector);
}
