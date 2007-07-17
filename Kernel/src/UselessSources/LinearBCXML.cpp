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
#include "LinearBCXML.h"
using namespace std;

LinearBCXML::LinearBCXML() : BoundaryConditionXML()
{
  this->omegaNode = NULL;
  this->omega0Node = NULL;
  this->omegaTNode = NULL;
}


LinearBCXML::LinearBCXML(xmlNode * linearBCNode) : BoundaryConditionXML(linearBCNode)
{
  this->loadLinearBCProperties();
}

LinearBCXML::~LinearBCXML()
{}

void LinearBCXML::loadLinearBCProperties()
{
  xmlNode *node;

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootBCNode, LINEARBC_OMEGA)) != NULL)
  {
    this->omegaNode = node;
  }
  else
  {
    XMLException::selfThrow("LinearBCXML - loadLinearBCProperties error : tag " + LINEARBC_OMEGA + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootBCNode, LINEARBC_OMEGA0)) != NULL)
  {
    this->omega0Node = node;
  }
  else
  {
    XMLException::selfThrow("LinearBCXML - loadLinearBCProperties error : tag " + LINEARBC_OMEGA0 + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootBCNode, LINEARBC_OMEGAT)) != NULL)
  {
    this->omegaTNode = node;
  }
  else
  {
    XMLException::selfThrow("LinearBCXML - loadLinearBCProperties error : tag " + LINEARBC_OMEGAT + " not found.");
  }
}

void LinearBCXML::updateBoundaryConditionXML(xmlNode* node) //, BoundaryCondition* bc)
{
  this->rootBCNode = node;
}

