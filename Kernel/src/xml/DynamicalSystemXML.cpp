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
#include "DynamicalSystemXML.hpp"
using namespace std;
using namespace DS;

DynamicalSystemXML::DynamicalSystemXML(xmlNodePtr DSNode, bool):
  rootNode(DSNode), stepsInMemoryNode(NULL), zNode(NULL)
{
  xmlNodePtr node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, "StepsInMemory")))
    stepsInMemoryNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, "z")))
    zNode = node;
}

const DS::TYPES DynamicalSystemXML::getType() const
{
  std::string res((char*)rootNode->name);
  if (res == "FirstOrderNonLinearDS")
    return FONLDS;
  else if (res == "FirstOrderLinearDS")
    return FOLDS;
  else if (res == "FirstOrderLinearTIDS")
    return FOLTIDS;
  else if (res == "LagrangianDS")
    return LNLDS;
  else if (res == "LagrangianLinearTIDS")
    return LLTIDS;
  else
  {
    XMLException::selfThrow("DynamicalSystemXML - getType: unknown type of DS.");
    return FONLDS;
  }
}

void DynamicalSystemXML::setStepsInMemory(const unsigned int& nb)
{
  if (!hasStepsInMemory())
    stepsInMemoryNode = SiconosDOMTreeTools::createIntegerNode(rootNode, "StepsInMemory", nb);
  else SiconosDOMTreeTools::setIntegerContentValue(stepsInMemoryNode, nb);
}

