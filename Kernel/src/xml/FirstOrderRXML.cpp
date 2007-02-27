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
#include "FirstOrderRXML.h"
using namespace std;

FirstOrderRXML::FirstOrderRXML(): RelationXML(), computeInputNode(NULL), computeOutputNode(NULL)
{}

FirstOrderRXML::FirstOrderRXML(xmlNodePtr relationNode): RelationXML(relationNode), computeInputNode(NULL), computeOutputNode(NULL)
{
  xmlNodePtr node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, "computeInput")) != NULL)
    computeInputNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, "computeOutput")) != NULL)
    computeOutputNode = node;
}

FirstOrderRXML::~FirstOrderRXML()
{}

void FirstOrderRXML::setComputeInputPlugin(const string&  plugin)
{
  if (computeInputNode == NULL)
  {
    computeInputNode = SiconosDOMTreeTools::createSingleNode(rootNode, "computeInput");
    xmlNewProp(computeInputNode, (xmlChar*)"plugin", (xmlChar*)plugin.c_str());
  }
  else SiconosDOMTreeTools::setStringAttributeValue(computeInputNode, "plugin", plugin);
}

void FirstOrderRXML::setComputeOutputPlugin(const string&  plugin)
{
  if (computeOutputNode == NULL)
  {
    computeOutputNode = SiconosDOMTreeTools::createSingleNode(rootNode, "computeOutput");
    xmlNewProp(computeOutputNode, (xmlChar*)"plugin", (xmlChar*)plugin.c_str());
  }
  else SiconosDOMTreeTools::setStringAttributeValue(computeOutputNode, "plugin", plugin);
}

string FirstOrderRXML::getComputeInputPlugin() const
{
  if (!isComputeInputPlugin())
    XMLException::selfThrow("FirstOrderRXML - getComputeInputPlugin : computeInput is not calculated from a plugin");
  return  SiconosDOMTreeTools::getStringAttributeValue(computeInputNode, "plugin");
}

string FirstOrderRXML::getComputeOutputPlugin() const
{
  if (!isComputeOutputPlugin())
    XMLException::selfThrow("FirstOrderRXML - getComputeOutputPlugin : computeOutput is not calculated from a plugin");
  return  SiconosDOMTreeTools::getStringAttributeValue(computeOutputNode, "plugin");
}
