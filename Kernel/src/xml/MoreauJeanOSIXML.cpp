/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
#include "MoreauJeanOSIXML.hpp"


MoreauJeanOSIXML::MoreauJeanOSIXML(): OneStepIntegratorXML(), thetaNode(NULL), WNode(NULL)
{}

MoreauJeanOSIXML::MoreauJeanOSIXML(xmlNode * MoreauJeanOSINode):
  OneStepIntegratorXML(MoreauJeanOSINode), thetaNode(NULL), WNode(NULL)
{
  xmlNodePtr node;
  if ((node = SiconosDOMTreeTools::findNodeChild(MoreauJeanOSINode, "W")))
    WNode = node;
  if ((node = SiconosDOMTreeTools::findNodeChild(MoreauJeanOSINode, "Theta")))
    thetaNode = node;
}

MoreauJeanOSIXML::~MoreauJeanOSIXML()
{}

void MoreauJeanOSIXML::getTheta(std::vector<double>& values) const
{
  if (hasThetaList())
    SiconosDOMTreeTools::getVector(thetaNode, values);
  else
    XMLException::selfThrow("MoreauJeanOSIXML::getTheta - No list of theta in xml tag.");
}

void MoreauJeanOSIXML::setTheta(const std::vector<double>& v)
{
  XMLException::selfThrow("OneStepIntegratorXML::setTheta - not yet implemented.");
  //   if( !hasTheta()) // create the node if it does not exist
  //     thetaNode = SiconosDOMTreeTools::createVectorNode(rootNode,  "Theta", v);
  //   else SiconosDOMTreeTools::setSiconosVectorNodeValue(thetaNode, v);
}

double MoreauJeanOSIXML::getSingleTheta() const
{
  if (!hasAllTheta())
    XMLException::selfThrow("MoreauJeanOSIXml getSingleTheta: the attribute all is not present in the tag theta");
  return SiconosDOMTreeTools::getAttributeValue<double>(thetaNode, "all");
}
