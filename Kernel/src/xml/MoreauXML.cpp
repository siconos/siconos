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
#include "MoreauXML.hpp"
using namespace std;

MoreauXML::MoreauXML(): OneStepIntegratorXML(), thetaNode(NULL), WNode(NULL)
{}

MoreauXML::MoreauXML(xmlNodePtr MoreauNode):
  OneStepIntegratorXML(MoreauNode), thetaNode(NULL), WNode(NULL)
{
  xmlNodePtr node;
  if ((node = SiconosDOMTreeTools::findNodeChild(MoreauNode, "W")))
    WNode = node;
  if ((node = SiconosDOMTreeTools::findNodeChild(MoreauNode, "Theta")))
    thetaNode = node;
}

MoreauXML::~MoreauXML()
{}

void MoreauXML::getTheta(vector<double>& values) const
{
  if (hasThetaList())
    SiconosDOMTreeTools::getVector(thetaNode, values);
  else
    XMLException::selfThrow("MoreauXML::getTheta - No list of theta in xml tag.");
}

void MoreauXML::setTheta(const vector<double>& v)
{
  XMLException::selfThrow("OneStepIntegratorXML::setTheta - not yet implemented.");
  //   if( !hasTheta()) // create the node if it does not exist
  //     thetaNode = SiconosDOMTreeTools::createVectorNode(rootNode,  "Theta", v);
  //   else SiconosDOMTreeTools::setSiconosVectorNodeValue(thetaNode, v);
}

const double MoreauXML::getSingleTheta() const
{
  if (!hasAllTheta())
    XMLException::selfThrow("MoreauXml getSingleTheta: the attribute all is not present in the tag theta");
  return SiconosDOMTreeTools::getAttributeValue<double>(thetaNode, "all");
}
