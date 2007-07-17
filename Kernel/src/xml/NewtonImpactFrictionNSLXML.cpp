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
#include "NewtonImpactFrictionNSLXML.h"
using namespace std;

NewtonImpactFrictionNSLXML::NewtonImpactFrictionNSLXML():
  NonSmoothLawXML(), enNode(NULL), etNode(NULL), muNode(NULL)
{}

NewtonImpactFrictionNSLXML::NewtonImpactFrictionNSLXML(xmlNodePtr nslNode):
  NonSmoothLawXML(nslNode), enNode(NULL), etNode(NULL), muNode(NULL)
{
  xmlNodePtr node;
  if ((node = SiconosDOMTreeTools::findNodeChild(nslNode, "en")) != NULL)
    enNode = node;
  else
    XMLException::selfThrow("NewtonImpactFrictionNSLXML - constructor error : tag en not found.");

  // et is not always required (friction2D)
  if ((node = SiconosDOMTreeTools::findNodeChild(nslNode, "et")) != NULL)
    etNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(nslNode, "mu")) != NULL)
    muNode = node;
  else
    XMLException::selfThrow("NewtonImpactFrictionNSLXML - constructor error : tag mu not found.");
}
