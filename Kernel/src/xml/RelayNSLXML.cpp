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

#include "RelayNSLXML.h"
using namespace std;

RelayNSLXML::RelayNSLXML(): CNode(NULL), DNode(NULL)
{}

RelayNSLXML::RelayNSLXML(xmlNode * relayNSLawNode)
  : NonSmoothLawXML(relayNSLawNode), CNode(NULL), DNode(NULL)
{
  xmlNode *node;

  if ((node = SiconosDOMTreeTools::findNodeChild(relayNSLawNode, RNSL_C)) != NULL)
    CNode = node;
  else
    XMLException::selfThrow("RelayNSLawXML - constructor error : tag " + RNSL_C + " not found.");

  if ((node = SiconosDOMTreeTools::findNodeChild(relayNSLawNode, RNSL_D)) != NULL)
    DNode = node;
  else
    XMLException::selfThrow("RelayNSLawXML - constructor error : tag " + RNSL_D + " not found.");
}
