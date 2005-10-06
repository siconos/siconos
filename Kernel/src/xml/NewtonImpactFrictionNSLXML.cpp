/* Siconos version 1.0, Copyright INRIA 2005.
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

NewtonImpactFrictionNSLXML::NewtonImpactFrictionNSLXML()
{
  this->enNode = NULL;
  this->etNode = NULL;
  this->muNode = NULL;
}

NewtonImpactFrictionNSLXML::NewtonImpactFrictionNSLXML(xmlNode * nslNode)
  : NonSmoothLawXML(nslNode)
{
  xmlNode *node;

  if ((node = SiconosDOMTreeTools::findNodeChild(nslNode, NEWTON_EN)) != NULL)
  {
    this->enNode = node;
  }
  else
  {
    XMLException::selfThrow("NewtonImpactFrictionNSLXML - constructor error : tag " + NEWTON_EN + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(nslNode, NEWTON_ET)) != NULL)
  {
    this->etNode = node;
  }
  else
  {
    XMLException::selfThrow("NewtonImpactFrictionNSLXML - constructor error : tag " + NEWTON_ET + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(nslNode, NEWTON_MU)) != NULL)
  {
    this->muNode = node;
  }
  else
  {
    XMLException::selfThrow("NewtonImpactFrictionNSLXML - constructor error : tag " + NEWTON_MU + " not found.");
  }

}
