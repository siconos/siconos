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
#include "LCPXML.h"
using namespace std;

LCPXML::LCPXML() : OneStepNSProblemXML(), MNode(NULL), qNode(NULL)
{}

LCPXML::LCPXML(xmlNode * OSNSNode):
  OneStepNSProblemXML(OSNSNode), MNode(NULL), qNode(NULL)
{
  // OSNSNode == "OneStepNSProblem"
  // problemTypeNode = "LCP"

  xmlNode* node;
  // dim
  //if ((node=SiconosDOMTreeTools::findNodeChild(problemTypeNode, "n")) !=NULL)
  //  dimNode=node;

  // M
  if ((node = SiconosDOMTreeTools::findNodeChild(problemTypeNode, "M")) != NULL)
    MNode = node;

  // q
  if ((node = SiconosDOMTreeTools::findNodeChild(problemTypeNode, "q")) != NULL)
    qNode = node;

}


LCPXML::~LCPXML() {}

void LCPXML::setM(const SiconosMatrix &m)
{
  if (!hasM())
    MNode = SiconosDOMTreeTools::createMatrixNode(problemTypeNode, "M", m);
  else SiconosDOMTreeTools::setSiconosMatrixNodeValue(MNode, m);
}

void LCPXML::setQ(const SiconosVector& q)
{
  if (!hasQ())
    qNode = SiconosDOMTreeTools::createVectorNode(problemTypeNode, "q", q);
  else SiconosDOMTreeTools::setSiconosVectorNodeValue(qNode, q);
}
