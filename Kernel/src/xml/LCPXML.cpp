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
#include "LCPXML.h"
using namespace std;

LCPXML::LCPXML() : OneStepNSProblemXML()
{
  this->MNode = NULL;
  this->qNode = NULL;
}

LCPXML::LCPXML(xmlNode * LCPNode, vector<int> definedInteractionNumbers)
  : OneStepNSProblemXML(LCPNode, definedInteractionNumbers)
{
  xmlNode *node, *lcpModelNode;

  lcpModelNode = SiconosDOMTreeTools::findNodeChild(LCPNode);
  if (lcpModelNode != NULL)
  {
    if (strcmp((char*)lcpModelNode->name, OSNSP_SOLVER.c_str()) == 0)
      lcpModelNode = SiconosDOMTreeTools::findFollowNode(lcpModelNode);
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(lcpModelNode, LCP_M)) != NULL)
  {
    this->MNode = node;
  }
  else
    this->MNode = NULL;

  if ((node = SiconosDOMTreeTools::findNodeChild(lcpModelNode, LCP_Q)) != NULL)
  {
    this->qNode = node;
  }
  else
    this->qNode = NULL;
}


LCPXML::~LCPXML() {}

void LCPXML::updateOneStepNSProblemXML(xmlNode* node, OneStepNSProblem* osnspb)
{
  IN("LCPXML::updateOneStepNSProblemXML\n");
  this->rootNode = node;
  this->rootNSProblemXMLNode = SiconosDOMTreeTools::findNodeChild(this->rootNode);
  //this->loadOneStepNSProblem( osnspb );
  OUT("LCPXML::updateOneStepNSProblemXML\n");
}

