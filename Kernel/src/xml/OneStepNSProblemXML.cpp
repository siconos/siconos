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

#include "OneStepNSProblemXML.h"
#include "NonSmoothSolverXML.h"

using namespace std;

OneStepNSProblemXML::OneStepNSProblemXML(xmlNodePtr oneStepNSProblemXMLNode):
  rootNode(oneStepNSProblemXMLNode), dimNode(NULL), interactionsConcernedNode(NULL),
  solverNode(NULL), solverXML(NULL), isSolverXMLAllocatedIn(false)
{
  // rootNode == formalisation type (LCP ...)

  xmlNodePtr node;
  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, "size")))
    dimNode = node;

  // get interactionsConcerned
  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, "Interactions_Concerned")))
    interactionsConcernedNode = node;

  // Solver
  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, "NonSmoothSolver")))
  {
    solverNode = node;
    solverXML = new NonSmoothSolverXML(solverNode);
    isSolverXMLAllocatedIn = true;
  }
  //  else -> nothing, it's up to OneStepNSProblem constructor to deal with the fact that no solver has been given
}

OneStepNSProblemXML::~OneStepNSProblemXML()
{
  if (isSolverXMLAllocatedIn) delete solverXML;
  solverXML = NULL;
}

void OneStepNSProblemXML::setDimNSProblem(int n)
{
  if (! hasDim())
    dimNode = SiconosDOMTreeTools::createIntegerNode(rootNode, "size", n);
  else SiconosDOMTreeTools::setIntegerContentValue(dimNode, n);
}

bool OneStepNSProblemXML::hasAllInteractions() const
{
  if (SiconosDOMTreeTools::hasAttributeValue(interactionsConcernedNode, ALL_ATTRIBUTE))
    return SiconosDOMTreeTools::getAttributeValue<bool>(interactionsConcernedNode, ALL_ATTRIBUTE);
  else return false;
}

void OneStepNSProblemXML::setAllInteractions(const bool& all)
{
  if (!hasAllInteractions())
  {
    if (all == true)
      xmlNewProp(interactionsConcernedNode, (xmlChar*)ALL_ATTRIBUTE.c_str(), (xmlChar*)"true");
  }
  else
  {
    if (all == false)
      xmlRemoveProp(xmlHasProp(interactionsConcernedNode, (xmlChar*)ALL_ATTRIBUTE.c_str()));
  }
}

void OneStepNSProblemXML::setNonSmoothSolverXMLPtr(NonSmoothSolverXML * solv)
{
  if (isSolverXMLAllocatedIn) delete solverXML;
  solverXML = solv;
  isSolverXMLAllocatedIn = false;
}

void OneStepNSProblemXML::updateOneStepNSProblemXML(xmlNode* node, SP::OneStepNSProblem osnspb)
{
  rootNode = node;
}

void OneStepNSProblemXML::getInteractionsNumbers(vector<int>& inNumbers)
{
  if (!hasAllInteractions())
    SiconosDOMTreeTools::getVector(interactionsConcernedNode, inNumbers);
  else
    XMLException::selfThrow("OneStepNSProblemXML::getInteractionsNumbers - The list of interactions is missing.");
}
