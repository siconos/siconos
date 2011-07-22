/* Siconos-Kernel, Copyright INRIA 2005-2011.
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

#include "OneStepNSProblemXML.hpp"


using namespace std;

OneStepNSProblemXML::OneStepNSProblemXML(xmlNodePtr oneStepNSProblemXMLNode):
  rootNode(oneStepNSProblemXMLNode), dimNode(NULL), interactionsConcernedNode(NULL), numericsSolverNameNode(NULL)
{
  // rootNode == formalisation type (LCP ...)

  xmlNodePtr node;
  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, "size")))
    dimNode = node;

  // get interactionsConcerned
  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, "Interactions_Concerned")))
    interactionsConcernedNode = node;

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, "NumericsSolverName")))
    numericsSolverNameNode = node;




  //  else -> nothing, it's up to OneStepNSProblem constructor to deal with the fact that no solver has been given
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


void OneStepNSProblemXML::getInteractionsNumbers(vector<int>& inNumbers)
{
  if (!hasAllInteractions())
    SiconosDOMTreeTools::getVector(interactionsConcernedNode, inNumbers);
  else
    XMLException::selfThrow("OneStepNSProblemXML::getInteractionsNumbers - The list of interactions is missing.");
}
