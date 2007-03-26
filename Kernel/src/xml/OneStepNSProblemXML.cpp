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

#include "OneStepNSProblemXML.h"
#include "SolverXML.h"

using namespace std;


OneStepNSProblemXML::OneStepNSProblemXML():
  rootNode(NULL), problemTypeNode(NULL), dimNode(NULL), interactionConcernedNode(NULL),
  interactionListNode(NULL), solverNode(NULL), solverXML(NULL), isSolverXMLAllocatedIn(false)
{}

OneStepNSProblemXML::OneStepNSProblemXML(xmlNode * oneStepNSProblemXMLNode):
  rootNode(oneStepNSProblemXMLNode), problemTypeNode(NULL), dimNode(NULL), interactionConcernedNode(NULL),
  interactionListNode(NULL), solverNode(NULL), solverXML(NULL), isSolverXMLAllocatedIn(false)
{
  // Two steps in OneStepNSProblem xml loading:
  //  - problem formalisation part ( LCP, FrictionContact ...) => partly done in derived class constructor
  //  - solver data loading => done in OneStepNS top class constructor.

  // rootNode == OneStepNSProblem
  // problemTypeNode == formalisation type (LCP ...)

  // === non smooth problem formalisation part ===
  problemTypeNode = SiconosDOMTreeTools::findNodeChild(rootNode);
  if (problemTypeNode == NULL)
    XMLException::selfThrow("OneStepNSProblemXML - constructor : tag NonSmooth Problem type not found");

  xmlNode *node;
  if ((node = SiconosDOMTreeTools::findNodeChild(problemTypeNode, "n")) != NULL)
    dimNode = node;

  // interactionConcerned
  if ((node = SiconosDOMTreeTools::findNodeChild(problemTypeNode, "Interaction_Concerned")) != NULL)
  {
    interactionConcernedNode = node;
    // Check if all interactions are concerned or not
    if (! hasAll())
    {
      // Get the indexList node
      if ((node = SiconosDOMTreeTools::findNodeChild(interactionConcernedNode, INDEX_LIST)) != NULL)
        interactionListNode = node;
      else
        XMLException::selfThrow("Tag indexList not found in Interaction.");
    }
  }
  // interaction list not required => all interactions are concerned by the problem.

  // === Solver part ===
  if ((node = SiconosDOMTreeTools::findNodeChild(problemTypeNode, "Solver")) != NULL)
  {
    solverNode = node;
    solverXML = new SolverXML(solverNode); // Is it really usefull to have a solverXML class?
    isSolverXMLAllocatedIn = true;
  }
  //  else -> solver = default one
}

OneStepNSProblemXML::~OneStepNSProblemXML()
{
  if (isSolverXMLAllocatedIn) delete solverXML;
  solverXML = NULL;
}

void OneStepNSProblemXML::setDimNSProblem(const int& n)
{
  if (! hasDim())
    dimNode = SiconosDOMTreeTools::createIntegerNode(problemTypeNode, "n", n);
  else SiconosDOMTreeTools::setIntegerContentValue(dimNode, n);
}

bool OneStepNSProblemXML::hasAll() const
{
  if (SiconosDOMTreeTools::hasAttributeValue(interactionConcernedNode, ALL_ATTRIBUTE))
    return SiconosDOMTreeTools::getAttributeValue<bool>(interactionConcernedNode, ALL_ATTRIBUTE);
  else return false;
}
void OneStepNSProblemXML::setAll(const bool& all)
{
  if (!hasAll())
  {
    if (all == true)
      xmlNewProp(interactionConcernedNode, (xmlChar*)ALL_ATTRIBUTE.c_str(), (xmlChar*)"true");
  }
  else
  {
    if (all == false)
      xmlRemoveProp(xmlHasProp(interactionConcernedNode, (xmlChar*)ALL_ATTRIBUTE.c_str()));
  }
}

void OneStepNSProblemXML::setSolverXMLPtr(SolverXML * solv)
{
  if (isSolverXMLAllocatedIn) delete solverXML;
  solverXML = solv;
  isSolverXMLAllocatedIn = false;
}

// \warning: following routine has to been checked and updated
void OneStepNSProblemXML::setSolver(const string& name, const string& normType,
                                    double tolerance, unsigned int maxIter, double searchDirection, double Rho)
{
  XMLException::selfThrow("OneStepNSProblemXML::setSolver, not yet implemented.");
  //   char tmpChar[128];
  //   xmlNode *node;
  //   xmlNodePtr solverAlgorithmNode;
  //   if( solverNode == NULL )
  //     {
  //       solverNode = xmlNewChild( rootNode, NULL, (xmlChar*)"Solver", NULL );
  //       node = xmlNewChild( solverNode, NULL, (xmlChar*)name.c_str(), NULL );
  //       solverAlgorithmNode = xmlNewChild( node, NULL, (xmlChar*)name.c_str(), NULL );

  //       sprintf(tmpChar, "%d", maxIter);
  //       if( maxIter != DEFAULT_ITER )
  //  xmlNewProp(solverAlgorithmNode, (xmlChar*)"maxIter", (xmlChar*)tmpChar );

  //       sprintf(tmpChar, "%f", tolerance);
  //       if( tolerance != DEFAULT_TOL )
  //  xmlNewProp(solverAlgorithmNode, (xmlChar*)"tolerance", (xmlChar*)tmpChar );

  //       /*
  //        *       /!\ normType not yet implemented in SICONOS/Numerics
  //        */
  //       if( normType != DEFAULT_NORMTYPE )
  //  xmlNewProp(solverAlgorithmNode, (xmlChar*)"normType", (xmlChar*)normType.c_str() );

  //       sprintf(tmpChar, "%f", searchDirection);
  //       if( searchDirection != DEFAULT_SEARCHDIR )
  //  xmlNewProp(solverAlgorithmNode, (xmlChar*)"searchDirection", (xmlChar*)tmpChar );
  //     }
  //   else
  //     {
  //       node = solverNode->next;
  //       if( node != NULL )
  //  node->name = (xmlChar*)name.c_str();
  //       else
  //  node = xmlNewChild( solverNode, NULL, (xmlChar*)name.c_str(), NULL );


  //       if( solverAlgorithmNode == NULL )
  //  {
  //    solverAlgorithmNode = xmlNewChild( node, NULL, (xmlChar*)name.c_str(), NULL );

  //    sprintf(tmpChar, "%d", maxIter);
  //    if( maxIter != DEFAULT_ITER )
  //      xmlNewProp(solverAlgorithmNode, (xmlChar*)"maxIter", (xmlChar*)tmpChar );

  //    sprintf(tmpChar, "%f", tolerance);
  //    if( tolerance != DEFAULT_TOL )
  //      xmlNewProp(solverAlgorithmNode, (xmlChar*)"tolerance", (xmlChar*)tmpChar );

  //    /*
  //     *       /!\ normType not yet implemented in SICONOS/Numerics
  //     */
  //    if( normType != DEFAULT_NORMTYPE )
  //      xmlNewProp(solverAlgorithmNode, (xmlChar*)"normType", (xmlChar*)normType.c_str() );

  //    sprintf(tmpChar, "%f", searchDirection);
  //    if( searchDirection != DEFAULT_SEARCHDIR )
  //      xmlNewProp(solverAlgorithmNode, (xmlChar*)"searchDirection", (xmlChar*)tmpChar );
  //  }
  //       else
  //  {
  //    xmlNode *node;
  //    node  = xmlNewNode( NULL, (xmlChar*)name.c_str() );

  //    sprintf(tmpChar, "%d", maxIter);
  //    if( maxIter != DEFAULT_ITER )
  //      xmlNewProp(solverAlgorithmNode, (xmlChar*)"maxIter", (xmlChar*)tmpChar );

  //    sprintf(tmpChar, "%f", tolerance);
  //    if( tolerance != DEFAULT_TOL )
  //      xmlNewProp(solverAlgorithmNode, (xmlChar*)"tolerance", (xmlChar*)tmpChar );

  //    /*
  //     *       /!\ normType not yet implemented in SICONOS/Numerics
  //     */
  //    if( normType != DEFAULT_NORMTYPE )
  //      xmlNewProp(solverAlgorithmNode, (xmlChar*)"normType", (xmlChar*)normType.c_str() );

  //    sprintf(tmpChar, "%f", searchDirection);
  //    if( searchDirection != DEFAULT_SEARCHDIR )
  //      xmlNewProp(solverAlgorithmNode, (xmlChar*)"searchDirection", (xmlChar*)tmpChar );

  //    xmlReplaceNode( solverAlgorithmNode, node );
  //  }
  //     }
}

void OneStepNSProblemXML::updateOneStepNSProblemXML(xmlNode* node, OneStepNSProblem* osnspb)
{
  rootNode = node;
  problemTypeNode = SiconosDOMTreeTools::findNodeChild(rootNode);
}
