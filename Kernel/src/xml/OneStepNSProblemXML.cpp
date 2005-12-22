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
#include "NSQPSolverXML.h"
 */

#include "OneStepNSProblemXML.h"
// To be removed thanks to factories:
#include "LemkeSolverXML.h"
#include "LexicoLemkeSolverXML.h"
#include "QPSolverXML.h"
#include "NSQPSolverXML.h"
#include "NLGSSolverXML.h"
#include "CPGSolverXML.h"
#include "LatinSolverXML.h"
using namespace std;


OneStepNSProblemXML::OneStepNSProblemXML():
  rootNode(NULL), problemTypeNode(NULL), dimNode(NULL), interactionConcernedNode(NULL),
  interactionListNode(NULL), solverNode(NULL), solverXML(NULL), isSolverXMLAllocatedIn(false)
{}

OneStepNSProblemXML::OneStepNSProblemXML(xmlNode * oneStepNSProblemXMLNode, vector<int> definedInteractionNumbers):
  rootNode(oneStepNSProblemXMLNode), problemTypeNode(NULL), dimNode(NULL), interactionConcernedNode(NULL),
  interactionListNode(NULL), solverNode(NULL), solverXML(NULL), isSolverXMLAllocatedIn(false)
{
  // Two steps in OneStepNSProblem xml loading:
  //  - problem formalisation part ( LCP, FrictionContact ...) => partly done in derived class constructor
  //  - solver data loading => done in OneStepNS top class constructor.

  // rootNode == OneStepNSProblem
  // problemTypeNode == formalisation type (LCP ...)

  // === Solver part ===
  solverNode = SiconosDOMTreeTools::findNodeChild(rootNode, "Solver");

  if (solverNode != NULL)
  {
    // \todo To be updated thanks to factories to avoid if ... else if ...
    xmlNodePtr solvingFormalisationNode = SiconosDOMTreeTools::findNodeChild(solverNode);
    xmlNodePtr solverAlgorithmNode = SiconosDOMTreeTools::findNodeChild(solvingFormalisationNode);
    if (solvingFormalisationNode == NULL || solverAlgorithmNode == NULL)
      XMLException::selfThrow("OneStepNSProblemXML - constructor : algorithm node name or solving form node name not found");

    string type = (char*)solverAlgorithmNode->name;
    if (type == "Lemke")
      solverXML = new LemkeSolverXML(solverNode, solvingFormalisationNode, solverAlgorithmNode);
    else if (type == "LexicoLemke")
      solverXML = new LexicoLemkeSolverXML(solverNode, solvingFormalisationNode, solverAlgorithmNode);
    else if (type == "QP")
      solverXML = new QPSolverXML(solverNode, solvingFormalisationNode, solverAlgorithmNode);
    else if (type == "NSQP")
      solverXML = new NSQPSolverXML(solverNode, solvingFormalisationNode, solverAlgorithmNode);
    else if (type == "NLGS")
      solverXML = new NLGSSolverXML(solverNode, solvingFormalisationNode, solverAlgorithmNode);
    else if (type == "CPG")
      solverXML = new CPGSolverXML(solverNode, solvingFormalisationNode, solverAlgorithmNode);
    else if (type == "Latin")
      solverXML = new LatinSolverXML(solverNode, solvingFormalisationNode, solverAlgorithmNode);
    else
      XMLException::selfThrow("OneStepNSProblemXML constructor, undefined solver algorithm type: " + type);
    isSolverXMLAllocatedIn = true;
  }
  else
    XMLException::selfThrow("OneStepNSProblemXML - constructor : tag Solver not found.");

  // === non smooth problem formalisation part ===
  problemTypeNode = SiconosDOMTreeTools::findNodeChild(rootNode);
  if (problemTypeNode == NULL)
    XMLException::selfThrow("OneStepNSProblemXML - constructor : tag NonSmooth Problem type not found");

  xmlNode *node;
  if ((node = SiconosDOMTreeTools::findNodeChild(problemTypeNode, "n")) != NULL)
    dimNode = node;

  // interactionConcerned
  if ((node = SiconosDOMTreeTools::findNodeChild(problemTypeNode, OSNSP_INTERACTION_CONCERNED)) != NULL)
  {
    interactionConcernedNode = node;
    // Check if all interactions are concerned or not
    if (! hasAll())
    {
      // Get the indexList node
      if ((node = SiconosDOMTreeTools::findNodeChild(interactionConcernedNode, INDEX_LIST)) != NULL)
        interactionListNode = node;
      else
        XMLException::selfThrow("tag indexList not found.");
    }
  }
  // interaction list not required => all interactions are concerned by the problem.
  // else
  //  XMLException::selfThrow("OneStepNSProblemXML - constructor : tag " + OSNSP_INTERACTION_CONCERNED + " not found.");
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
void OneStepNSProblemXML::setSolver(const string& name, const string& methodName, const string& normType,
                                    const double& tolerance, const unsigned int& maxIter, const double& searchDirection)
{
  char tmpChar[128];
  xmlNode *node;
  xmlNodePtr solverAlgorithmNode;
  if (solverNode == NULL)
  {
    solverNode = xmlNewChild(rootNode, NULL, (xmlChar*)OSNSP_SOLVER.c_str(), NULL);
    node = xmlNewChild(solverNode, NULL, (xmlChar*)name.c_str(), NULL);
    solverAlgorithmNode = xmlNewChild(node, NULL, (xmlChar*)methodName.c_str(), NULL);

    sprintf(tmpChar, "%d", maxIter);
    if (maxIter != DefaultAlgoMaxIter)
      xmlNewProp(solverAlgorithmNode, (xmlChar*)OSNSP_MAXITER.c_str(), (xmlChar*)tmpChar);

    sprintf(tmpChar, "%f", tolerance);
    if (tolerance != DefaultAlgoTolerance)
      xmlNewProp(solverAlgorithmNode, (xmlChar*)OSNSP_TOLERANCE.c_str(), (xmlChar*)tmpChar);

    /*
     *       /!\ normType not yet implemented in SICONOS/Numerics
     */
    if (normType != DefaultAlgoNormType)
      xmlNewProp(solverAlgorithmNode, (xmlChar*)OSNSP_NORMTYPE.c_str(), (xmlChar*)normType.c_str());

    sprintf(tmpChar, "%f", searchDirection);
    if (searchDirection != DefaultAlgoSearchDirection)
      xmlNewProp(solverAlgorithmNode, (xmlChar*)OSNSP_SEARCHDIRECTION.c_str(), (xmlChar*)tmpChar);
  }
  else
  {
    node = solverNode->next;
    if (node != NULL)
      node->name = (xmlChar*)name.c_str();
    else
      node = xmlNewChild(solverNode, NULL, (xmlChar*)name.c_str(), NULL);


    if (solverAlgorithmNode == NULL)
    {
      solverAlgorithmNode = xmlNewChild(node, NULL, (xmlChar*)methodName.c_str(), NULL);

      sprintf(tmpChar, "%d", maxIter);
      if (maxIter != DefaultAlgoMaxIter)
        xmlNewProp(solverAlgorithmNode, (xmlChar*)OSNSP_MAXITER.c_str(), (xmlChar*)tmpChar);

      sprintf(tmpChar, "%f", tolerance);
      if (tolerance != DefaultAlgoTolerance)
        xmlNewProp(solverAlgorithmNode, (xmlChar*)OSNSP_TOLERANCE.c_str(), (xmlChar*)tmpChar);

      /*
       *       /!\ normType not yet implemented in SICONOS/Numerics
       */
      if (normType != DefaultAlgoNormType)
        xmlNewProp(solverAlgorithmNode, (xmlChar*)OSNSP_NORMTYPE.c_str(), (xmlChar*)normType.c_str());

      sprintf(tmpChar, "%f", searchDirection);
      if (searchDirection != DefaultAlgoSearchDirection)
        xmlNewProp(solverAlgorithmNode, (xmlChar*)OSNSP_SEARCHDIRECTION.c_str(), (xmlChar*)tmpChar);
    }
    else
    {
      xmlNode *node;
      node  = xmlNewNode(NULL, (xmlChar*)methodName.c_str());

      sprintf(tmpChar, "%d", maxIter);
      if (maxIter != DefaultAlgoMaxIter)
        xmlNewProp(solverAlgorithmNode, (xmlChar*)OSNSP_MAXITER.c_str(), (xmlChar*)tmpChar);

      sprintf(tmpChar, "%f", tolerance);
      if (tolerance != DefaultAlgoTolerance)
        xmlNewProp(solverAlgorithmNode, (xmlChar*)OSNSP_TOLERANCE.c_str(), (xmlChar*)tmpChar);

      /*
       *       /!\ normType not yet implemented in SICONOS/Numerics
       */
      if (normType != DefaultAlgoNormType)
        xmlNewProp(solverAlgorithmNode, (xmlChar*)OSNSP_NORMTYPE.c_str(), (xmlChar*)normType.c_str());

      sprintf(tmpChar, "%f", searchDirection);
      if (searchDirection != DefaultAlgoSearchDirection)
        xmlNewProp(solverAlgorithmNode, (xmlChar*)OSNSP_SEARCHDIRECTION.c_str(), (xmlChar*)tmpChar);

      xmlReplaceNode(solverAlgorithmNode, node);
    }
  }
}

void OneStepNSProblemXML::updateOneStepNSProblemXML(xmlNode* node, OneStepNSProblem* osnspb)
{
  rootNode = node;
  problemTypeNode = SiconosDOMTreeTools::findNodeChild(rootNode);
}
