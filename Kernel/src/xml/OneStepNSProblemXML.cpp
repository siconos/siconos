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

#include "OneStepNSProblemXML.h"
using namespace std;


OneStepNSProblemXML::OneStepNSProblemXML():
  rootNSProblemXMLNode(NULL), rootNode(NULL), interactionConcernedNode(NULL),
  solverNode(NULL), solverAlgorithmNode(NULL), nNode(NULL)
{}

OneStepNSProblemXML::OneStepNSProblemXML(xmlNode * oneStepNSProblemXMLNode, vector<int> definedInteractionNumbers)
{
  xmlNode *node;

  rootNode = oneStepNSProblemXMLNode;
  //rootNSProblemXMLNode = oneStepNSProblemXMLNode;

  /*
   * node definition of the OneStepNSProblem
   */
  node = SiconosDOMTreeTools::findNodeChild(rootNode);
  if (node != NULL)
  {
    if (strcmp((char*)node->name, OSNSP_SOLVER.c_str()) == 0)
    {
      /*
       * in this case, we got the Solver tag in first
       */
      solverNode = node;
      node = SiconosDOMTreeTools::findNodeChild(node);
      solverAlgorithmNode = SiconosDOMTreeTools::findNodeChild(node);

      rootNSProblemXMLNode = SiconosDOMTreeTools::findFollowNode(solverNode);
    }
    else
    {
      /*
       * in that case, it must be solving model definition in first
       */
      rootNSProblemXMLNode = SiconosDOMTreeTools::findNodeChild(rootNode);

      if ((node = SiconosDOMTreeTools::findNodeChild(rootNode, OSNSP_SOLVER)) != NULL)
      {
        solverNode = node;
        node = SiconosDOMTreeTools::findNodeChild(node);
        solverAlgorithmNode = SiconosDOMTreeTools::findNodeChild(node);
      }
      else
      {
        solverNode = NULL;
        solverAlgorithmNode = NULL;
        cout << "Warning : optional tag not found : OneStepNSProblemXML - constructor, tag " << OSNSP_SOLVER << " not found." << endl;
      }
    }
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNSProblemXMLNode, OSNSP_N)) != NULL)
    nNode = node;
  else
  {
    //XMLException::selfThrow("OneStepNSProblemXML - constructor : tag " + OSNSP_N + " not found.");
    nNode = NULL;
    cout << "Warning : optional tag not found : OneStepNSProblemXML - constructor : tag " << OSNSP_N << " not found." << endl;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(rootNSProblemXMLNode, OSNSP_INTERACTION_CONCERNED)) != NULL)
  {
    interactionConcernedNode = node;
    loadOneStepNSProblemConcernedInteraction(node, definedInteractionNumbers);
  }
  else
    XMLException::selfThrow("OneStepNSProblemXML - constructor : tag " + OSNSP_INTERACTION_CONCERNED + " not found.");

  //  if ((node=SiconosDOMTreeTools::findNodeChild(rootNode, OSNSP_SOLVER)) !=NULL)
  //    {
  //      cout<<" node Solver found !!!!!!!!!!"<<endl<<"<<press enter>>"<<endl;
  //      getchar();
  //      solverNode = node;
  //      node = SiconosDOMTreeTools::findNodeChild( node );
  //      solverAlgorithmNode = SiconosDOMTreeTools::findNodeChild( node );
  //    }
  //    else
  //      cout<<"Warning : optional tag not found : OneStepNSProblemXML - constructor, tag "<<OSNSP_SOLVER<<" not found."<<endl;
  //    //XMLException::selfThrow("OneStepNSProblemXML - constructor : tag " + OSNSP_SOLVER + " not found.");
}

OneStepNSProblemXML::~OneStepNSProblemXML()
{}


void OneStepNSProblemXML::loadOneStepNSProblemConcernedInteraction(xmlNode * interactionConcernedNode, vector<int> definedInteractionNumbers)
{
  xmlNode *interactionNode;
  int number;
  int size = 0, i = 0;
  unsigned int j = 0;

  interactionNumbersVector.clear();

  if (((interactionNode = SiconosDOMTreeTools::findNodeChild((const xmlNode*)interactionConcernedNode, INTERACTION_TAG)) == NULL) && (!SiconosDOMTreeTools::getBooleanAttributeValue(interactionConcernedNode, ALL_ATTRIBUTE)))
    XMLException::selfThrow("OneStepNSProblemXML - loadOneStepNSProblemConcernedInteraction error : at least one " + INTERACTION_TAG + " tag must be declared in " + OSNSP_INTERACTION_CONCERNED + " tag.");

  size = SiconosDOMTreeTools::getNodeChildrenNumber(interactionConcernedNode);
  //  cout<<"number of children = "<<size<<endl;
  //  cout<<"# Press <<ENTER>> !!"<<endl;
  //  getchar();

  while ((interactionNode != NULL) && (i < size))
  {
    number = SiconosDOMTreeTools::getIntegerAttributeValue(interactionNode, NUMBER_ATTRIBUTE);

    //Verifying Interaction number exists
    j = 0;
    while ((j < definedInteractionNumbers.size()) && (number != definedInteractionNumbers[j]))
      j++;

    if (j == definedInteractionNumbers.size())
    {
      char errorMsg[1024];
      sprintf(errorMsg, "OneStepNSProblemXML - loadOneStepNSProblemConcernedInteraction error : in a tag %s you define an Interaction number who doesn't exist : %d.", OSNSP_INTERACTION_CONCERNED.c_str(), number);
      XMLException::selfThrow(errorMsg);
    }

    interactionNumbersVector.push_back(number);

    interactionNode = SiconosDOMTreeTools::findFollowNode(interactionNode, INTERACTION_TAG);
    i++;
  }

  if (i < size)
  {
    XMLException::selfThrow("OneStepIntegratorXML - loadOneStepNSProblemConcernedInteraction error : the size attribute given in the tag " + OSNSP_INTERACTION_CONCERNED + " is not correct.");
  }
}

void OneStepNSProblemXML::setInteractionConcerned(vector<int> v, bool all)
{
  unsigned int i;
  xmlNode* node;
  xmlNode* interactionNode;
  char num[32];

  if (all)
  {
    /*
     * in that case, the attribute "all" of the Interaction_Concerned must be used
     */
    node = xmlNewNode(NULL, (xmlChar*)OSNSP_INTERACTION_CONCERNED.c_str());
    xmlSetProp(node, (xmlChar*)ALL_ATTRIBUTE.c_str(), (xmlChar*)"true");

    if (interactionConcernedNode != NULL)
      xmlReplaceNode(interactionConcernedNode, node);
    else xmlAddChild(rootNSProblemXMLNode, node);
    interactionConcernedNode = node;
    cout << "#### OneStepNSProblemXML::setInteractionConcerned ALL" << endl;
  }
  else
  {
    cout << "#### OneStepNSProblemXML::setInteractionConcerned !ALL" << endl;
    if (interactionConcernedNode == NULL)
    {
      node = xmlNewChild(rootNSProblemXMLNode, NULL, (xmlChar*)OSNSP_INTERACTION_CONCERNED.c_str(), NULL);
      interactionConcernedNode = node;
      //sprintf(num, "%d", v.size());
      //xmlNewProp(node, (xmlChar*)OSNSP_SIZE.c_str(), (xmlChar*)num );

      for (i = 0; i < v.size(); i++)
      {
        interactionNode = xmlNewChild(node, NULL, (xmlChar*)INTERACTION_TAG.c_str(), NULL);
        sprintf(num, "%d", v[i]);
        xmlNewProp(interactionNode, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
      }
      interactionNumbersVector = v;
    }
    else
    {
      /*
       * it must check if a nmber of the vector is not already in the DOM tree
       */
      vector<int>::iterator iter;
      for (i = 0; i < v.size(); i++)
      {
        iter = find(interactionNumbersVector.begin(), interactionNumbersVector.end(), v[i]);
        if (iter == interactionNumbersVector.end())
        {
          /*
           * in this case, we must add in the DOMtree the number of this Interaction
           */
          interactionNode = xmlNewChild(interactionConcernedNode, NULL, (xmlChar*)INTERACTION_TAG.c_str(), NULL);
          sprintf(num, "%d", v[i]);
          xmlNewProp(interactionNode, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
        }
      }
    }
  }
}

void OneStepNSProblemXML::setSolver(string name, string methodName, string normType,
                                    double tolerance, int maxIter, double searchDirection)
{
  char tmpChar[128];
  xmlNode *node;

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

//void OneStepNSProblemXML::loadOneStepNSProblem( /*OneStepNSProblem *osnsp*/ )
//{
//  /*
//   * At this time, the OneStepNSProblemXML is just created and contains only
//   * one node for the solving model (LCP, QP, Relay, ...)
//   *
//   * the rootNSProblemXMLNode is set to the value of the node containing the solving model
//   */
//  rootNSProblemXMLNode = SiconosDOMTreeTools::findNodeChild( rootNode);
//}

