
#include "OneStepNSProblemXML.h"

#include "check.h"

OneStepNSProblemXML::OneStepNSProblemXML()
{
  this->rootNSProblemXMLNode = NULL;
  this->nNode = NULL;
  this->interactionConcernedNode = NULL;
  this->solverNode = NULL;
  this->solverAlgorithmNode = NULL;
}

OneStepNSProblemXML::OneStepNSProblemXML(xmlNode * oneStepNSProblemXMLNode, vector<int> definedInteractionNumbers)
{
  xmlNode *node;

  this->rootNode = oneStepNSProblemXMLNode;
  //this->rootNSProblemXMLNode = oneStepNSProblemXMLNode;

  /*
   * node definition of the OneStepNSProblem
   */
  node = SiconosDOMTreeTools::findNodeChild(this->rootNode);
  if (node != NULL)
  {
    if (strcmp((char*)node->name, OSNSP_SOLVER.c_str()) == 0)
    {
      /*
       * in this case, we got the Solver tag in first
       */
      this->solverNode = node;
      node = SiconosDOMTreeTools::findNodeChild(node);
      this->solverAlgorithmNode = SiconosDOMTreeTools::findNodeChild(node);

      this->rootNSProblemXMLNode = SiconosDOMTreeTools::findFollowNode(this->solverNode);
    }
    else
    {
      /*
       * in that case, it must be solving model definition in first
       */
      this->rootNSProblemXMLNode = SiconosDOMTreeTools::findNodeChild(this->rootNode);

      if ((node = SiconosDOMTreeTools::findNodeChild(this->rootNode, OSNSP_SOLVER)) != NULL)
      {
        this->solverNode = node;
        node = SiconosDOMTreeTools::findNodeChild(node);
        this->solverAlgorithmNode = SiconosDOMTreeTools::findNodeChild(node);
      }
      else
      {
        this->solverNode = NULL;
        this->solverAlgorithmNode = NULL;
        cout << "Warning : optional tag not found : OneStepNSProblemXML - constructor, tag " << OSNSP_SOLVER << " not found." << endl;
      }
    }
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootNSProblemXMLNode, OSNSP_N)) != NULL)
    this->nNode = node;
  else
  {
    //XMLException::selfThrow("OneStepNSProblemXML - constructor : tag " + OSNSP_N + " not found.");
    this->nNode = NULL;
    cout << "Warning : optional tag not found : OneStepNSProblemXML - constructor : tag " << OSNSP_N << " not found." << endl;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootNSProblemXMLNode, OSNSP_INTERACTION_CONCERNED)) != NULL)
  {
    this->interactionConcernedNode = node;
    this->loadOneStepNSProblemConcernedInteraction(node, definedInteractionNumbers);
  }
  else
    XMLException::selfThrow("OneStepNSProblemXML - constructor : tag " + OSNSP_INTERACTION_CONCERNED + " not found.");

  //  if ((node=SiconosDOMTreeTools::findNodeChild(this->rootNode, OSNSP_SOLVER)) !=NULL)
  //    {
  //      cout<<" node Solver found !!!!!!!!!!"<<endl<<"<<press enter>>"<<endl;
  //      getchar();
  //      this->solverNode = node;
  //      node = SiconosDOMTreeTools::findNodeChild( node );
  //      this->solverAlgorithmNode = SiconosDOMTreeTools::findNodeChild( node );
  //    }
  //    else
  //      cout<<"Warning : optional tag not found : OneStepNSProblemXML - constructor, tag "<<OSNSP_SOLVER<<" not found."<<endl;
  //    //XMLException::selfThrow("OneStepNSProblemXML - constructor : tag " + OSNSP_SOLVER + " not found.");
}

void OneStepNSProblemXML::loadOneStepNSProblemConcernedInteraction(xmlNode * interactionConcernedNode, vector<int> definedInteractionNumbers)
{
  xmlNode *interactionNode;
  int number;
  int size = 0;
  int i = 0, j = 0;

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

    this->interactionNumbersVector.push_back(number);

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
  int i;
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

    if (this->interactionConcernedNode != NULL)
      xmlReplaceNode(this->interactionConcernedNode, node);
    else xmlAddChild(this->rootNSProblemXMLNode, node);
    this->interactionConcernedNode = node;
    cout << "#### OneStepNSProblemXML::setInteractionConcerned ALL" << endl;
  }
  else
  {
    cout << "#### OneStepNSProblemXML::setInteractionConcerned !ALL" << endl;
    if (this->interactionConcernedNode == NULL)
    {
      node = xmlNewChild(this->rootNSProblemXMLNode, NULL, (xmlChar*)OSNSP_INTERACTION_CONCERNED.c_str(), NULL);
      this->interactionConcernedNode = node;
      //sprintf(num, "%d", v.size());
      //xmlNewProp(node, (xmlChar*)OSNSP_SIZE.c_str(), (xmlChar*)num );

      for (i = 0; i < v.size(); i++)
      {
        interactionNode = xmlNewChild(node, NULL, (xmlChar*)INTERACTION_TAG.c_str(), NULL);
        sprintf(num, "%d", v[i]);
        xmlNewProp(interactionNode, (xmlChar*)NUMBER_ATTRIBUTE.c_str(), (xmlChar*)num);
      }
      this->interactionNumbersVector = v;
    }
    else
    {
      /*
       * it must check if a nmber of the vector is not already in the DOM tree
       */
      vector<int>::iterator iter;
      for (i = 0; i < v.size(); i++)
      {
        iter = find(this->interactionNumbersVector.begin(), this->interactionNumbersVector.end(), v[i]);
        if (iter == this->interactionNumbersVector.end())
        {
          /*
           * in this case, we must add in the DOMtree the number of this Interaction
           */
          interactionNode = xmlNewChild(this->interactionConcernedNode, NULL, (xmlChar*)INTERACTION_TAG.c_str(), NULL);
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
  int tmpInt;
  double tmpDouble;
  char tmpChar[128];
  xmlNode *node;

  if (this->solverNode == NULL)
  {
    this->solverNode = xmlNewChild(this->rootNode, NULL, (xmlChar*)OSNSP_SOLVER.c_str(), NULL);
    node = xmlNewChild(this->solverNode, NULL, (xmlChar*)name.c_str(), NULL);
    this->solverAlgorithmNode = xmlNewChild(node, NULL, (xmlChar*)methodName.c_str(), NULL);

    sprintf(tmpChar, "%d", maxIter);
    if (maxIter != DefaultAlgoMaxIter)
      xmlNewProp(this->solverAlgorithmNode, (xmlChar*)OSNSP_MAXITER.c_str(), (xmlChar*)tmpChar);

    sprintf(tmpChar, "%f", tolerance);
    if (tolerance != DefaultAlgoTolerance)
      xmlNewProp(this->solverAlgorithmNode, (xmlChar*)OSNSP_TOLERANCE.c_str(), (xmlChar*)tmpChar);

    /*
     *       /!\ normType not yet implemented in SICONOS/Numerics
     */
    if (normType != DefaultAlgoNormType)
      xmlNewProp(this->solverAlgorithmNode, (xmlChar*)OSNSP_NORMTYPE.c_str(), (xmlChar*)normType.c_str());

    sprintf(tmpChar, "%f", searchDirection);
    if (searchDirection != DefaultAlgoSearchDirection)
      xmlNewProp(this->solverAlgorithmNode, (xmlChar*)OSNSP_SEARCHDIRECTION.c_str(), (xmlChar*)tmpChar);
  }
  else
  {
    node = this->solverNode->next;
    if (node != NULL)
      node->name = (xmlChar*)name.c_str();
    else
      node = xmlNewChild(this->solverNode, NULL, (xmlChar*)name.c_str(), NULL);


    if (this->solverAlgorithmNode == NULL)
    {
      this->solverAlgorithmNode = xmlNewChild(node, NULL, (xmlChar*)methodName.c_str(), NULL);

      sprintf(tmpChar, "%d", maxIter);
      if (maxIter != DefaultAlgoMaxIter)
        xmlNewProp(this->solverAlgorithmNode, (xmlChar*)OSNSP_MAXITER.c_str(), (xmlChar*)tmpChar);

      sprintf(tmpChar, "%f", tolerance);
      if (tolerance != DefaultAlgoTolerance)
        xmlNewProp(this->solverAlgorithmNode, (xmlChar*)OSNSP_TOLERANCE.c_str(), (xmlChar*)tmpChar);

      /*
       *       /!\ normType not yet implemented in SICONOS/Numerics
       */
      if (normType != DefaultAlgoNormType)
        xmlNewProp(this->solverAlgorithmNode, (xmlChar*)OSNSP_NORMTYPE.c_str(), (xmlChar*)normType.c_str());

      sprintf(tmpChar, "%f", searchDirection);
      if (searchDirection != DefaultAlgoSearchDirection)
        xmlNewProp(this->solverAlgorithmNode, (xmlChar*)OSNSP_SEARCHDIRECTION.c_str(), (xmlChar*)tmpChar);
    }
    else
    {
      xmlNode *node;
      node  = xmlNewNode(NULL, (xmlChar*)methodName.c_str());

      sprintf(tmpChar, "%d", maxIter);
      if (maxIter != DefaultAlgoMaxIter)
        xmlNewProp(this->solverAlgorithmNode, (xmlChar*)OSNSP_MAXITER.c_str(), (xmlChar*)tmpChar);

      sprintf(tmpChar, "%f", tolerance);
      if (tolerance != DefaultAlgoTolerance)
        xmlNewProp(this->solverAlgorithmNode, (xmlChar*)OSNSP_TOLERANCE.c_str(), (xmlChar*)tmpChar);

      /*
       *       /!\ normType not yet implemented in SICONOS/Numerics
       */
      if (normType != DefaultAlgoNormType)
        xmlNewProp(this->solverAlgorithmNode, (xmlChar*)OSNSP_NORMTYPE.c_str(), (xmlChar*)normType.c_str());

      sprintf(tmpChar, "%f", searchDirection);
      if (searchDirection != DefaultAlgoSearchDirection)
        xmlNewProp(this->solverAlgorithmNode, (xmlChar*)OSNSP_SEARCHDIRECTION.c_str(), (xmlChar*)tmpChar);

      xmlReplaceNode(this->solverAlgorithmNode, node);
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
//  this->rootNSProblemXMLNode = SiconosDOMTreeTools::findNodeChild( this->rootNode);
//}

