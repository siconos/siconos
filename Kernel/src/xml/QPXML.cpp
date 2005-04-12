#include "QPXML.h"

QPXML::QPXML() : OneStepNSProblemXML()
{
  this->QNode = NULL;
  this->pNode = NULL;
}

QPXML::QPXML(xmlNode * QPNode, vector<int> definedInteractionNumbers)
  : OneStepNSProblemXML(QPNode, definedInteractionNumbers)
{
  xmlNode *node, *qpModelNode;

  qpModelNode = SiconosDOMTreeTools::findNodeChild(QPNode);
  if (qpModelNode != NULL)
  {
    if (strcmp((char*)qpModelNode->name, OSNSP_SOLVER.c_str()) == 0)
      qpModelNode = SiconosDOMTreeTools::findFollowNode(qpModelNode);
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(qpModelNode, QP_Q)) != NULL)
  {
    this->QNode = node;
  }
  else
  {
    XMLException::selfThrow("QPXML - constructor : tag " + QP_Q + " not found.");
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(qpModelNode, QP_P)) != NULL)
  {
    this->pNode = node;
  }
  else
  {
    XMLException::selfThrow("QPXML - constructor : tag " + QP_P + " not found.");
  }

}

void QPXML::updateOneStepNSProblemXML(xmlNode* node, OneStepNSProblem* osnspb)
{
  IN("LCPXML::updateOneStepNSProblemXML\n");
  this->rootNode = node;
  this->rootNSProblemXMLNode = SiconosDOMTreeTools::findNodeChild(this->rootNode);
  //this->loadOneStepNSProblem( onsnpb );
  OUT("LCPXML::updateOneStepNSProblemXML\n");
}

