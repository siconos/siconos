#include "LCPXML.h"

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
  {
    //XMLException::selfThrow("LCPXML - constructor : tag " + LCP_M + " not found.");
    this->MNode = NULL;
    cout << "LCPXML - constructor : tag " << LCP_M << " not found. Optional attribute." << endl;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(lcpModelNode, LCP_Q)) != NULL)
  {
    this->qNode = node;
  }
  else
  {
    //XMLException::selfThrow("LCPXML - constructor : tag " + LCP_Q + " not found.");
    this->qNode = NULL;
    cout << "LCPXML - constructor : tag " << LCP_Q << " not found. Optional attribute." << endl;
  }
}

void LCPXML::updateOneStepNSProblemXML(xmlNode* node, OneStepNSProblem* osnspb)
{
  IN("LCPXML::updateOneStepNSProblemXML\n");
  this->rootNode = node;
  this->rootNSProblemXMLNode = SiconosDOMTreeTools::findNodeChild(this->rootNode);
  //this->loadOneStepNSProblem( osnspb );
  OUT("LCPXML::updateOneStepNSProblemXML\n");
}

