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

