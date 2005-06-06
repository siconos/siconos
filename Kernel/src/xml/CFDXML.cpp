#include "CFDXML.h"
using namespace std;

CFDXML::CFDXML() : OneStepNSProblemXML()
{
  this->MNode = NULL;
  this->qNode = NULL;
}

CFDXML::CFDXML(xmlNode * CFDNode, vector<int> definedInteractionNumbers)
  : OneStepNSProblemXML(CFDNode, definedInteractionNumbers)
{
  xmlNode *node, *cfdModelNode;

  cfdModelNode = SiconosDOMTreeTools::findNodeChild(CFDNode);
  if (cfdModelNode != NULL)
  {
    if (strcmp((char*)cfdModelNode->name, OSNSP_SOLVER.c_str()) == 0)
      cfdModelNode = SiconosDOMTreeTools::findFollowNode(cfdModelNode);
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(cfdModelNode, CFD_M)) != NULL)
  {
    this->MNode = node;
  }
  else
  {
    //XMLException::selfThrow("CFDXML - constructor : tag " + CFD_M + " not found.");
    this->MNode = NULL;
    cout << "CFDXML - constructor : tag " << CFD_M << " not found. Optional attribute." << endl;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(cfdModelNode, CFD_Q)) != NULL)
  {
    this->qNode = node;
  }
  else
  {
    //XMLException::selfThrow("CFDXML - constructor : tag " + CFD_Q + " not found.");
    this->qNode = NULL;
    cout << "CFDXML - constructor : tag " << CFD_Q << " not found. Optional attribute." << endl;
  }
}

void CFDXML::updateOneStepNSProblemXML(xmlNode* node, OneStepNSProblem* osnspb)
{
  IN("CFDXML::updateOneStepNSProblemXML\n");
  this->rootNode = node;
  this->rootNSProblemXMLNode = SiconosDOMTreeTools::findNodeChild(this->rootNode);
  //this->loadOneStepNSProblem( osnspb );
  OUT("CFDXML::updateOneStepNSProblemXML\n");
}

