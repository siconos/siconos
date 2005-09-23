#include "DFC_2DXML.h"
using namespace std;

DFC_2DXML::DFC_2DXML() : OneStepNSProblemXML()
{
  this->MNode = NULL;
  this->qNode = NULL;
}

DFC_2DXML::DFC_2DXML(xmlNode * DFC_2DNode, vector<int> definedInteractionNumbers)
  : OneStepNSProblemXML(DFC_2DNode, definedInteractionNumbers)
{
  xmlNode *node, *dfc_2DModelNode;

  dfc_2DModelNode = SiconosDOMTreeTools::findNodeChild(DFC_2DNode);
  if (dfc_2DModelNode != NULL)
  {
    if (strcmp((char*)dfc_2DModelNode->name, OSNSP_SOLVER.c_str()) == 0)
      dfc_2DModelNode = SiconosDOMTreeTools::findFollowNode(dfc_2DModelNode);
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(dfc_2DModelNode, DFC_2D_M)) != NULL)
  {
    this->MNode = node;
  }
  else
  {
    //XMLException::selfThrow("DFC_2DXML - constructor : tag " + DFC_2D_M + " not found.");
    this->MNode = NULL;
    cout << "DFC_2DXML - constructor : tag " << DFC_2D_M << " not found. Optional attribute." << endl;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(dfc_2DModelNode, DFC_2D_Q)) != NULL)
  {
    this->qNode = node;
  }
  else
  {
    //XMLException::selfThrow("DFC_2DXML - constructor : tag " + DFC_2D_Q + " not found.");
    this->qNode = NULL;
    cout << "DFC_2DXML - constructor : tag " << DFC_2D_Q << " not found. Optional attribute." << endl;
  }
}

void DFC_2DXML::updateOneStepNSProblemXML(xmlNode* node, OneStepNSProblem* osnspb)
{
  IN("DFC_2DXML::updateOneStepNSProblemXML\n");
  this->rootNode = node;
  this->rootNSProblemXMLNode = SiconosDOMTreeTools::findNodeChild(this->rootNode);
  //this->loadOneStepNSProblem( osnspb );
  OUT("DFC_2DXML::updateOneStepNSProblemXML\n");
}

