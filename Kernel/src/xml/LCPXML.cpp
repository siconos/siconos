//$Id: LCPXML.cpp,v 1.10 2005/02/04 14:52:44 jbarbier Exp $
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

//$Log: LCPXML.cpp,v $
//Revision 1.10  2005/02/04 14:52:44  jbarbier
//- Rolling balls in progress (contact is detected)
//
//- time data is given in parameter for computInput and Output in the Relation. Parameter is sent by methods of the OneStepNSProblem
//
//Revision 1.9  2005/01/25 09:27:18  jbarbier
//- save of Solver tag in the OneStepNSProblem tag available when saving without XML input file and with partial XML input file
//
//Revision 1.8  2005/01/24 14:33:03  jbarbier
//- OneStepNSProblem > Solver tag is available and managed in the XML part
//
//- tests added on OneStepNSProblem > Solver tag
//
//Revision 1.7  2004/09/14 13:49:56  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.6  2004/06/29 15:12:02  acary
//Change in the naming comvention for the LCP
//The LCP Matrix is now denoted by M.
//The LCP Vector is now denoted by q.
//