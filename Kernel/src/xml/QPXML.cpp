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

//$Log: QPXML.cpp,v $
//Revision 1.11  2005/01/25 09:27:18  jbarbier
//- save of Solver tag in the OneStepNSProblem tag available when saving without XML input file and with partial XML input file
//
//Revision 1.10  2005/01/24 14:33:03  jbarbier
//- OneStepNSProblem > Solver tag is available and managed in the XML part
//
//- tests added on OneStepNSProblem > Solver tag
//
//Revision 1.9  2004/09/27 13:27:14  jbarbier
//
//- Siconos schema renamed : SiconosModelSchema-V1.0.xsd
//
//- new required tags of the model : title, author, description, date, xmlSchema.
//They replace previous attributes author, description and date of the Model.
//
//Revision 1.8  2004/09/14 13:49:59  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.7  2004/07/29 14:25:45  jbarbier
