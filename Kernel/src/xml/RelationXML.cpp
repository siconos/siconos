
#include "RelationXML.h"


RelationXML::RelationXML()
{
  this->computeInputNode = NULL;
  this->computeOutputNode = NULL;
}

RelationXML::RelationXML(xmlNode *relationNode)
{
  IN("RelationXML::RelationXML(xmlNode*)\n");
  xmlNode *node;
  string type((char*)relationNode->name);
  this->rootRelationXMLNode = relationNode;

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootRelationXMLNode, COMPUTE_INPUT_TAG)) != NULL)
  {
    this->computeInputNode = node;
  }
  else
  {
    //if( SiconosDOMTreeTools::getStringAttributeValue(this->rootRelationXMLNode, RELATION_TYPE) == RELATION_LNL )
    if (LAGRANGIAN_NON_LINEAR_RELATION_TAG == type)
    {
      XMLException::selfThrow("RelationXML - RelationXML::RelationXML(xmlNode *relationNode) error : tag " + COMPUTE_INPUT_TAG + " not found.");
    }
    //else cout<< "RelationXML - RelationXML::RelationXML(xmlNode *relationNode) Warning : tag " << RELATION_INPUT << " not found. This is attribute is optional for this relation" <<endl;
    this->computeInputNode = NULL;
  }

  if ((node = SiconosDOMTreeTools::findNodeChild(this->rootRelationXMLNode, COMPUTE_OUTPUT_TAG)) != NULL)
  {
    this->computeOutputNode = node;
  }
  else
  {
    //if( SiconosDOMTreeTools::getStringAttributeValue(this->rootRelationXMLNode, RELATION_TYPE) == RELATION_LNL )
    if (LAGRANGIAN_NON_LINEAR_RELATION_TAG == type)
    {
      XMLException::selfThrow("RelationXML - RelationXML::RelationXML(xmlNode *relationNode) error : tag " + COMPUTE_OUTPUT_TAG + " not found.");
    }
    //else cout<< "RelationXML - RelationXML::RelationXML(xmlNode *relationNode) Warning : tag " << RELATION_OUTPUT << " not found. This is attribute is optional for this relation" <<endl;
    this->computeOutputNode = NULL;
  }
  OUT("RelationXML::RelationXML(xmlNode*)\n");
}

RelationXML::~RelationXML()
{}

void RelationXML::updateRelationXML(xmlNode* node, Relation* rel)
{
  IN("RelationXML::updateRelationXML\n");
  this->rootRelationXMLNode = node;
  OUT("RelationXML::updateRelationXML\n");
}

//$Log: RelationXML.cpp,v $
//Revision 1.13  2005/03/07 13:17:21  jbarbier
//- new test : Ball2D, with a ball moving in a 2D system
//
//- another constant variables moved/refactored in XMLTagsName
//- making uniform the name of the constant variables
//
//Revision 1.12  2004/12/08 12:49:39  jbarbier
//- changes in the XML Schema, respect of the recommandations of the W3C
//version 1.1
//
//- changes in all balises DS, Relation, NSLaw, OneStepIntegrator, OneStepNSProblem
//in the XML files into specific names like LagrangianNLDS, LinearSystemDS, ...
//for the DS
//
//Revision 1.11  2004/09/23 14:09:24  jbarbier
//- modification of the integrators, the attribute r is always optional.
//
//- modification of the LagrangianNonLinearR. computeInput and computeOutput are
//required.
//
//Revision 1.10  2004/09/14 13:49:59  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.9  2004/08/03 12:07:12  jbarbier
//- all test on th eModel are successfull
//
//- new tests on the Model with the opening of XML file
//
//- link TimeDiscretisation -> Strategy
//
//- attribute T of the Model is now optional
//
//Revision 1.8  2004/07/29 14:25:45  jbarbier
