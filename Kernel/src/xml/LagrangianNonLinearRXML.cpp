//$Id: LagrangianNonLinearRXML.cpp,v 1.2 2004/09/14 13:49:58 jbarbier Exp $

#include "LagrangianNonLinearRXML.h"


LagrangianNonLinearRXML::LagrangianNonLinearRXML(): RelationXML()
{}

LagrangianNonLinearRXML::LagrangianNonLinearRXML(xmlNode * LNLRelationNode)
  : RelationXML(LNLRelationNode)
{}

LagrangianNonLinearRXML::~LagrangianNonLinearRXML()
{}

//$Log: LagrangianNonLinearRXML.cpp,v $
//Revision 1.2  2004/09/14 13:49:58  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.1  2004/08/12 11:55:19  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
