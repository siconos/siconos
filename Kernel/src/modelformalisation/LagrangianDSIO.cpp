//$Id: LagrangianDSIO.cpp,v 1.5 2005/03/14 16:05:27 jbarbier Exp $

#include "LagrangianDSIO.h"
#include "check.h"

LagrangianDSIO::LagrangianDSIO()
{
  this->dsioType = LAGRANGIANDSIO;
}
LagrangianDSIO::~LagrangianDSIO()
{}


//void LagrangianDSIO::saveDSInputOutputToXML()
//{}
//
void LagrangianDSIO::createDSInputOutput(DSInputOutputXML * dsioXML, int number,
    SiconosMatrix *H)
{
  if (dsioXML != NULL)
  {
    //    this->init();
    this->dsioxml = dsioXML;
    this->dsioType = LAGRANGIANDSIO;
    this->fillDSInputOutputWithDSInputOutputXML();
  }
  else
  {
    //this->dsioxml = dsioXML;
    this->dsioType = LAGRANGIANDSIO;
    this->number = number;
    this->H = *H;
    // computeInput
    //    this->setComputeInputFunction(this->cShared.getPluginName( computeInput ), this->cShared.getPluginFunctionName( computeInput ));
    //
    //    // computeOutput
    //    this->setComputeOutputFunction(this->cShared.getPluginName( computeOutput ), this->cShared.getPluginFunctionName( computeOutput ));
  }
}

//$Log: LagrangianDSIO.cpp,v $
//Revision 1.5  2005/03/14 16:05:27  jbarbier
//- manual creation of DSInputOutput saving OK
//
//- in progress for EqualityConstraint
//
//Revision 1.4  2005/03/10 12:55:20  jbarbier
//- implmentation of the EqualityConstraint and DSInputOutput classes in progress
//    attributes H (DSIO) et G (EC) added in XML and managed in XML objects
//
//Revision 1.3  2005/03/09 15:30:28  jbarbier
//- add of LagrangianEC class
//
//- in progress : implementation of the EqualityConstraint and DSInputOutput - create methods
//
//Revision 1.2  2005/03/08 12:41:36  jbarbier
//- constant variables files modified :
//Some constants added in SiconosConst
//
//all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//
//Revision 1.1  2005/01/17 10:56:24  jbarbier
//- classes EqualityConstraint and DSInputOutput added with inherited classes
//
//- classes EqualityConstraintXML and DSInputOutputXML added with inherited classes
//
