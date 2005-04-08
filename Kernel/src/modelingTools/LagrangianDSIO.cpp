
#include "LagrangianDSIO.h"
#include "check.h"

LagrangianDSIO::LagrangianDSIO()
{
  this->dsioType = LAGRANGIANDSIO;
}
LagrangianDSIO::~LagrangianDSIO()
{}


//void LagrangianDSIO::saveDSInputOutputToXML()
//{
//  DSInputOutput::saveDSInputOutputToXML();
//}

//void DSInputOutput::createDSInputOutput(DSInputOutputXML * dsioXML, int number,
//                    string computeInput, string computeOutput)
//{
//  if( dsioXML != NULL )
//  {
//////    this->init();
//    this->dsioxml = dsioXML;
//    this->dsioType = NLINEARDSIO;
//    this->fillDSInputOutputWithDSInputOutputXML();
//  }
//  else
//  {
//    this->dsioxml = NULL;
//    this->dsioType = NLINEARDSIO;
//    // computeInput
//    this->setComputeInputFunction(this->cShared.getPluginName( computeInput ), this->cShared.getPluginFunctionName( computeInput ));
//
//    // computeOutput
//    this->setComputeOutputFunction(this->cShared.getPluginName( computeOutput ), this->cShared.getPluginFunctionName( computeOutput ));
//    this->number = number;
//  }
//}

