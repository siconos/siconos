
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

