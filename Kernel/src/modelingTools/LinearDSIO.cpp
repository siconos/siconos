
#include "LinearDSIO.h"
#include "check.h"

LinearDSIO::LinearDSIO(): DSInputOutput()
{
  this->dsioType = LINEARDSIO;
}

LinearDSIO::LinearDSIO(DSInputOutputXML* dsioxml): DSInputOutput(dsioxml)
{
  this->dsioType = LINEARDSIO;
}

LinearDSIO::~LinearDSIO()
{}


//void LinearDSIO::fillDSInputOutputWithDSInputOutputXML()
//{}
//
//void LinearDSIO::saveDSInputOutputToXML()
//{}

void LinearDSIO::createDSInputOutput(DSInputOutputXML * dsioXML, int number,
                                     SiconosMatrix *H)
{
  if (dsioXML != NULL)
  {
    //    this->init();
    this->dsioxml = dsioXML;
    this->dsioType = LINEARDSIO;
    this->fillDSInputOutputWithDSInputOutputXML();
  }
  else
  {
    //this->dsioxml = dsioXML;
    this->dsioType = LINEARDSIO;
    this->number = number;
    this->H = *H;
    // computeInput
    //    this->setComputeInputFunction(this->cShared.getPluginName( computeInput ), this->cShared.getPluginFunctionName( computeInput ));
    //
    //    // computeOutput
    //    this->setComputeOutputFunction(this->cShared.getPluginName( computeOutput ), this->cShared.getPluginFunctionName( computeOutput ));
  }
}

