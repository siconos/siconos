//$Id: DSInputOutput.cpp,v 1.7 2005/03/14 16:05:26 jbarbier Exp $
#include "DSInputOutput.h"

#include "check.h"


DSInputOutput::DSInputOutput()
{
  this->init();
  this->dsioxml = NULL;
  this->dsioType = NLINEARDSIO;
}

DSInputOutput::DSInputOutput(DSInputOutputXML* dsioxml)
{
  this->init();
  this->dsioType = NLINEARDSIO;
  this->dsioxml = dsioxml;
}

DSInputOutput::~DSInputOutput()
{}

void DSInputOutput::display() const
{
  cout << "-----------------------------------------------------" << endl;
  cout << "____ data of the DSInputOutput " << endl;
  cout << "| id : " << this->id << endl;
  cout << "| number : " << this->number << endl;
  cout << "| H : " << endl;
  this->H.display();
  cout << "-----------------------------------------------------" << endl << endl;
}

void DSInputOutput::fillDSInputOutputWithDSInputOutputXML()
{
  IN("DSInputOutput::fillDSInputOutputWithDSInputOutputXML\n");
  if (this->dsioxml != NULL)
  {
    //    string plugin;
    //
    //    // computeInput
    //    if( this->dsioxml->hasComputeInput() )
    //    {
    //      cout<<"DSInputOutputPluginType == "<< this->dsioType <<endl;
    //      plugin = (this->dsioxml)->getComputeInputPlugin();
    //      this->setComputeInputFunction(this->cShared.getPluginName( plugin ), this->cShared.getPluginFunctionName( plugin ));
    //    }
    //    else cout<<"Warning - No computeInput method is defined in a DSInputOutput "<< this->getType() <<endl;
    //
    //    // computeOutput
    //    if( this->dsioxml->hasComputeOutput() )
    //    {
    //      cout<<"DSInputOutputPluginType == "<< this->dsioType <<endl;
    //      plugin = (this->dsioxml)->getComputeOutputPlugin();
    //      this->setComputeOutputFunction(this->cShared.getPluginName( plugin ), this->cShared.getPluginFunctionName( plugin ));
    //    }
    //    else cout<<"Warning - No computeOutput method is defined in a DSInputOutput "<< this->getType() <<endl;

    this->number = this->dsioxml->getNumber();
    this->H = this->dsioxml->getH();
  }
  else RuntimeException::selfThrow("DSInputOutput::fillDSInputOutputWithDSInputOutputXML - object DSInputOutputXML does not exist");

  OUT("DSInputOutput::fillDSInputOutputWithDSInputOutputXML\n");
}

void DSInputOutput::init()
{
  IN("DSInputOutput::init\n");
  this->number = 0;
  this->id = "none";
  this->dsioxml = NULL;

  //  this->setComputeOutputFunction("BasicPlugin.so", "computeOutput");
  //  this->setComputeInputFunction("BasicPlugin.so", "computeInput");
  OUT("DSInputOutput::init\n");
}

void DSInputOutput::saveDSInputOutputToXML()
{
  IN("DSInputOutput::saveDSInputOutputToXML\n");
  if (this->dsioxml != NULL)
  {
    /*
     * these attributes are only required for LagrangianNonLinear DSInputOutput !
     */
    //    this->disoxml->setComputeInputPlugin( this->computeInputName );
    //    this->dsioxml->setComputeOutputPlugin( this->computeOutputName );
    this->dsioxml->setH(&(this->H));
  }
  else RuntimeException::selfThrow("DSInputOutput::saveDSInputOutputToXML - object DSInputOutputXML does not exist");
  OUT("DSInputOutput::saveDSInputOutputToXML\n");
}

void DSInputOutput::createDSInputOutput(DSInputOutputXML * dsioXML, int number, SiconosMatrix *H)
{
  if (dsioXML != NULL)
  {
    ////    this->init();
    this->dsioxml = dsioXML;
    this->dsioType = NLINEARDSIO;
    this->fillDSInputOutputWithDSInputOutputXML();
  }
  //  else
  {
    this->dsioxml = NULL;
    this->dsioType = NLINEARDSIO;
    //    // computeInput
    ////    this->setComputeInputFunction(this->cShared.getPluginName( computeInput ), this->cShared.getPluginFunctionName( computeInput ));
    ////
    ////    // computeOutput
    ////    this->setComputeOutputFunction(this->cShared.getPluginName( computeOutput ), this->cShared.getPluginFunctionName( computeOutput ));
    this->number = number;
    this->H = *H;
  }
}
