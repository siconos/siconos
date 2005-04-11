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
  cout << "# # # # # # # # # DSInputOutput(DSInputOutputXML* dsioxml)" << endl;
  this->init();
  this->dsioType = NLINEARDSIO;
  this->dsioxml = dsioxml;
}

DSInputOutput::~DSInputOutput()
{}

void DSInputOutput::setComputeOutputFunction(std::string pluginPath, std::string functionName)
{
  this->computeOutputPtr = NULL;
  cShared.setFunction(&computeOutputPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->computeOutputName = plugin + ":" + functionName;
}

void DSInputOutput::setComputeInputFunction(std::string pluginPath, std::string functionName)
{
  this->computeInputPtr = NULL;
  cShared.setFunction(&computeInputPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->computeInputName = plugin + ":" + functionName;
}

void DSInputOutput::display() const
{
  cout << "-----------------------------------------------------" << endl;
  cout << "____ data of the DSInputOutput " << endl;
  cout << "| id : " << this->id << endl;
  cout << "| number : " << this->number << endl;
  cout << "| computeInput plugin name : " << this->computeInputName << endl;
  cout << "| computeOutput plugin name : " << this->computeOutputName << endl;
  cout << "-----------------------------------------------------" << endl << endl;
}

void DSInputOutput::fillDSInputOutputWithDSInputOutputXML()
{
  IN("DSInputOutput::fillDSInputOutputWithDSInputOutputXML\n");
  if (this->dsioxml != NULL)
  {
    string plugin;

    // computeInput
    if (this->dsioxml->hasComputeInput())
    {
      cout << "DSInputOutputPluginType == " << this->dsioType << endl;
      plugin = (this->dsioxml)->getComputeInputPlugin();
      this->setComputeInputFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
    }
    else cout << "Warning - No computeInput method is defined in a DSInputOutput " << this->getType() << endl;

    // computeOutput
    if (this->dsioxml->hasComputeOutput())
    {
      cout << "DSInputOutputPluginType == " << this->dsioType << endl;
      plugin = (this->dsioxml)->getComputeOutputPlugin();
      this->setComputeOutputFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
    }
    else cout << "Warning - No computeOutput method is defined in a DSInputOutput " << this->getType() << endl;

    this->number = this->dsioxml->getNumber();
    //this->H = this->dsioxml->getH();
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

  this->setComputeOutputFunction("BasicPlugin.so", "computeOutput");
  this->setComputeInputFunction("BasicPlugin.so", "computeInput");
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
    this->dsioxml->setComputeInputPlugin(this->computeInputName);
    this->dsioxml->setComputeOutputPlugin(this->computeOutputName);
  }
  else RuntimeException::selfThrow("DSInputOutput::saveDSInputOutputToXML - object DSInputOutputXML does not exist");
  OUT("DSInputOutput::saveDSInputOutputToXML\n");
}

void DSInputOutput::createDSInputOutput(DSInputOutputXML * dsioXML, int number,
                                        string computeInput, string computeOutput)
{
  if (dsioXML != NULL)
  {
    ////    this->init();
    this->dsioxml = dsioXML;
    this->dsioType = NLINEARDSIO;
    this->fillDSInputOutputWithDSInputOutputXML();
  }
  else
  {
    this->dsioxml = NULL;
    this->dsioType = NLINEARDSIO;
    // computeInput
    this->setComputeInputFunction(this->cShared.getPluginName(computeInput), this->cShared.getPluginFunctionName(computeInput));

    // computeOutput
    this->setComputeOutputFunction(this->cShared.getPluginName(computeOutput), this->cShared.getPluginFunctionName(computeOutput));
    this->number = number;
  }
}
