#include "Relation.h"
using namespace std;

Relation::Relation(): relationType("none"), interaction(NULL), relationxml(NULL),
  computeInputName("none"), computeOutputName("none")
{
  setComputeOutputFunction("BasicPlugin.so", "computeOutput");
  setComputeInputFunction("BasicPlugin.so", "computeInput");
}

Relation::Relation(RelationXML* relxml): relationType("none"), interaction(NULL),
  relationxml(relxml), computeInputName("none"),
  computeOutputName("none")

{
  setComputeOutputFunction("BasicPlugin.so", "computeOutput");
  setComputeInputFunction("BasicPlugin.so", "computeInput");
  if (relationxml != NULL)
  {
    string plugin;

    // computeInput
    if (relationxml->hasComputeInput())
    {
      plugin = (relationxml)->getComputeInputPlugin();
      setComputeInputFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    // computeOutput
    if (relationxml->hasComputeOutput())
    {
      plugin = (relationxml)->getComputeOutputPlugin();
      setComputeOutputFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
  }
  else RuntimeException::selfThrow("Relation::fillRelationWithRelationXML - object RelationXML does not exist");
}

Relation::~Relation()
{}


vector<DSInputOutput*> Relation::getDSInputOutputs(void)
{
  return dsioVector;
}

DSInputOutput* Relation::getDSInputOutput(const unsigned int& i)
{
  if (i >= dsioVector.size())
    RuntimeException::selfThrow("Relation - getDSInputOutput : \'i\' is out of range");
  return dsioVector[i];
}

void Relation::setDSInputOutputs(vector<DSInputOutput*> dsioVect)
{
  dsioVector = dsioVect;
}

void Relation::addDSInputOutput(DSInputOutput* dsio)
{
  /*
   *  in EqualityConstraint class, we don't create new objects in the DSInputOutput vector
   *    => we only save a link (pointer) on the DSInputOutputs of the DynamicalSystems !!
   */
  dsioVector.push_back(dsio);
}



void Relation::computeOutput(const double& time)
{
  //to do
  RuntimeException::selfThrow("Relation - computeOutput: not yet implemented for relation of type" + getType());
}

void Relation::computeFreeOutput(const double& time)
{
  RuntimeException::selfThrow("Relation - computeFreeOutput: not yet implemented for relation of type" + getType());
  //to do
}

void Relation::computeInput(const double& time)
{
  //to do
  RuntimeException::selfThrow("Relation - computeIntput: not yet implemented for relation of type" + getType());

}

void Relation::setComputeOutputFunction(const string& pluginPath, const string& functionName)
{
  computeOutputPtr = NULL;
  cShared.setFunction(&computeOutputPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  computeOutputName = plugin + ":" + functionName;
}

void Relation::setComputeInputFunction(const string& pluginPath, const string& functionName)
{
  computeInputPtr = NULL;
  cShared.setFunction(&computeInputPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  computeInputName = plugin + ":" + functionName;
}

