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

DSInputOutput* Relation::getDSInputOutput(const int& i)
{
  if (i < dsioVector.size())
    return dsioVector[i];
  else RuntimeException::selfThrow("Relation - getDSInputOutput : \'i\' is out of range");
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
  if (computeOutputPtr == NULL) RuntimeException::selfThrow("computeOutput() is not linked to a plugin function");

  //to do
  //computeOutputPtr(&x(0), &time, &lambdaPtr(0), &y(0));
  //  vector<DynamicalSystem*> vDS = interaction->getDynamicalSystems();
  //
  //  DynamicalSystem *ds1 ,*ds2;
  //  SiconosVector *y = interaction->getYPtr();
  //  SiconosVector *yDot = interaction->getYDotPtr();
  //  if (vDS.size() == 2)
  //  {
  //      ds1=vDS[0];
  //      ds2=vDS[1];
  //      if (((ds1->getType() == LNLDS) || (ds1->getType() == LTIDS)) && ((ds2->getType() == LNLDS) || (ds2->getType() == LTIDS)))
  //      {
  //        LagrangianDS *d1 = static_cast<LagrangianDS*> (ds1);
  //        LagrangianDS *d2 = static_cast<LagrangianDS*> (ds2);
  //
  //        CompositeVector q;
  //        q.add(*(d1->getQPtr()));
  //      q.add(*(d2->getQPtr()));
  //        //*y = (h * q) + b;
  //
  //      CompositeVector vel;
  //      vel.add(*(d1->getVelocityPtr()));
  //      vel.add(*(d2->getVelocityPtr()));
  //      *yDot = (h * vel);
  //
  //      computeOutputPtr(*q, 0.0, lambda, y);
  //      }
  //    else
  //    {
  //      // To be Finished
  //      RuntimeException::selfThrow("LagrangianLinearR::computeOutput not yet implemented for this type of dynamical system "+vDS[0]->getType());
  //    }
}

void Relation::computeFreeOutput(const double& time)
{
  if (computeOutputPtr == NULL) RuntimeException::selfThrow("computeFreeOutput() is not linked to a plugin function");

  //to do
  //computeOutputPtr(&xFree(0), &time, &lambdaPtr(0), &y(0));
}

void Relation::computeInput(const double& time)
{
  if (computeInputPtr == NULL) RuntimeException::selfThrow("computeInput() is not linked to a plugin function");

  //to do
  //computeInputPtr(&x(0), &time, &lambdaPtr(0), &r(0));
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

