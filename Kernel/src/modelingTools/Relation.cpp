#include "Relation.h"

#include "check.h"


Relation::Relation()
{
  this->init();
  this->relationxml = NULL;
  this->interaction = NULL;
}

Relation::Relation(RelationXML* relxml)
{
  this->init();
  this->relationxml = relxml;
  this->interaction = NULL;
}

Relation::~Relation()
{}


vector<DSInputOutput*> Relation::getDSInputOutputs(void)
{
  return dsioVector;
}

DSInputOutput* Relation::getDSInputOutput(int i)
{
  if (i < this->dsioVector.size())
  {
    return this->dsioVector[i];
  }
  RuntimeException::selfThrow("Relation - getDSInputOutput : \'i\' is out of range");
}

void Relation::setDSInputOutputs(vector<DSInputOutput*> dsioVect)
{
  this->dsioVector = dsioVect;
}

void Relation::addDSInputOutput(DSInputOutput* dsio)
{
  /*
   *  in EqualityConstraint class, we don't create new objects in the DSInputOutput vector
   *    => we only save a link (pointer) on the DSInputOutputs of the DynamicalSystems !!
   */
  this->dsioVector.push_back(dsio);
}



void Relation::computeOutput(double time)
{
  if (computeOutputPtr == NULL) RuntimeException::selfThrow("computeOutput() is not linked to a plugin function");

  //to do
  //this->computeOutputPtr(&x(0), &time, &lambdaPtr(0), &y(0));
  //  vector<DynamicalSystem*> vDS = this->interaction->getDynamicalSystems();
  //
  //  DynamicalSystem *ds1 ,*ds2;
  //  SiconosVector *y = this->interaction->getYPtr();
  //  SiconosVector *yDot = this->interaction->getYDotPtr();
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
  //        //*y = (this->h * q) + this->b;
  //
  //      CompositeVector vel;
  //      vel.add(*(d1->getVelocityPtr()));
  //      vel.add(*(d2->getVelocityPtr()));
  //      *yDot = (this->h * vel);
  //
  //      this->computeOutputPtr(*q, 0.0, this->lambda, y);
  //      }
  //    else
  //    {
  //      // To be Finished
  //      RuntimeException::selfThrow("LagrangianLinearR::computeOutput not yet implemented for this type of dynamical system "+vDS[0]->getType());
  //    }
}
void Relation::computeFreeOutput(double time)
{
  if (computeOutputPtr == NULL) RuntimeException::selfThrow("computeFreeOutput() is not linked to a plugin function");

  //to do
  //this->computeOutputPtr(&xFree(0), &time, &lambdaPtr(0), &y(0));
}

void Relation::computeInput(double time)
{
  if (computeInputPtr == NULL) RuntimeException::selfThrow("computeInput() is not linked to a plugin function");

  //to do
  //this->computeInputPtr(&x(0), &time, &lambdaPtr(0), &r(0));
}

void Relation::setComputeOutputFunction(std::string pluginPath, std::string functionName)
{
  this->computeOutputPtr = NULL;
  cShared.setFunction(&computeOutputPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->computeOutputName = plugin + ":" + functionName;
}

void Relation::setComputeInputFunction(std::string pluginPath, std::string functionName)
{
  this->computeInputPtr = NULL;
  cShared.setFunction(&computeInputPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->computeInputName = plugin + ":" + functionName;
}

void Relation::fillRelationWithRelationXML()
{
  IN("Relation::fillRelationWithRelationXML\n");
  if (this->relationxml != NULL)
  {
    string plugin;

    // computeInput
    if (this->relationxml->hasComputeInput())
    {
      cout << "RelationPluginType == " << this->relationType << endl;
      plugin = (this->relationxml)->getComputeInputPlugin();
      this->setComputeInputFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
    }
    else cout << "Warning - No computeInput method is defined in a Relation " << this->getType() << endl;

    // computeOutput
    if (this->relationxml->hasComputeOutput())
    {
      cout << "RelationPluginType == " << this->relationType << endl;
      plugin = (this->relationxml)->getComputeOutputPlugin();
      this->setComputeOutputFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
    }
    else cout << "Warning - No computeOutput method is defined in a Relation " << this->getType() << endl;
  }
  else RuntimeException::selfThrow("Relation::fillRelationWithRelationXML - object RelationXML does not exist");

  OUT("Relation::fillRelationWithRelationXML\n");
}

void Relation::init()
{
  IN("Relation::init\n");
  this->setComputeOutputFunction("BasicPlugin.so", "computeOutput");
  this->setComputeInputFunction("BasicPlugin.so", "computeInput");
  OUT("Relation::init\n");
}

void Relation::saveRelationToXML()
{
  IN("Relation::saveRelationToXML\n");
  if (this->relationxml != NULL)
  {
    /*
     * these attributes are only required for LagrangianNonLinear relation !
     */
  }
  else RuntimeException::selfThrow("Relation::saveRelationToXML - object RelationXML does not exist");
  OUT("Relation::saveRelationToXML\n");
}
