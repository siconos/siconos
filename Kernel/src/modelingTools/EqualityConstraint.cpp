#include "EqualityConstraint.h"

EqualityConstraint::EqualityConstraint()
{
  this->ecXML = NULL;
  this->type = NLINEAREC;
}

EqualityConstraint::EqualityConstraint(EqualityConstraintXML* ecXML)
{
  this->ecXML = ecXML;
  this->type = NLINEAREC;
}

EqualityConstraint::~EqualityConstraint()
{}



vector<DSInputOutput*> EqualityConstraint::getDSInputOutputs(void)
{
  return dsioVector;
}

DSInputOutput* EqualityConstraint::getDSInputOutput(int i)
{
  if (i < this->dsioVector.size())
  {
    return this->dsioVector[i];
  }
  RuntimeException::selfThrow("EqualityConstraint - getDSInputOutput : \'i\' is out of range");
}

void EqualityConstraint::setDSInputOutputs(vector<DSInputOutput*> dsioVect)
{
  this->dsioVector = dsioVect;
}

void EqualityConstraint::addDSInputOutput(DSInputOutput* dsio)
{
  //  DSInputOutput* dsioTmp;
  //  dsioTmp = new DSInputOutput();
  //  *dsioTmp = *dsio;

  /*
   *  in EqualityConstraint class, we don't create new objects in the DSInputOutput vector
   *    => we only save a link (pointer) on the DSInputOutputs of the DynamicalSystems !!
   */
  this->dsioVector.push_back(dsio);
}

void EqualityConstraint::saveEqualityConstraintToXML()
{
  IN("EqualityConstraint::saveEqualityConstraintToXML\n");
  if (this->ecXML != NULL)
  {
    /*
     * these attributes are only required for LagrangianNonLinear DSInputOutput !
     */
    //    this->disoxml->setComputeInputPlugin( this->computeInputName );
    //    this->dsioxml->setComputeOutputPlugin( this->computeOutputName );
    this->ecXML->setG(&this->G);
  }
  else RuntimeException::selfThrow("EqualityConstraint::saveEqualityConstraintToXML - object EqualityConstraintXML does not exist");
  OUT("EqualityConstraint::saveEqualityConstraintToXML\n");
}

void EqualityConstraint::display() const
{
  cout << "-----------------------------------------------------" << endl;
  cout << "____ data of the EqualityConstraint " << endl;
  cout << "| id : " << this->id << endl;
  cout << "| number : " << this->number << endl;
  cout << "| G : " << endl;
  this->G.display();
  cout << "-----------------------------------------------------" << endl << endl;
}

void EqualityConstraint::init()
{
  this->number = 0;
  this->id = "none";
  this->ecXML = NULL;
}

////////////////////////////////
void EqualityConstraint::computeOutput(double time)
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

void EqualityConstraint::computeInput(double time)
{
  if (computeInputPtr == NULL) RuntimeException::selfThrow("computeInput() is not linked to a plugin function");

  //to do
  //this->computeInputPtr(&x(0), &time, &lambdaPtr(0), &r(0));
}

void EqualityConstraint::setComputeOutputFunction(std::string pluginPath, std::string functionName)
{
  this->computeOutputPtr = NULL;
  cShared.setFunction(&computeOutputPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->computeOutputName = plugin + ":" + functionName;
}

void EqualityConstraint::setComputeInputFunction(std::string pluginPath, std::string functionName)
{
  this->computeInputPtr = NULL;
  cShared.setFunction(&computeInputPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->computeInputName = plugin + ":" + functionName;
}
///////////////////////////////

void EqualityConstraint::fillEqualityConstraintWithEqualityConstraintXML()
{
  if (this->ecXML != NULL)
  {
    string plugin;
    // computeInput
    if (this->ecXML->hasComputeInput())
    {
      cout << "EqualityConstraintPluginType == " << this->type << endl;
      plugin = (this->ecXML)->getComputeInputPlugin();
      this->setComputeInputFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
    }
    else cout << "Warning - No computeInput method is defined in a EqualityConstraint " << this->getType() << endl;

    // computeOutput
    if (this->ecXML->hasComputeOutput())
    {
      cout << "EqualityConstraintPluginType == " << this->type << endl;
      plugin = (this->ecXML)->getComputeOutputPlugin();
      this->setComputeOutputFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
    }
    else cout << "Warning - No computeOutput method is defined in a Relation " << this->getType() << endl;

    this->number = this->ecXML->getNumber();
    this->G = this->ecXML->getG();
  }
  //else RuntimeException::selfThrow("EqualityConstraint::fillEqualityConstraintWithEqualityConstraintXML - object EqualityConstraintXML does not exist");
}

void EqualityConstraint::createEqualityConstraint(EqualityConstraintXML *ecXML,
    int number,  SiconosMatrix *G,
    vector<DSInputOutput*> *dsioVector)
{
  if (ecXML != NULL)
  {
    this->ecXML = ecXML;
    this->type = NLINEAREC;
    this->fillEqualityConstraintWithEqualityConstraintXML();
  }
  else
  {
    this->ecXML = NULL;
    this->type = NLINEAREC;
    this->number = number;
    this->G = *G;
    this->dsioVector = *dsioVector;
  }
}

