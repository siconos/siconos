//$Id: EqualityConstraint.cpp,v 1.10 2005/03/21 16:48:02 jbarbier Exp $
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
  //        LagrangianNLDS *d1 = static_cast<LagrangianNLDS*> (ds1);
  //        LagrangianNLDS *d2 = static_cast<LagrangianNLDS*> (ds2);
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

//$Log: EqualityConstraint.cpp,v $
//Revision 1.10  2005/03/21 16:48:02  jbarbier
//- EqualityConstraint : computeInput and computeOutput functions added (plugin funcitons)
//
//- link OneStepNSProblem - EqualityConstraint established
//
//- modification of OneStepNSProblem save according to change to normType[64] in SiconosNumerics.h
//
//Revision 1.9  2005/03/15 14:44:03  jbarbier
//- pySiconos.i edited to remove local paths
//
//- checkCoherency checks whether the DSInputOutputs and EqualityConstraints have unique numbers
//
//Revision 1.8  2005/03/15 09:57:47  jbarbier
//- EqualityConstraint save OK
//
//Revision 1.7  2005/03/14 16:05:27  jbarbier
//- manual creation of DSInputOutput saving OK
//
//- in progress for EqualityConstraint
//
//Revision 1.6  2005/03/11 15:06:20  jbarbier
//- save to XML methods of EqualityConstraint and DSInputOutput added
//
//- XML loading process modified : Model loads NSDS, then NSDS loads the DynamicalSystems, EqualityConstraints, Interactions; Modle loads Strategy, then Strategy loads TimeDiscretisation, then the Integrators, then the OneStepNSProblem
//
//Revision 1.5  2005/03/10 12:55:19  jbarbier
//- implmentation of the EqualityConstraint and DSInputOutput classes in progress
//    attributes H (DSIO) et G (EC) added in XML and managed in XML objects
//
//Revision 1.4  2005/03/09 15:30:25  jbarbier
//- add of LagrangianEC class
//
//- in progress : implementation of the EqualityConstraint and DSInputOutput - create methods
//
//Revision 1.3  2005/01/25 14:51:46  jbarbier
//- attributes id, type and XML object added to EqualityConstraint
//
//Revision 1.2  2005/01/25 13:56:03  jbarbier
//- link DynamicalSystem-DSInputOutput, NonSmoothDynamicalSystem-EqualityConstraint, EquaityConstraint-DSInputOutput and Relation-DSInputOutput available
//
//Revision 1.1  2005/01/17 10:56:24  jbarbier
//- classes EqualityConstraint and DSInputOutput added with inherited classes
//
//- classes EqualityConstraintXML and DSInputOutputXML added with inherited classes
//