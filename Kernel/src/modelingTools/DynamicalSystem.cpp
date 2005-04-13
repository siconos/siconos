#include "DynamicalSystem.h"

#include "LinearBC.h"
#include "NLinearBC.h"
#include "PeriodicBC.h"

#include "LinearDSIO.h"
#include "LagrangianDSIO.h"
#include "LagrangianLinearDSIO.h"


DynamicalSystem::DynamicalSystem()
{
  this->init();
  this->DSType = NLSDS;
  this->dsxml = NULL;
}

DynamicalSystem::DynamicalSystem(DSXML *dsxml)
{
  IN("DynamicalSystem::DynamicalSystem(DSXML *dsxml)\n");
  this->init();
  this->DSType = NLSDS;
  this->dsxml = dsxml;
  OUT("DynamicalSystem::DynamicalSystem(DSXML *dsxml)\n");
}

DynamicalSystem::~DynamicalSystem()
{

  IN("DynamicalSystem::~DynamicalSystem()\n");

  if (this->x != NULL) delete x;
  if (this->x0 != NULL) delete x0;
  if (this->xFree != NULL) delete xFree;
  for (int i = 0; i < this->dsioVector.size(); i++)
    if (this->dsioVector[i] != NULL) delete this->dsioVector[i];

  OUT("DynamicalSystem::~DynamicalSystem()\n");
}
////////////////////////////


vector<DSInputOutput*> DynamicalSystem::getDSInputOutputs(void)
{
  return dsioVector;
}

DSInputOutput* DynamicalSystem::getDSInputOutput(int i)
{
  if (i < this->dsioVector.size())
  {
    return this->dsioVector[i];
  }
  //cout<<"# i = "<<i<<endl;
  RuntimeException::selfThrow("EqualityConstraint - getDSInputOutput : \'i\' is out of range");
}

void DynamicalSystem::setDSInputOutputs(vector<DSInputOutput*> dsioVect)
{
  this->dsioVector = dsioVect;
}

void DynamicalSystem::addDSInputOutput(DSInputOutput* dsio)
{
  //  DSInputOutput* dsioTmp;
  //  dsioTmp = new DSInputOutput();
  //  *dsioTmp = *dsio;
  //  this->dsioVector.push_back( dsioTmp );
  this->dsioVector.push_back(dsio);
}

////////////////////////////

void DynamicalSystem::setVectorFieldFunction(string pluginPath, string functionName)
{
  this->vectorFieldPtr = NULL;
  cShared.setFunction(&vectorFieldPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->vectorFieldFunctionName = plugin + ":" + functionName;
}

void DynamicalSystem::setComputeJacobianXFunction(string pluginPath, string functionName)
{
  this->computeJacobianXPtr = NULL;
  cShared.setFunction(&computeJacobianXPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->computeJacobianXFunctionName = plugin + ":" + functionName;
}

void DynamicalSystem::vectorField(double time)
{
  if (vectorFieldPtr == NULL)
    RuntimeException::selfThrow("vectorField() is not linked to a plugin function");

  int size = x->size();
  this->vectorFieldPtr(&size, &time, &(*x)(0) , &xDot(0));
}

void DynamicalSystem::computeJacobianX(double time)
{
  if (computeJacobianXPtr == NULL)
    RuntimeException::selfThrow("computeJacobianX() is not linked to a plugin function");

  int size = x->size();
  this->computeJacobianXPtr(&size, &time, &(*x)(0), &jacobianX(0, 0));
}


void DynamicalSystem::swapInMemory(void)
{
  IN("DynamicalSystem::swapInMemory\n ");
  xMemory.swap(this->x);
  xDotMemory.swap(&this->xDot);
  rMemory.swap(&r);
  //cout<<"DynamicalSystem::swapInMemory Done..."<<endl;
  OUT("DynamicalSystem::swapInMemory\n ");
}


void DynamicalSystem::fillDSWithDSXML()
{
  IN("DynamicalSystem::fillDSWithDSXML\n");
  if (this->dsxml != NULL)
  {
    this->number = this->dsxml->getNumber();

    if (this->dsxml->hasId() == true) this->id = this->dsxml->getId();
    else cout << "Warning : Id is not defined in the XML " << endl;

    if (this->dsxml->hasN() == true) this->n = this->dsxml->getN();
    else cout << "Warning : n is not defined in the XML " << endl;

    if (this->dsxml->hasX0() == true)
      *(this->x0) = this->dsxml->getX0();
    else cout << "Warning : x0 is not defined in the XML " << endl;

    if (this->dsxml->hasX() == true)
    {
      *(this->x) = this->dsxml->getX();
    }
    else cout << "Warning : x is not defined in the XML " << endl;

    if (this->dsxml->hasXDot() == true)(this->xDot) = this->dsxml->getXDot();
    else cout << "Warning : xDot is not defined in the XML " << endl;

    if (this->dsxml->hasXMemory() == true) this->xMemory = SiconosMemory::SiconosMemory(this->dsxml->getXMemoryXML()); //this->dsxml->getXMemory();
    else cout << "Warning : xMemory is not defined in the XML " << endl;

    if (this->dsxml->hasXDotMemory() == true) this->xDotMemory = SiconosMemory::SiconosMemory(this->dsxml->getXDotMemoryXML()); //this->dsxml->getXDotMemory();
    else cout << "Warning : xDotMemory is not defined in the XML " << endl;

    if (this->dsxml->hasStepsInMemory() == true) this->stepsInMemory = this->dsxml->getStepsInMemory();
    else cout << "Warning : stepsInMemory is not defined in the XML " << endl;

    string plugin;
    // vectorField
    if (this->dsxml->hasVectorFieldPlugin() == true)
    {
      plugin = this->dsxml->getVectorFieldPlugin();
      this->setVectorFieldFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
    }
    else cout << "Warning : vectorFieldPlugin is not defined in the XML " << endl;


    // computeJacobianX
    if (this->dsxml->hasComputeJacobianXPlugin() == true)
    {
      plugin = this->dsxml->getComputeJacobianXPlugin();
      this->setComputeJacobianXFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
    }
    else cout << "Warning : computeJacobianXPlugin is not defined in the XML " << endl;

    this->r = SimpleVector(this->n);
  }
  else RuntimeException::selfThrow("DynamicalSystem::fillDSWithDSXML - DSXML object not exists");
  OUT("DynamicalSystem::fillDSWithDSXML\n");
}

void DynamicalSystem::display() const
{
  IN("DynamicalSystem::display\n");
  cout << "____ data of the Dynamical System " << endl;
  cout << "| number : " << this->number << endl;
  cout << "| id : " << this->id << endl;
  cout << "| n : " << this->n << endl;
  cout << "| x " << endl;
  this->x->display();
  cout << "| x0 " << endl;
  this->x0->display();
  cout << "| xFree " << endl;
  this->xFree->display();
  cout << "| xDot " << endl;
  this->xDot.display();
  cout << "| stepsInMemory : " << this->stepsInMemory << endl;
  cout << "| r " << endl;
  this->r.display();
  OUT("DynamicalSystem::display\n");

}

void DynamicalSystem::linkDSXML()
{
  IN("DynamicalSystem::linkDSXML\n");
  if (this->dsxml->getBoundaryConditionXML() != NULL)
  {
    //cout<<"#DynamicalSystem::linkDSXML - BC type == "<< this->dsxml->getBoundaryConditionXML()->getType() <<endl;
    if (this->dsxml->getBoundaryConditionXML()->getType() == LINEARBC_TAG)
    {
      // creation of the LinearBC with this constructor and call of a method to fill
      this->BC = new LinearBC();
      static_cast<LinearBC*>(this->BC)->createBoundaryCondition(this->dsxml->getBoundaryConditionXML());
    }

    else if (this->dsxml->getBoundaryConditionXML()->getType() == NON_LINEARBC_TAG)
    {
      // creation of the NLinearBC with this constructor and call of a method to fill
      this->BC = new NLinearBC();
      static_cast<NLinearBC*>(this->BC)->createBoundaryCondition(this->dsxml->getBoundaryConditionXML());
    }

    else if (this->dsxml->getBoundaryConditionXML()->getType() == PERIODICBC_TAG)
    {
      // creation of the PeriodicBC with this constructor and call of a method to fill
      this->BC = new PeriodicBC();
      static_cast<PeriodicBC*>(this->BC)->createBoundaryCondition(this->dsxml->getBoundaryConditionXML());
    }
    else RuntimeException::selfThrow("DynamicalSystem::linkDSXML - bad kind of BoundaryCondition : " + this->dsxml->getBoundaryConditionXML()->getType());
  }
  else this->BC = NULL;

  DSInputOutput *dsio;
  vector<int> nbDSIOtab = this->dsxml->getDSInputOutputNumbers();
  //cout<<"DS == "<<this->DSType<<" || size of DSIputOutput == "<<nbDSIOtab.size()<<endl;
  for (int i = 0; i < nbDSIOtab.size(); i++)
  {
    cout << "DynamicalSystem => linkDS, DSIputOutputNumbers == " << nbDSIOtab[i] << endl;
    if (this->dsxml->getDSInputOutputXML(nbDSIOtab[i])->getType() == LINEAR_DSIO_TAG)
    {
      dsio = new LinearDSIO();
      this->dsioVector.push_back(dsio);
      static_cast<LinearDSIO*>(dsio)->createDSInputOutput(this->dsxml->getDSInputOutputXML(nbDSIOtab[i]));
    }
    else if (this->dsxml->getDSInputOutputXML(nbDSIOtab[i])->getType() == NON_LINEAR_DSIO_TAG)
    {
      dsio = new DSInputOutput();
      this->dsioVector.push_back(dsio);
      static_cast<DSInputOutput*>(dsio)->createDSInputOutput(this->dsxml->getDSInputOutputXML(nbDSIOtab[i]));
    }
    else if (this->dsxml->getDSInputOutputXML(nbDSIOtab[i])->getType() == LAGRANGIAN_DSIO_TAG)
    {
      dsio = new LagrangianDSIO();
      this->dsioVector.push_back(dsio);
      static_cast<LagrangianDSIO*>(dsio)->createDSInputOutput(this->dsxml->getDSInputOutputXML(nbDSIOtab[i]));
    }
    else if (this->dsxml->getDSInputOutputXML(nbDSIOtab[i])->getType() == LAGRANGIAN_LINEAR_DSIO_TAG)
    {
      dsio = new LagrangianDSIO();
      this->dsioVector.push_back(dsio);
      static_cast<LagrangianLinearDSIO*>(dsio)->createDSInputOutput(this->dsxml->getDSInputOutputXML(nbDSIOtab[i]));
    }
    else RuntimeException::selfThrow("DynamicalSystem::linkDSXML - bad kind of BoundaryCondition : " + this->dsxml->getBoundaryConditionXML()->getType());
  }
  OUT("DynamicalSystem::linkDSXML\n");
}


void DynamicalSystem::init()
{
  IN("DynamicalSystem::init\n");
  //this->nsds = NULL;
  this->number = 0;
  this->id = "none";
  this->n = 0;

  this->x0 = new SimpleVector();
  this->x = new SimpleVector();
  this->xDot = SimpleVector::SimpleVector();
  this->xFree = new SimpleVector();

  this->r = SimpleVector::SimpleVector();
  this->BC = NULL;

  this->jacobianX = SiconosMatrix::SiconosMatrix();

  this->stepsInMemory = 1;
  this->setVectorFieldFunction("BasicPlugin.so", "vectorField");
  this->setComputeJacobianXFunction("BasicPlugin.so", "computeJacobianX");
  this->dsxml = NULL;

  OUT("DynamicalSystem::init\n");
}


void DynamicalSystem::initMemory(int steps)
{
  IN("DynamicalSystem::initMemory\n");
  if (steps < 0)
    RuntimeException::selfThrow("DynamicalSystem::initMemory(int steps) - steps < 0");
  else
  {
    this->stepsInMemory = steps;

    /*
     ** we made the initialization of the memories
     *
     * for rMemory, we don't need to load data for the DOM tree because there are no data saved in the XML for r
     *
     * the other memories are resized with the first parameter 'steps', and data are reloaded from the DOM tree
     * only if there are data in the DOM tree
     */

    this->rMemory = SiconosMemory::SiconosMemory(steps);
    this->xMemory = SiconosMemory::SiconosMemory(steps, this->xMemory.getSiconosMemoryXML());
    this->xDotMemory = SiconosMemory::SiconosMemory(steps, this->xDotMemory.getSiconosMemoryXML());
  }

  OUT("DynamicalSystem::initMemory\n");
}


void DynamicalSystem::saveDSToXML()
{
  IN("DynamicalSystem::saveDSToXML\n");

  /*
   * save of the BoundariesConditions
   */
  if (this->BC != NULL)
  {
    if (this->BC->getType() == LINEARBC)
      (static_cast<LinearBC*>(this->BC))->saveBCToXML();
    else if (this->BC->getType() == NLINEARBC)
      (static_cast<NLinearBC*>(this->BC))->saveBCToXML();
    else if (this->BC->getType() == PERIODICBC)
      (static_cast<PeriodicBC*>(this->BC))->saveBCToXML();
    else RuntimeException::selfThrow("DynamicalSystem::saveDSToXML - bad kind of BoundaryCondition");
  }

  if (this->dsioVector.size() != 0)
  {
    for (int i = 0; i < this->dsioVector.size(); i++)
    {
      if (this->dsioVector[i]->getType() == LINEARDSIO)
        (static_cast<LinearDSIO*>(this->dsioVector[i]))->saveDSInputOutputToXML();
      else if (this->dsioVector[i]->getType() == NLINEARDSIO)
        (static_cast<DSInputOutput*>(this->dsioVector[i]))->saveDSInputOutputToXML();
      else if (this->dsioVector[i]->getType() == LAGRANGIANDSIO)
        (static_cast<LagrangianDSIO*>(this->dsioVector[i]))->saveDSInputOutputToXML();
      else RuntimeException::selfThrow("DynamicalSystem::saveDSToXML - bad kind of DSInputOuput");
    }
  }

  if (this->dsxml != NULL)
  {
    this->dsxml->setId(this->id);
    if ((this->dsxml->getType() != LNLDS)
        && (this->dsxml->getType() != LTIDS))
      this->dsxml->setN(this->n);

    this->dsxml->setX0(this->x0);

    this->dsxml->setX(this->x);
    this->dsxml->setXMemory(&(this->xMemory));

    this->dsxml->setXDot(&this->xDot);
    this->dsxml->setXDotMemory(&(this->xDotMemory));

    this->dsxml->setStepsInMemory(this->stepsInMemory);

    this->dsxml->setR(&(this->r));

    /*
     * vectorField and computeJacobianX function must be saved only for NonLinearSystemDS
     */
    if (this->DSType == NLSDS)
    {
      this->dsxml->setVectorFieldPlugin(this->vectorFieldFunctionName);
      this->dsxml->setComputeJacobianXPlugin(this->computeJacobianXFunctionName);
    }
  }
  else RuntimeException::selfThrow("DynamicalSystem::saveDSToXML - The DSXML object doesn't exists");
  OUT("DynamicalSystem::saveDSToXML\n");
}

void DynamicalSystem::createDynamicalSystem(DSXML * dsXML, int number, int n,
    SiconosVector* x0, string vectorFieldPlugin)//, NonSmoothDynamicalSystem * nsds, BoundaryCondition* bc)
{
  IN("DynamicalSystem::createDynamicalSystem\n");
  if (dsXML != NULL)
  {
    this->DSType = NLSDS;
    //this->init();
    this->dsxml = dsXML;

    this->fillDSWithDSXML();
    this->linkDSXML();
  }
  else
  {
    this->DSType = NLSDS;
    this->number = number;
    this->n = n;

    this->x0 = new SimpleVector(n);
    this->x = new SimpleVector(n);
    this->xDot = /*new*/ SimpleVector::SimpleVector(n);
    this->xFree = new SimpleVector(n);

    *(this->x0) = *x0;
    this->setVectorFieldFunction(this->cShared.getPluginName(vectorFieldPlugin), this->cShared.getPluginFunctionName(vectorFieldPlugin));

  }
  OUT("DynamicalSystem::createDynamicalSystem\n");
}

BoundaryCondition* DynamicalSystem::createPeriodicBC()
{
  this->BC = new PeriodicBC();
  static_cast<PeriodicBC*>(this->BC)->createBoundaryCondition(NULL);
  return this->BC;
}

BoundaryCondition* DynamicalSystem::createLinearBC(SiconosVector* omega, SiconosMatrix* omega0, SiconosMatrix* omegaT)
{
  this->BC = new LinearBC();
  static_cast<LinearBC*>(this->BC)->createBoundaryCondition(NULL, omega, omega0, omegaT);
  return this->BC;
}

BoundaryCondition* DynamicalSystem::createNLinearBC()
{
  this->BC = new NLinearBC();
  static_cast<NLinearBC*>(this->BC)->createBoundaryCondition(NULL);
  return this->BC;
}

