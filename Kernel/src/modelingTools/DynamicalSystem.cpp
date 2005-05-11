#include "DynamicalSystem.h"

#include "LinearBC.h"
#include "NLinearBC.h"
#include "PeriodicBC.h"

#include "LinearDSIO.h"
#include "LagrangianDSIO.h"
#include "LagrangianLinearDSIO.h"

// --- Constructors ---


// From XML file
DynamicalSystem::DynamicalSystem(DSXML * dsXML):  DSType(NLDS), nsds(0), number(0), id("none"), n(0), x0(0), x(0), xDot(0), xFree(0), r(0),
  stepsInMemory(1), jacobianX(0), BC(0), dsxml(dsXML), vectorFieldPtr(0)
{
  IN("DynamicalSystem::DynamicalSystem - XML constructor\n");
  // --- get values in xml file ---
  if (dsXML != 0)
  {
    this->number = this->dsxml->getNumber();
    if (this->dsxml->hasId() == true) this->id = this->dsxml->getId();
    if (this->dsxml->hasN() == true) this->n = this->dsxml->getN();

    // --- Memory allocation for vector and matrix members ---
    this->x0 = new SimpleVector(n);
    this->x = new SimpleVector(n);
    this->xDot = new SimpleVector(n);
    this->xFree = new SimpleVector(n);
    this->r = new SimpleVector(this->n);
    this->jacobianX = new SiconosMatrix(this->n, this->n);

    // xml loading of vector and matrix members
    if (this->dsxml->hasX0() == true) *(this->x0) = this->dsxml->getX0();
    if (this->dsxml->hasX() == true) *(this->x) = this->dsxml->getX();
    if (this->dsxml->hasXMemory() == true) this->xMemory = SiconosMemory::SiconosMemory(this->dsxml->getXMemoryXML()); //this->dsxml->getXMemory();
    if (this->dsxml->hasXDot() == true) *(this->xDot) = this->dsxml->getXDot();
    if (this->dsxml->hasXDotMemory() == true) this->xDotMemory = SiconosMemory::SiconosMemory(this->dsxml->getXDotMemoryXML()); //this->dsxml->getXDotMemory();
    if (this->dsxml->hasStepsInMemory() == true) this->stepsInMemory = this->dsxml->getStepsInMemory();

    // --- Plugins ---
    string plugin;
    // VectorField
    if (this->dsxml->hasVectorFieldPlugin() == true)
    {
      plugin = this->dsxml->getVectorFieldPlugin();
      this->setVectorFieldFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
    }
    // JacobianX
    if (this->dsxml->hasComputeJacobianXPlugin() == true)
    {
      plugin = this->dsxml->getComputeJacobianXPlugin();
      this->setComputeJacobianXFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
    }
    // --- Boundary conditions ---
    fillBoundaryConditionsFromXml();

    // --- DS input-output ---
    fillDsioFromXml();
  }
  else
  {
    cout << "DynamicalSystem::DynamicalSystem - DSXML paramater must not be 0" << endl;
  }
  OUT("DynamicalSystem::DynamicalSystem - XML constructor\n");
}

// From a minimum set of data
DynamicalSystem::DynamicalSystem(int number, int n, SiconosVector* x0, string vectorFieldPlugin):
  DSType(NLDS), nsds(0), number(0), id("none"), n(n), x0(0), x(0), xDot(0), xFree(0), r(0),
  stepsInMemory(1), jacobianX(0), BC(0), dsxml(0), vectorFieldPtr(0)
{
  IN("DynamicalSystem::DynamicalSystem - Minimum data constructor\n");
  this->x0 = new SimpleVector(n);
  this->x = new SimpleVector(n);
  this->xDot = new SimpleVector(n);
  this->xFree = new SimpleVector(n);
  this->r = new SimpleVector(n);
  this->jacobianX = new SiconosMatrix(this->n, this->n);
  *(this->x0) = *x0;
  this->setVectorFieldFunction(this->cShared.getPluginName(vectorFieldPlugin), this->cShared.getPluginFunctionName(vectorFieldPlugin));
  this->setComputeJacobianXFunction("BasicPlugin.so", "computeJacobianX");
  OUT("DynamicalSystem::DynamicalSystem - Minimum data constructor\n");
}

// --- Destructor ---
DynamicalSystem::~DynamicalSystem()
{
  IN("DynamicalSystem::~DynamicalSystem()\n");
  delete x0;
  x0 = 0 ;
  delete x;
  x = 0;
  delete xDot;
  xDot = 0;
  delete xFree;
  xFree = 0;
  delete r;
  r = 0;
  delete jacobianX;
  jacobianX = 0;
  for (int i = 0; i < this->dsioVector.size(); i++)
  {
    delete this->dsioVector[i];
    dsioVector[i] = 0;
  }
  delete BC;
  BC = 0;
  OUT("DynamicalSystem::~DynamicalSystem()\n");
}

// Boundary conditions built-in (called from constructors)
void DynamicalSystem::fillBoundaryConditionsFromXml()
{
  IN("DynamicalSystem::fillBoundaryConditionsFromXml\n");
  if (this->dsxml->getBoundaryConditionXML() != 0)
  {
    if (this->dsxml->getBoundaryConditionXML()->getType() == LINEARBC_TAG)
    {
      //  Linear BC
      this->BC = new LinearBC();
      static_cast<LinearBC*>(this->BC)->createBoundaryCondition(this->dsxml->getBoundaryConditionXML());
    }
    else if (this->dsxml->getBoundaryConditionXML()->getType() == NON_LINEARBC_TAG)
    {
      // Non linear BC
      this->BC = new NLinearBC();
      static_cast<NLinearBC*>(this->BC)->createBoundaryCondition(this->dsxml->getBoundaryConditionXML());
    }

    else if (this->dsxml->getBoundaryConditionXML()->getType() == PERIODICBC_TAG)
    {
      // Periodic BC
      this->BC = new PeriodicBC();
      static_cast<PeriodicBC*>(this->BC)->createBoundaryCondition(this->dsxml->getBoundaryConditionXML());
    }
    else RuntimeException::selfThrow("DynamicalSystem::linkDSXML - bad kind of BoundaryCondition : " + this->dsxml->getBoundaryConditionXML()->getType());
  }
  OUT("DynamicalSystem::fillBoundaryConditionsFromXml\n");
}

// DSIO built-in (called from constructors)
void DynamicalSystem::fillDsioFromXml()
{
  IN("DynamicalSystem::fillDsioFromXml\n");
  DSInputOutput *dsio;
  // get the numbers of DSIO
  vector<int> nbDSIOtab = this->dsxml->getDSInputOutputNumbers();
  for (int i = 0; i < nbDSIOtab.size(); i++)
  {
    if (this->dsxml->getDSInputOutputXML(nbDSIOtab[i])->getType() == LINEAR_DSIO_TAG)
    {
      // Linear DSIO
      dsio = new LinearDSIO();
      this->dsioVector.push_back(dsio);
      static_cast<LinearDSIO*>(dsio)->createDSInputOutput(this->dsxml->getDSInputOutputXML(nbDSIOtab[i]));
    }
    else if (this->dsxml->getDSInputOutputXML(nbDSIOtab[i])->getType() == NON_LINEAR_DSIO_TAG)
    {
      // Non linear DSIO
      dsio = new DSInputOutput();
      this->dsioVector.push_back(dsio);
      static_cast<DSInputOutput*>(dsio)->createDSInputOutput(this->dsxml->getDSInputOutputXML(nbDSIOtab[i]));
    }
    else if (this->dsxml->getDSInputOutputXML(nbDSIOtab[i])->getType() == LAGRANGIAN_DSIO_TAG)
    {
      // Lagrangian DSIO
      dsio = new LagrangianDSIO();
      this->dsioVector.push_back(dsio);
      static_cast<LagrangianDSIO*>(dsio)->createDSInputOutput(this->dsxml->getDSInputOutputXML(nbDSIOtab[i]));
    }
    else if (this->dsxml->getDSInputOutputXML(nbDSIOtab[i])->getType() == LAGRANGIAN_LINEAR_DSIO_TAG)
    {
      // Linear lagrangian DSIO
      dsio = new LagrangianDSIO();
      this->dsioVector.push_back(dsio);
      static_cast<LagrangianLinearDSIO*>(dsio)->createDSInputOutput(this->dsxml->getDSInputOutputXML(nbDSIOtab[i]));
    }
    else RuntimeException::selfThrow("DynamicalSystem::linkDSXML - bad kind of DSInputOutput: " + this->dsxml->getDSInputOutputXML(nbDSIOtab[i])->getType());
  }
  OUT("DynamicalSystem::fillDsioFromXml\n");
}

DSInputOutput* DynamicalSystem::getDSInputOutput(int i)
{
  if (i < this->dsioVector.size())
  {
    return this->dsioVector[i];
  }
  else RuntimeException::selfThrow("DS - getDSInputOutput : \'i\' is out of range");
}

////////////////////////////

void DynamicalSystem::setVectorFieldFunction(const string& pluginPath, const string& functionName)
{
  this->vectorFieldPtr = 0;
  cShared.setFunction(&vectorFieldPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->vectorFieldFunctionName = plugin + ":" + functionName;
}

void DynamicalSystem::setComputeJacobianXFunction(const string& pluginPath, const string& functionName)
{
  this->computeJacobianXPtr = 0;
  cShared.setFunction(&computeJacobianXPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->computeJacobianXFunctionName = plugin + ":" + functionName;
}

void DynamicalSystem::vectorField(const double& time)
{
  if (vectorFieldPtr == 0)
    RuntimeException::selfThrow("vectorField() is not linked to a plugin function");

  int size = x->size();
  // const_cast to be deleted: problem const in C function signature? To see ...
  this->vectorFieldPtr(&size, const_cast<double*>(&time), &(*x)(0) , &(*xDot)(0));
}

void DynamicalSystem::computeJacobianX(const double& time)
{
  if (computeJacobianXPtr == 0)
    RuntimeException::selfThrow("computeJacobianX() is not linked to a plugin function");

  int size = x->size();
  // const_cast to be deleted: problem const in C function signature? To see ...
  this->computeJacobianXPtr(&size, const_cast<double*>(&time), &(*x)(0), &(*jacobianX)(0, 0));
}


void DynamicalSystem::swapInMemory(void)
{
  IN("DynamicalSystem::swapInMemory\n ");
  xMemory.swap(this->x);
  xDotMemory.swap(this->xDot);
  rMemory.swap(this->r);
  OUT("DynamicalSystem::swapInMemory\n ");
}

void DynamicalSystem::display() const
{
  IN("DynamicalSystem::display\n");
  cout << "____ data of the Dynamical System " << endl;
  cout << "| number : " << this->number << endl;
  cout << "| id : " << this->id << endl;
  cout << "| n : " << this->n << endl;
  cout << "| x " << endl;
  if (x != 0) this->x->display();
  else cout << "-> 0" << endl;
  cout << "| x0 " << endl;
  if (x0 != 0) this->x0->display();
  else cout << "-> 0" << endl;
  cout << "| xFree " << endl;
  if (xFree != 0) this->xFree->display();
  else cout << "-> 0" << endl;
  cout << "| xDot " << endl;
  if (xDot != 0) this->xDot->display();
  else cout << "-> 0" << endl;
  cout << "| stepsInMemory : " << this->stepsInMemory << endl;
  cout << "| r " << endl;
  if (r != 0) this->r->display();
  else cout << "-> 0" << endl;
  cout << "| VectorField plugin: " << this->vectorFieldFunctionName << endl;
  cout << "| JacobianX plugin: " << this->computeJacobianXFunctionName << endl;
  OUT("DynamicalSystem::display\n");

}

void DynamicalSystem::initMemory(const int& steps)
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

  // --- general DS data ---
  saveDSDataToXML();

  // --- other data ---
  if (this->dsxml != 0)
  {
    this->dsxml->setN(this->n);
    this->dsxml->setVectorFieldPlugin(this->vectorFieldFunctionName);
    if (this->computeJacobianXFunctionName != "")
      this->dsxml->setComputeJacobianXPlugin(this->computeJacobianXFunctionName);
    else
      this->dsxml->setComputeJacobianXPlugin("BasicPlugin:computeJacobianX");
  }
  else RuntimeException::selfThrow("DynamicalSystem::saveDSToXML - The DSXML object doesn't exists");
  OUT("DynamicalSystem::saveDSToXML\n");
}

// Save data common to each system into the xml file
void DynamicalSystem::saveDSDataToXML()
{
  //--- Boundary conditions ---
  saveBCToXML();
  // --- DS input-output ---
  saveDSIOToXML();
  if (this->dsxml != 0)
  {
    this->dsxml->setId(this->id);
    this->dsxml->setX0(this->x0);
    this->dsxml->setX(this->x);
    this->dsxml->setXMemory(&(this->xMemory));
    this->dsxml->setXDot(this->xDot);
    this->dsxml->setXDotMemory(&(this->xDotMemory));
    this->dsxml->setStepsInMemory(this->stepsInMemory);
    this->dsxml->setR(this->r);
  }
  else RuntimeException::selfThrow("DynamicalSystem::saveDSToXML - The DSXML object doesn't exists");
  OUT("DynamicalSystem::saveDSToXML\n");

}

// Save boundary conditions to xml file
void DynamicalSystem::saveBCToXML()
{
  if (this->BC != 0)
  {
    if (this->BC->getType() == LINEARBC)
      (static_cast<LinearBC*>(this->BC))->saveBCToXML();
    else if (this->BC->getType() == NLINEARBC)
      (static_cast<NLinearBC*>(this->BC))->saveBCToXML();
    else if (this->BC->getType() == PERIODICBC)
      (static_cast<PeriodicBC*>(this->BC))->saveBCToXML();
    else RuntimeException::selfThrow("DynamicalSystem::saveDSToXML - bad kind of BoundaryCondition");
  }
}

// Save DS Input-Output to xml file
void DynamicalSystem::saveDSIOToXML()
{
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
}

BoundaryCondition* DynamicalSystem::createPeriodicBC()
{
  this->BC = new PeriodicBC();
  static_cast<PeriodicBC*>(this->BC)->createBoundaryCondition(0);
  return this->BC;
}

BoundaryCondition* DynamicalSystem::createLinearBC(SiconosVector* omega, SiconosMatrix* omega0, SiconosMatrix* omegaT)
{
  this->BC = new LinearBC();
  static_cast<LinearBC*>(this->BC)->createBoundaryCondition(0, omega, omega0, omegaT);
  return this->BC;
}

BoundaryCondition* DynamicalSystem::createNLinearBC()
{
  this->BC = new NLinearBC();
  static_cast<NLinearBC*>(this->BC)->createBoundaryCondition(0);
  return this->BC;
}

double DynamicalSystem::dsConvergenceIndicator() const
{
  RuntimeException::selfThrow("DynamicalSystem:dsConvergenceIndicator - not yet implemented for Dynamical system type :" + this->DSType);
}

// Default constructor
DynamicalSystem::DynamicalSystem(): DSType(NLDS), nsds(0), number(0), id("none"), n(0), x0(0), x(0), xDot(0), xFree(0), r(0),
  stepsInMemory(1), jacobianX(0), BC(0), dsxml(0), vectorFieldPtr(0)
{
  IN("DynamicalSystem::DynamicalSystem - Default constructor\n");
  // --- plugins -> connected to  "false" plugin
  this->setVectorFieldFunction("BasicPlugin.so", "vectorField");
  this->setComputeJacobianXFunction("BasicPlugin.so", "computeJacobianX");
  OUT("DynamicalSystem::DynamicalSystem - Default constructor\n");
}
