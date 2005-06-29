#include "DynamicalSystem.h"

// includes to be deleted thanks to factories
#include "LinearBC.h"
#include "NLinearBC.h"
#include "PeriodicBC.h"
#include "LinearDSIO.h"
#include "LagrangianDSIO.h"
#include "LagrangianLinearDSIO.h"

using namespace std;

// --- Constructors ---


// From XML file
DynamicalSystem::DynamicalSystem(DSXML * dsXML):
  DSType(NLDS), nsds(NULL), number(0), id("none"), n(0), x0(NULL), x(NULL), xDot(NULL), xFree(NULL), r(NULL),
  stepsInMemory(1), jacobianX(NULL), BC(NULL), dsxml(dsXML), vectorFieldPtr(NULL),
  isX0AllocatedIn(true), isXAllocatedIn(true), isXMemoryAllocatedIn(false), isXDotAllocatedIn(true), isXDotMemoryAllocatedIn(false),
  isXFreeAllocatedIn(true), isRAllocatedIn(true), isRMemoryAllocatedIn(false), isJacobianXAllocatedIn(true), isBCAllocatedIn(false)
{
  IN("DynamicalSystem::DynamicalSystem - XML constructor\n");
  // --- get values in xml file ---
  if (dsXML != NULL)
  {
    number = dsxml->getNumber();
    if (dsxml->hasId() == true) id = dsxml->getId();
    if (dsxml->hasN() == true) n = dsxml->getN();

    // --- Memory allocation for vector and matrix members ---
    x0 = new SimpleVector(n);
    x = new SimpleVector(n);
    xDot = new SimpleVector(n);
    xFree = new SimpleVector(n);
    r = new SimpleVector(n);
    jacobianX = new SiconosMatrix(n, n);

    // xml loading of vector and matrix members
    if (dsxml->hasX0() == true) *(x0) = dsxml->getX0();
    if (dsxml->hasX() == true) *(x) = dsxml->getX();
    if (dsxml->hasXDot() == true) *(xDot) = dsxml->getXDot();
    if (dsxml->hasStepsInMemory() == true) stepsInMemory = dsxml->getStepsInMemory();
    if (dsxml->hasXMemory() == true)
    {
      xMemory = new SiconosMemory(dsxml->getXMemoryXML());
      isXMemoryAllocatedIn = true;
    }
    if (dsxml->hasXDotMemory() == true)
    {
      xDotMemory = new SiconosMemory(dsxml->getXDotMemoryXML());
      isXDotMemoryAllocatedIn = true;
    }

    // --- Plugins ---
    string plugin;
    // VectorField
    if (dsxml->hasVectorFieldPlugin() == true)
    {
      plugin = dsxml->getVectorFieldPlugin();
      setVectorFieldFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    // JacobianX
    if (dsxml->hasComputeJacobianXPlugin() == true)
    {
      plugin = dsxml->getComputeJacobianXPlugin();
      setComputeJacobianXFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    // --- Boundary conditions ---
    fillBoundaryConditionsFromXml();

    // --- DS input-output ---
    fillDsioFromXml();
  }
  else
    cout << "DynamicalSystem::DynamicalSystem - DSXML paramater must not be NULL" << endl;
  OUT("DynamicalSystem::DynamicalSystem - XML constructor\n");
}

// From a minimum set of data
DynamicalSystem::DynamicalSystem(int newNumber, int newN, SiconosVector* newX0, string vectorFieldPlugin):
  DSType(NLDS), nsds(NULL), number(newNumber), id("none"), n(newN), x0(NULL), x(NULL), xDot(NULL), xFree(NULL), r(NULL),
  stepsInMemory(1), jacobianX(NULL), BC(NULL), dsxml(NULL), vectorFieldPtr(NULL),
  isX0AllocatedIn(true), isXAllocatedIn(true), isXMemoryAllocatedIn(false), isXDotAllocatedIn(true), isXDotMemoryAllocatedIn(false),
  isXFreeAllocatedIn(true), isRAllocatedIn(true), isRMemoryAllocatedIn(false), isJacobianXAllocatedIn(true), isBCAllocatedIn(false)
{
  IN("DynamicalSystem::DynamicalSystem - Minimum data constructor\n");
  x0 = new SimpleVector(n);
  x = new SimpleVector(n);
  xDot = new SimpleVector(n);
  xFree = new SimpleVector(n);
  r = new SimpleVector(n);
  jacobianX = new SiconosMatrix(n, n);
  *(x0) = *newX0;
  setVectorFieldFunction(cShared.getPluginName(vectorFieldPlugin), cShared.getPluginFunctionName(vectorFieldPlugin));
  setComputeJacobianXFunction("BasicPlugin.so", "computeJacobianX");
  OUT("DynamicalSystem::DynamicalSystem - Minimum data constructor\n");
}

// --- Destructor ---
DynamicalSystem::~DynamicalSystem()
{
  IN("DynamicalSystem::~DynamicalSystem()\n");
  if (isX0AllocatedIn)
  {
    delete x0;
    x0 = NULL ;
  }
  if (isXAllocatedIn)
  {
    delete x;
    x = NULL;
  }
  if (isXMemoryAllocatedIn)
  {
    delete xMemory;
    xMemory = NULL;
  }
  if (isXDotAllocatedIn)
  {
    delete xDot;
    xDot = NULL;
  }
  if (isXDotMemoryAllocatedIn)
  {
    delete xDotMemory;
    xDotMemory = NULL;
  }
  if (isXFreeAllocatedIn)
  {
    delete xFree;
    xFree = NULL;
  }
  if (isRAllocatedIn)
  {
    delete r;
    r = NULL;
  }
  if (isRMemoryAllocatedIn)
  {
    delete rMemory;
    rMemory = NULL;
  }
  if (isJacobianXAllocatedIn)
  {
    delete jacobianX;
    jacobianX = NULL;
  }
  for (unsigned int i = 0; i < dsioVector.size(); i++)
  {
    if (isDsioAllocatedIn[i])
    {
      delete dsioVector[i];
      dsioVector[i] = NULL;
    }
  }
  if (isBCAllocatedIn)
  {
    delete BC;
    BC = NULL;
  }
  OUT("DynamicalSystem::~DynamicalSystem()\n");
}

// Boundary conditions built-in (called from constructors)
void DynamicalSystem::fillBoundaryConditionsFromXml()
{
  IN("DynamicalSystem::fillBoundaryConditionsFromXml\n");
  if (dsxml->getBoundaryConditionXML() != 0)
  {
    if (dsxml->getBoundaryConditionXML()->getType() == LINEARBC_TAG)
    {
      //  Linear BC
      BC = new LinearBC();
      static_cast<LinearBC*>(BC)->createBoundaryCondition(dsxml->getBoundaryConditionXML());
      isBCAllocatedIn = true;
    }
    else if (dsxml->getBoundaryConditionXML()->getType() == NON_LINEARBC_TAG)
    {
      // Non linear BC
      BC = new NLinearBC();
      static_cast<NLinearBC*>(BC)->createBoundaryCondition(dsxml->getBoundaryConditionXML());
      isBCAllocatedIn = true;
    }

    else if (dsxml->getBoundaryConditionXML()->getType() == PERIODICBC_TAG)
    {
      // Periodic BC
      BC = new PeriodicBC();
      static_cast<PeriodicBC*>(BC)->createBoundaryCondition(dsxml->getBoundaryConditionXML());
      isBCAllocatedIn = true;
    }
    else RuntimeException::selfThrow("DynamicalSystem::linkDSXML - bad kind of BoundaryCondition : " + dsxml->getBoundaryConditionXML()->getType());
  }
  OUT("DynamicalSystem::fillBoundaryConditionsFromXml\n");
}

// DSIO built-in (called from constructors)
void DynamicalSystem::fillDsioFromXml()
{
  IN("DynamicalSystem::fillDsioFromXml\n");
  DSInputOutput *dsio;
  // get the numbers of DSIO
  vector<int> nbDSIOtab = dsxml->getDSInputOutputNumbers();
  for (unsigned int i = 0; i < nbDSIOtab.size(); i++)
  {
    if (dsxml->getDSInputOutputXML(nbDSIOtab[i])->getType() == LINEAR_DSIO_TAG)
    {
      // Linear DSIO
      dsio = new LinearDSIO();
      dsioVector.push_back(dsio);
      isDsioAllocatedIn.push_back(true);
      static_cast<LinearDSIO*>(dsio)->createDSInputOutput(dsxml->getDSInputOutputXML(nbDSIOtab[i]));
    }
    else if (dsxml->getDSInputOutputXML(nbDSIOtab[i])->getType() == NON_LINEAR_DSIO_TAG)
    {
      // Non linear DSIO
      dsio = new DSInputOutput();
      dsioVector.push_back(dsio);
      isDsioAllocatedIn.push_back(true);
      static_cast<DSInputOutput*>(dsio)->createDSInputOutput(dsxml->getDSInputOutputXML(nbDSIOtab[i]));
    }
    else if (dsxml->getDSInputOutputXML(nbDSIOtab[i])->getType() == LAGRANGIAN_DSIO_TAG)
    {
      // Lagrangian DSIO
      dsio = new LagrangianDSIO();
      dsioVector.push_back(dsio);
      isDsioAllocatedIn.push_back(true);
      static_cast<LagrangianDSIO*>(dsio)->createDSInputOutput(dsxml->getDSInputOutputXML(nbDSIOtab[i]));
    }
    else if (dsxml->getDSInputOutputXML(nbDSIOtab[i])->getType() == LAGRANGIAN_LINEAR_DSIO_TAG)
    {
      // Linear lagrangian DSIO
      dsio = new LagrangianDSIO();
      dsioVector.push_back(dsio);
      isDsioAllocatedIn.push_back(true);
      static_cast<LagrangianLinearDSIO*>(dsio)->createDSInputOutput(dsxml->getDSInputOutputXML(nbDSIOtab[i]));
    }
    else RuntimeException::selfThrow("DynamicalSystem::linkDSXML - bad kind of DSInputOutput: " + dsxml->getDSInputOutputXML(nbDSIOtab[i])->getType());
  }
  OUT("DynamicalSystem::fillDsioFromXml\n");
}


void DynamicalSystem::setX0Ptr(SiconosVector* newPtr)
{
  if (isX0AllocatedIn) delete x0;
  x0 = newPtr;
  isX0AllocatedIn = false;
}

void DynamicalSystem::setXPtr(SiconosVector* newPtr)
{
  if (isXAllocatedIn) delete x;
  x = newPtr;
  isXAllocatedIn = false;
}

void DynamicalSystem::setXMemoryPtr(SiconosMemory * newPtr)
{
  if (isXMemoryAllocatedIn) delete xMemory;
  xMemory = newPtr;
  isXMemoryAllocatedIn = false;
}

void DynamicalSystem::setXDotPtr(SiconosVector* newPtr)
{
  if (isXDotAllocatedIn) delete xDot;
  xDot = newPtr;
  isXDotAllocatedIn = false;
}

void DynamicalSystem::setXDotMemoryPtr(SiconosMemory * newPtr)
{
  if (isXDotMemoryAllocatedIn) delete xDotMemory;
  xDotMemory = newPtr;
  isXDotMemoryAllocatedIn = false;
}

void DynamicalSystem::setXFreePtr(SiconosVector* newPtr)
{
  if (isXFreeAllocatedIn) delete xFree;
  xFree = newPtr;
  isXFreeAllocatedIn = false;
}

void DynamicalSystem::setRPtr(SimpleVector *newPtr)
{
  if (isRAllocatedIn) delete r;
  r = newPtr;
  isRAllocatedIn = false;
}

void DynamicalSystem::setRMemoryPtr(SiconosMemory * newPtr)
{
  if (isRMemoryAllocatedIn) delete rMemory;
  rMemory = newPtr;
  isRMemoryAllocatedIn = false;
}

void DynamicalSystem::setJacobianXPtr(SiconosMatrix *newPtr)
{
  if (isJacobianXAllocatedIn) delete jacobianX;
  jacobianX = newPtr;
  isJacobianXAllocatedIn = false;
}

DSInputOutput* DynamicalSystem::getDSInputOutput(const unsigned int& i)
{
  if (i >= dsioVector.size())
    RuntimeException::selfThrow("DS - getDSInputOutput : \'i\' is out of range");
  return dsioVector[i];
}

void DynamicalSystem::setBoundaryConditionPtr(BoundaryCondition *newBC)
{
  if (isBCAllocatedIn) delete BC;
  BC = newBC;
  isBCAllocatedIn = false;
}

void DynamicalSystem::setVectorFieldFunction(const string& pluginPath, const string& functionName)
{
  vectorFieldPtr = NULL;
  cShared.setFunction(&vectorFieldPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  vectorFieldFunctionName = plugin + ":" + functionName;
}

void DynamicalSystem::setComputeJacobianXFunction(const string& pluginPath, const string& functionName)
{
  computeJacobianXPtr = NULL;
  cShared.setFunction(&computeJacobianXPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  computeJacobianXFunctionName = plugin + ":" + functionName;
}

void DynamicalSystem::vectorField(const double& time)
{
  if (vectorFieldPtr == NULL)
    RuntimeException::selfThrow("vectorField() is not linked to a plugin function");

  int size = x->size();
  // const_cast to be deleted: problem const in C function signature? To see ...
  vectorFieldPtr(&size, const_cast<double*>(&time), &(*x)(0) , &(*xDot)(0));
}

void DynamicalSystem::computeJacobianX(const double& time)
{
  if (computeJacobianXPtr == NULL)
    RuntimeException::selfThrow("computeJacobianX() is not linked to a plugin function");

  int size = x->size();
  // const_cast to be deleted: problem const in C function signature? To see ...
  computeJacobianXPtr(&size, const_cast<double*>(&time), &(*x)(0), &(*jacobianX)(0, 0));
}


void DynamicalSystem::swapInMemory(void)
{
  IN("DynamicalSystem::swapInMemory\n ");
  xMemory->swap(*x);
  xDotMemory->swap(*xDot);
  rMemory->swap(*r);
  OUT("DynamicalSystem::swapInMemory\n ");
}

void DynamicalSystem::display() const
{
  IN("DynamicalSystem::display\n");
  cout << "____ data of the Dynamical System " << endl;
  cout << "| number : " << number << endl;
  cout << "| id : " << id << endl;
  cout << "| n : " << n << endl;
  cout << "| x " << endl;
  if (x != NULL) x->display();
  else cout << "-> NULL" << endl;
  cout << "| x0 " << endl;
  if (x0 != NULL) x0->display();
  else cout << "-> NULL" << endl;
  cout << "| xFree " << endl;
  if (xFree != NULL) xFree->display();
  else cout << "-> NULL" << endl;
  cout << "| xDot " << endl;
  if (xDot != NULL) xDot->display();
  else cout << "-> NULL" << endl;
  cout << "| stepsInMemory : " << stepsInMemory << endl;
  cout << "| r " << endl;
  if (r != NULL) r->display();
  else cout << "-> NULL" << endl;
  cout << "| VectorField plugin: " << vectorFieldFunctionName << endl;
  cout << "| JacobianX plugin: " << computeJacobianXFunctionName << endl;
  OUT("DynamicalSystem::display\n");

}

void DynamicalSystem::initMemory(const unsigned int& steps)
{
  IN("DynamicalSystem::initMemory\n");
  if (steps == 0)
    cout << "Warning : DynamicalSystem::initMemory with size equal to zero" << endl;
  else
  {
    // \warning FP: explicit call to xml; what about other cases (without xml)? To be reviewed
    stepsInMemory = steps;
    if (isXMemoryAllocatedIn) delete rMemory;
    //xMemory = new SiconosMemory(steps, xMemory.getSiconosMemoryXML() );
    xMemory = new SiconosMemory(steps);
    isXMemoryAllocatedIn = true;
    if (isXDotMemoryAllocatedIn) delete rMemory;
    //xDotMemory = new SiconosMemory(steps, xDotMemory.getSiconosMemoryXML());
    xDotMemory = new SiconosMemory(steps);
    isXDotMemoryAllocatedIn = true;
    if (isRMemoryAllocatedIn) delete rMemory;
    rMemory = new SiconosMemory(steps);
    isRMemoryAllocatedIn = true;
  }

  OUT("DynamicalSystem::initMemory\n");
}


void DynamicalSystem::saveDSToXML()
{
  IN("DynamicalSystem::saveDSToXML\n");

  // --- general DS data ---
  saveDSDataToXML();

  // --- other data ---
  if (dsxml != NULL)
  {
    dsxml->setN(n);
    dsxml->setVectorFieldPlugin(vectorFieldFunctionName);
    if (computeJacobianXFunctionName != "")
      dsxml->setComputeJacobianXPlugin(computeJacobianXFunctionName);
    else
      dsxml->setComputeJacobianXPlugin("BasicPlugin:computeJacobianX");
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
  if (dsxml != NULL)
  {
    dsxml->setId(id);
    dsxml->setX0(x0);
    dsxml->setX(x);
    dsxml->setXMemory(xMemory);
    dsxml->setXDot(xDot);
    dsxml->setXDotMemory(xDotMemory);
    dsxml->setStepsInMemory(stepsInMemory);
    dsxml->setR(r);
  }
  else RuntimeException::selfThrow("DynamicalSystem::saveDSToXML - The DSXML object doesn't exists");
  OUT("DynamicalSystem::saveDSToXML\n");

}

// Save boundary conditions to xml file
void DynamicalSystem::saveBCToXML()
{
  if (BC != NULL)
  {
    if (BC->getType() == LINEARBC)
      (static_cast<LinearBC*>(BC))->saveBCToXML();
    else if (BC->getType() == NLINEARBC)
      (static_cast<NLinearBC*>(BC))->saveBCToXML();
    else if (BC->getType() == PERIODICBC)
      (static_cast<PeriodicBC*>(BC))->saveBCToXML();
    else RuntimeException::selfThrow("DynamicalSystem::saveDSToXML - bad kind of BoundaryCondition");
  }
}

// Save DS Input-Output to xml file
void DynamicalSystem::saveDSIOToXML()
{
  if (dsioVector.size() != 0)
  {
    for (unsigned int i = 0; i < dsioVector.size(); i++)
    {
      if (dsioVector[i]->getType() == LINEARDSIO)
        (static_cast<LinearDSIO*>(dsioVector[i]))->saveDSInputOutputToXML();
      else if (dsioVector[i]->getType() == NLINEARDSIO)
        (static_cast<DSInputOutput*>(dsioVector[i]))->saveDSInputOutputToXML();
      else if (dsioVector[i]->getType() == LAGRANGIANDSIO)
        (static_cast<LagrangianDSIO*>(dsioVector[i]))->saveDSInputOutputToXML();
      else RuntimeException::selfThrow("DynamicalSystem::saveDSToXML - bad kind of DSInputOuput");
    }
  }
}

BoundaryCondition* DynamicalSystem::createPeriodicBC()
{
  BC = new PeriodicBC();
  isBCAllocatedIn = true;
  static_cast<PeriodicBC*>(BC)->createBoundaryCondition(0);
  return BC;
}

BoundaryCondition* DynamicalSystem::createLinearBC(SiconosVector* omega, SiconosMatrix* omega0, SiconosMatrix* omegaT)
{
  BC = new LinearBC();
  isBCAllocatedIn = true;
  static_cast<LinearBC*>(BC)->createBoundaryCondition(0, omega, omega0, omegaT);
  return BC;
}

BoundaryCondition* DynamicalSystem::createNLinearBC()
{
  BC = new NLinearBC();
  isBCAllocatedIn = true;
  static_cast<NLinearBC*>(BC)->createBoundaryCondition(0);
  return BC;
}

double DynamicalSystem::dsConvergenceIndicator()
{
  RuntimeException::selfThrow("DynamicalSystem:dsConvergenceIndicator - not yet implemented for Dynamical system type :" + DSType);
  return 0;
}

// Default constructor
DynamicalSystem::DynamicalSystem():
  DSType(NLDS), nsds(NULL), number(0), id("none"), n(0), x0(NULL), x(NULL), xDot(NULL), xFree(NULL), r(NULL),
  stepsInMemory(1), jacobianX(NULL), BC(NULL), dsxml(NULL), vectorFieldPtr(NULL), isX0AllocatedIn(false),
  isXAllocatedIn(false), isXMemoryAllocatedIn(false), isXDotAllocatedIn(false), isXDotMemoryAllocatedIn(false),
  isXFreeAllocatedIn(false), isRAllocatedIn(false), isRMemoryAllocatedIn(false), isJacobianXAllocatedIn(false),
  isBCAllocatedIn(false)

{
  // --- plugins -> connected to  "false" plugin
  setVectorFieldFunction("BasicPlugin.so", "vectorField");
  setComputeJacobianXFunction("BasicPlugin.so", "computeJacobianX");
}
