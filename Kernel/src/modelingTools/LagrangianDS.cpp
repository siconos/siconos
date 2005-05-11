#include "LagrangianDS.h"
#include "check.h"
#include "LinearBC.h"
#include "NLinearBC.h"
#include "PeriodicBC.h"
#include "LinearDSIO.h"
#include "LagrangianDSIO.h"
#include "LagrangianLinearDSIO.h"


// --- Constructor from an xml file ---
LagrangianDS::LagrangianDS(DSXML * dsXML): DynamicalSystem(), ndof(0), q(0), q0(0), qFree(0), velocity(0), velocity0(0),
  velocityFree(0), p(0), mass(0), fInt(0), fExt(0), QNLInertia(0), jacobianQFInt(0),
  jacobianVelocityFInt(0), jacobianQQNLInertia(0), jacobianVelocityQNLInertia(0),
  computeMassPtr(0), computeFIntPtr(0), computeFExtPtr(0), computeQNLInertiaPtr(0), computeJacobianQFIntPtr(0),
  computeJacobianVelocityFIntPtr(0), computeJacobianQQNLInertiaPtr(0), computeJacobianVelocityQNLInertiaPtr(0)
{
  IN("LagrangianDS::LagrangianDS() - Xml constructor\n");
  if (dsXML != 0)
  {
    // --- DS BASE-CLASS MEMBERS ---
    // --- Settings and xml load ---
    this->DSType = LNLDS;
    this->dsxml = dsXML;
    this->number = this->dsxml->getNumber();
    if (this->dsxml->hasId() == true) this->id = this->dsxml->getId();
    // -- Memory allocation for vector and matrix members --
    this->x = new CompositeVector();
    this->x0 = new CompositeVector();
    this->xDot = new CompositeVector();
    this->xFree = new CompositeVector();
    //
    this->r = new SimpleVector(this->n);
    // this->jacobianX = new SiconosMatrix(this->n,this->n);
    if (this->dsxml->hasStepsInMemory() == true) this->stepsInMemory = this->dsxml->getStepsInMemory();

    // -- plugins --
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
    // -- Boundary conditions --
    fillBoundaryConditionsFromXml();
    // -- DS input-output --
    fillDsioFromXml();

    // --- LAGRANGIAN INHERITED CLASS MEMBERS ---
    // --- Settings and xml load ---
    // -- Lagrangian  xml object --
    LagrangianDSXML* lgptr = static_cast <LagrangianDSXML*>(this->dsxml);

    // -- Size of the system and number of degrees of freedom --
    this->ndof = lgptr->getNdof();
    this->n = 2 * this->ndof;

    // -- Vector and matrix members memory allocation --
    this->q = new SimpleVector(this->ndof);
    this->q0 = new SimpleVector(this->ndof);
    this->qFree = new SimpleVector(this->ndof);
    this->velocity = new SimpleVector(this->ndof);
    this->velocity0 = new SimpleVector(this->ndof);
    this->velocityFree = new SimpleVector(this->ndof);
    this->p = new SimpleVector(this->ndof);
    this->mass = new SiconosMatrix(this->ndof, this->ndof);
    this->fInt = new SimpleVector(this->ndof);
    this->fExt = new SimpleVector(this->ndof);
    this->QNLInertia = new SimpleVector(this->ndof);
    this->jacobianQFInt = new SiconosMatrix(this->ndof, this->ndof);
    this->jacobianVelocityFInt = new SiconosMatrix(this->ndof, this->ndof);
    this->jacobianQQNLInertia = new SiconosMatrix(this->ndof, this->ndof);
    this->jacobianVelocityQNLInertia = new SiconosMatrix(this->ndof, this->ndof);

    // -- Plugins --
    // Mass
    if ((static_cast <LagrangianDSXML*>(this->dsxml))->isMPlugin())
    {
      plugin = lgptr->getMPlugin();
      this->setComputeMassFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
    }
    // \warning : VA:  It is a very good idea to take the constant Mass Matrix, but for the moment a constant
    //  Mass Matrix is only read by a LagrangianLinearTIDS
    else *this->mass = lgptr->getMMatrix();

    // q0, q and qMemory
    *this->q0 = lgptr->getQ0();
    if (lgptr->hasQ()) *this->q = lgptr->getQ();
    else *this->q = *this->q0;
    if (lgptr->hasQMemory()) this->qMemory = SiconosMemory::SiconosMemory(lgptr->getQMemoryXML());

    // velocity0, velocity and velocityMemory
    *this->velocity0 = lgptr->getVelocity0();
    if (lgptr->hasVelocity()) *this->velocity = lgptr->getVelocity();
    else *this->velocity = *this->velocity0;
    if (lgptr->hasVelocityMemory()) this->velocityMemory = SiconosMemory::SiconosMemory(lgptr->getVelocityMemoryXML());

    // fill in x, x0 and xFree with q and velocity
    static_cast<CompositeVector*>(this->x)->add(this->q);
    static_cast<CompositeVector*>(this->x)->add(this->velocity);
    static_cast<CompositeVector*>(this->x0)->add(this->q0);
    static_cast<CompositeVector*>(this->x0)->add(this->velocity0);
    static_cast<CompositeVector*>(this->xFree)->add(this->qFree);
    static_cast<CompositeVector*>(this->xFree)->add(this->velocityFree);

    // FExt
    if (lgptr->isFextPlugin())
    {
      plugin = lgptr->getFextPlugin();
      this->setComputeFExtFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
    }
    else *this->fExt = lgptr->getFextVector();

    // -- FInt and its jacobian --
    if (lgptr->hasFint())
    {
      if (lgptr->isFintPlugin())
      {
        plugin = lgptr->getFintPlugin();
        this->setComputeFIntFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
      }
      else *this->fInt = lgptr->getFintVector();
    }
    // Jacobian Q FInt
    if (lgptr->isJacobianQFintPlugin())
    {
      plugin = lgptr->getJacobianQFintPlugin();
      this->setComputeJacobianQFIntFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
    }
    else *this->jacobianQFInt = lgptr->getJacobianQFintMatrix();
    // Jacobian Velocity FInt
    if (lgptr->isJacobianVelocityFintPlugin())
    {
      plugin = lgptr->getJacobianVelocityFintPlugin();
      this->setComputeJacobianVelocityFIntFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
    }
    else *this->jacobianVelocityFInt = lgptr->getJacobianVelocityFintMatrix();

    // -- QNL Inertia and its jacobian --
    if (lgptr->isQNLInertiaPlugin())
    {
      plugin = lgptr->getQNLInertiaPlugin();
      this->setComputeQNLInertiaFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
    }
    else *this->QNLInertia = lgptr->getQNLInertiaVector();
    // Jacobian Q QNLInertia
    if (lgptr->isJacobianQQNLInertiaPlugin())
    {
      plugin = lgptr->getJacobianQQNLInertiaPlugin();
      this->setComputeJacobianQQNLInertiaFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
    }
    else *this->jacobianQQNLInertia = lgptr->getJacobianQQNLInertiaMatrix();
    // Jacobian Velocity QNLInertia
    if (lgptr->isJacobianVelocityQNLInertiaPlugin())
    {
      plugin = lgptr->getJacobianVelocityQNLInertiaPlugin();
      this->setComputeJacobianVelocityQNLInertiaFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
    }
    else *this->jacobianVelocityQNLInertia = lgptr->getJacobianVelocityQNLInertiaMatrix();

  }
  else RuntimeException::selfThrow("LagrangianDS::LagrangianDS - DSXML paramater must not be 0");
  OUT("LagrangianDS::LagrangianDS() - Xml constructor\n");
}

// From a minimum set of data
LagrangianDS::LagrangianDS(int newNumber, int newNdof,
                           SiconosVector* newQ0, SiconosVector* newVelocity0, string mass,
                           string fInt, string fExt, string jacobianQFInt, string jacobianVelocityFInt,
                           string jacobianQQNLInertia, string jacobianVelocityQNLInertia, string QNLInertia):
  DynamicalSystem(), ndof(newNdof), q(0), q0(0), qFree(0), velocity(0), velocity0(0), velocityFree(0), p(0), mass(0),
  fInt(0), fExt(0), QNLInertia(0), jacobianQFInt(0), jacobianVelocityFInt(0), jacobianQQNLInertia(0), jacobianVelocityQNLInertia(0),
  computeMassPtr(0), computeFIntPtr(0), computeFExtPtr(0), computeQNLInertiaPtr(0), computeJacobianQFIntPtr(0),
  computeJacobianVelocityFIntPtr(0), computeJacobianQQNLInertiaPtr(0), computeJacobianVelocityQNLInertiaPtr(0)

{
  IN("LagrangianDS::LagrangianDS - From a minimum set of data\n");
  // --- DS BASE-CLASS MEMBERS ---
  // --- Settings and xml load ---
  this->DSType = LNLDS;
  this->number = newNumber;
  this->n = 2 * ndof;

  // -- Memory allocation for vector and matrix members --
  this->x = new CompositeVector();
  this->x0 = new CompositeVector();
  this->xDot = new CompositeVector();
  this->xFree = new CompositeVector();
  //
  // this->r = new SimpleVector(this->n);
  // this->jacobianX = new SiconosMatrix(this->n,this->n);
  // -- plugins --
  string plugin;
  // VectorField
  // JacobianX

  // --- LAGRANGIAN INHERITED CLASS MEMBERS ---
  // -- Memory allocation for vector and matrix members --
  this->mass = new SiconosMatrix(this->ndof, this->ndof);
  this->q = new SimpleVector(this->ndof);
  this->q0 = new SimpleVector(this->ndof);
  this->qFree = new SimpleVector(this->ndof);
  this->velocity = new SimpleVector(this->ndof);
  this->velocity0 = new SimpleVector(this->ndof);
  this->velocityFree = new SimpleVector(this->ndof);
  this->p = new SimpleVector(this->ndof);
  this->fInt = new SimpleVector(this->ndof);
  this->fExt = new SimpleVector(this->ndof);
  this->QNLInertia = new SimpleVector(this->ndof);
  this->jacobianQFInt = new SiconosMatrix(this->ndof, this->ndof);
  this->jacobianVelocityFInt = new SiconosMatrix(this->ndof, this->ndof);
  this->jacobianQQNLInertia = new SiconosMatrix(this->ndof, this->ndof);
  this->jacobianVelocityQNLInertia = new SiconosMatrix(this->ndof, this->ndof);

  // --- initial state filling ---
  *this->q0 = *newQ0;
  *this->q = *newQ0;
  *this->velocity0 = *newVelocity0;
  *this->velocity = *newVelocity0;

  // --- x, xDot and xFree update ---

  static_cast<CompositeVector*>(this->x)->add(this->q);
  static_cast<CompositeVector*>(this->x)->add(this->velocity);
  static_cast<CompositeVector*>(this->x0)->add(this->q0);
  static_cast<CompositeVector*>(this->x0)->add(this->velocity0);
  static_cast<CompositeVector*>(this->xFree)->add(this->qFree);
  static_cast<CompositeVector*>(this->xFree)->add(this->velocityFree);

  //   --- plugins ---
  this->setComputeMassFunction(this->cShared.getPluginName(mass), this->cShared.getPluginFunctionName(mass));
  this->setComputeFIntFunction(this->cShared.getPluginName(fInt), this->cShared.getPluginFunctionName(fInt));
  this->setComputeFExtFunction(this->cShared.getPluginName(fExt), this->cShared.getPluginFunctionName(fExt));
  this->setComputeJacobianQFIntFunction(this->cShared.getPluginName(jacobianQFInt), this->cShared.getPluginFunctionName(jacobianQFInt));
  this->setComputeJacobianVelocityFIntFunction(this->cShared.getPluginName(jacobianVelocityFInt), this->cShared.getPluginFunctionName(jacobianQFInt));
  this->setComputeQNLInertiaFunction(this->cShared.getPluginName(QNLInertia), this->cShared.getPluginFunctionName(QNLInertia));
  this->setComputeJacobianQQNLInertiaFunction(this->cShared.getPluginName(jacobianQQNLInertia), this->cShared.getPluginFunctionName(jacobianQQNLInertia));
  this->setComputeJacobianVelocityQNLInertiaFunction(this->cShared.getPluginName(jacobianVelocityQNLInertia), this->cShared.getPluginFunctionName(jacobianVelocityQNLInertia));
  OUT("LagrangianDS::LagrangianDS - From a minimum set of data\n");
}

LagrangianDS::~LagrangianDS()
{
  IN("LagrangianDS::~LagrangianDS()\n");
  delete q;
  q = 0;
  delete q0;
  q0 = 0;
  delete qFree;
  qFree = 0;
  delete velocity ;
  velocity = 0;
  delete velocity0 ;
  velocity0 = 0;
  delete velocityFree ;
  velocityFree = 0;
  delete p ;
  p = 0;
  delete mass;
  mass = 0;
  delete fInt ;
  fInt = 0;
  delete fExt ;
  fExt = 0;
  delete QNLInertia ;
  QNLInertia = 0;
  delete jacobianQFInt  ;
  jacobianQFInt = 0;
  delete  jacobianVelocityFInt ;
  jacobianVelocityFInt = 0;
  delete jacobianQQNLInertia ;
  jacobianQQNLInertia = 0;
  delete jacobianVelocityQNLInertia ;
  jacobianVelocityQNLInertia = 0;
  OUT("LagrangianDS::~LagrangianDS()\n");
}


// --- Plugins related functions ---
void LagrangianDS::computeMass(double time)
{
  if (computeMassPtr == 0)
    RuntimeException::selfThrow("computeMass() is not linked to a plugin function");
  int size = q->size();
  this->computeMassPtr(&size, &time, &(*q)(0), &(*mass)(0, 0));
}

void LagrangianDS::computeMass(double time, SimpleVector *q)
{
  if (computeMassPtr == 0)
    RuntimeException::selfThrow("computeMass() is not linked to a plugin function");

  int size = q->size();
  this->computeMassPtr(&size, &time, &(*q)(0), &(*mass)(0, 0));
}

void LagrangianDS::computeFInt(double time)
{
  if (computeFIntPtr == 0)
    RuntimeException::selfThrow("computeFInt() is not linked to a plugin function");

  int size = this->q->size();
  this->computeFIntPtr(&size, &time, &(*this->q)(0), &(*this->velocity)(0), &(*this->fInt)(0));
}
void LagrangianDS::computeFInt(double time, SimpleVector *q, SimpleVector *velocity)
{
  if (computeFIntPtr == 0)
    RuntimeException::selfThrow("computeFInt() is not linked to a plugin function");

  int size = q->size();
  this->computeFIntPtr(&size, &time, &(*q)(0), &(*velocity)(0), &(*this->fInt)(0));
}

void LagrangianDS::computeFExt(double time)
{
  IN("LagrangianDS::computeFExt(double time)\n");
  if (computeFExtPtr == 0)
    RuntimeException::selfThrow("computeFExt() is not linked to a plugin function");

  int size = q->size();

  this->computeFExtPtr(&size, &time, &(*fExt)(0));
  OUT("LagrangianDS::computeFExt(double time)\n");

}

void LagrangianDS::computeQNLInertia()
{
  if (computeQNLInertiaPtr == 0)
    RuntimeException::selfThrow("computeQ() is not linked to a plugin function");

  int size = q->size();
  this->computeQNLInertiaPtr(&size, &(*q)(0), &(*velocity)(0), &(*QNLInertia)(0));
}

void LagrangianDS::computeQNLInertia(SimpleVector *q, SimpleVector *velocity)
{
  if (computeQNLInertiaPtr == 0)
    RuntimeException::selfThrow("computeQ() is not linked to a plugin function");

  int size = q->size();
  this->computeQNLInertiaPtr(&size, &(*q)(0), &(*velocity)(0), &(*QNLInertia)(0));
}

void LagrangianDS::computeJacobianQFInt(double time)
{
  if (computeJacobianQFIntPtr == 0)
    RuntimeException::selfThrow("computeJacobianQFInt() is not linked to a plugin function");

  int size = q->size();
  this->computeJacobianQFIntPtr(&size, &time, &(*q)(0), &(*velocity)(0), &(*jacobianQFInt)(0, 0));
}

void LagrangianDS::computeJacobianQFInt(double time, SimpleVector *q, SimpleVector *velocity)
{
  if (computeJacobianQFIntPtr == 0)
    RuntimeException::selfThrow("computeJacobianQFInt() is not linked to a plugin function");

  int size = q->size();
  this->computeJacobianQFIntPtr(&size, &time, &(*q)(0), &(*velocity)(0), &(*jacobianQFInt)(0, 0));
}

void LagrangianDS::computeJacobianVelocityFInt(double time)
{
  if (computeJacobianVelocityFIntPtr == 0)
    RuntimeException::selfThrow("computeJacobianVelocityFInt() is not linked to a plugin function");

  int size = this->q->size();
  this->computeJacobianVelocityFIntPtr(&size, &time, &(*this->q)(0), &(*this->velocity)(0), &(*jacobianVelocityFInt)(0, 0));
}
void LagrangianDS::computeJacobianVelocityFInt(double time, SimpleVector *q, SimpleVector *velocity)
{
  if (computeJacobianVelocityFIntPtr == 0)
    RuntimeException::selfThrow("computeJacobianVelocityFInt() is not linked to a plugin function");

  int size = q->size();
  this->computeJacobianVelocityFIntPtr(&size, &time, &(*q)(0), &(*velocity)(0), &(*jacobianVelocityFInt)(0, 0));
}


void LagrangianDS::computeJacobianQQNLInertia(double time)
{
  if (computeJacobianQQNLInertiaPtr == 0)
    RuntimeException::selfThrow("computeJacobianQQNLInertia() is not linked to a plugin function");

  int size = this->q->size();
  this->computeJacobianQQNLInertiaPtr(&size, &(*this->q)(0), &(*this->velocity)(0), &(*this->jacobianQQNLInertia)(0, 0));
}

void LagrangianDS::computeJacobianQQNLInertia(double time, SimpleVector *q, SimpleVector *velocity)
{
  if (computeJacobianQQNLInertiaPtr == 0)
    RuntimeException::selfThrow("computeJacobianQQNLInertia() is not linked to a plugin function");

  int size = q->size();
  this->computeJacobianQQNLInertiaPtr(&size, &(*q)(0), &(*velocity)(0), &(*jacobianQQNLInertia)(0, 0));
}

void LagrangianDS::computeJacobianVelocityQNLInertia(double time)
{
  if (computeJacobianVelocityQNLInertiaPtr == 0)
    RuntimeException::selfThrow("computeJacobianVelocityQNLInertia() is not linked to a plugin function");

  int size = this->q->size();
  this->computeJacobianVelocityQNLInertiaPtr(&size, &(*this->q)(0), &(*this->velocity)(0), &(*this->jacobianVelocityQNLInertia)(0, 0));
}
void LagrangianDS::computeJacobianVelocityQNLInertia(double time, SimpleVector *q, SimpleVector *velocity)
{
  if (computeJacobianVelocityQNLInertiaPtr == 0)
    RuntimeException::selfThrow("computeJacobianVelocityQNLInertia() is not linked to a plugin function");

  int size = q->size();
  this->computeJacobianVelocityQNLInertiaPtr(&size, &(*q)(0), &(*velocity)(0), &(*jacobianVelocityQNLInertia)(0, 0));
}

void LagrangianDS::setComputeMassFunction(const string& pluginPath, const string& functionName)
{
  IN("LagrangianDS::setComputeMassFunction\n");
  this->computeMassPtr = 0;
  cShared.setFunction(&computeMassPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->massFunctionName = plugin + ":" + functionName;

  OUT("LagrangianDS::setComputeMassFunction\n");

}

void LagrangianDS::setComputeFIntFunction(const string& pluginPath, const string& functionName)
{
  IN("LagrangianDS::setComputeFIntFunction\n");
  this->computeFIntPtr = 0;
  cShared.setFunction(&computeFIntPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->fIntFunctionName = plugin + ":" + functionName;

  OUT("LagrangianDS::setComputeFIntFunction\n");
}

void LagrangianDS::setComputeFExtFunction(const string& pluginPath, const string& functionName)
{
  IN("LagrangianDS::setComputeFExtFunction\n");
  this->computeFExtPtr = 0;
  cShared.setFunction(&computeFExtPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->fExtFunctionName = plugin + ":" + functionName;

  OUT("LagrangianDS::setComputeFExtFunction\n");
}

void LagrangianDS::setComputeQNLInertiaFunction(const string& pluginPath, const string& functionName)
{
  IN("LagrangianDS::setComputeQNLInertiaFunction\n");
  this->computeQNLInertiaPtr = 0;
  cShared.setFunction(&computeQNLInertiaPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->QNLInertiaFunctionName = plugin + ":" + functionName;

  OUT("LagrangianDS::setComputeQNLInertiaFunction\n");
}

void LagrangianDS::setComputeJacobianQFIntFunction(const string& pluginPath, const string& functionName)
{
  IN("LagrangianDS::setComputeJacobianQFIntFunction\n");
  this->computeJacobianQFIntPtr = 0;
  cShared.setFunction(&computeJacobianQFIntPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->jacobianQFIntFunctionName = plugin + ":" + functionName;

  OUT("LagrangianDS::setComputeJacobianQFIntFunction\n");
}

void LagrangianDS::setComputeJacobianVelocityFIntFunction(const string& pluginPath, const string& functionName)
{
  IN("LagrangianDS::setComputeJacobianVelocityFIntFunction\n");
  this->computeJacobianVelocityFIntPtr = 0;
  cShared.setFunction(&computeJacobianVelocityFIntPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->jacobianVelocityFIntFunctionName = plugin + ":" + functionName;

  OUT("LagrangianDS::setComputeJacobianVelocityFIntFunction\n");
}

void LagrangianDS::setComputeJacobianQQNLInertiaFunction(const string& pluginPath, const string& functionName)
{
  IN("LagrangianDS::setComputeJacobianQQNLInertiaFunction\n");
  this->computeJacobianQQNLInertiaPtr = 0;
  cShared.setFunction(&computeJacobianQQNLInertiaPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->jacobianQQNLInertiaFunctionName = plugin + ":" + functionName;

  OUT("LagrangianDS::setComputeJacobianQQNLInertiaFunction\n");
}

void LagrangianDS::setComputeJacobianVelocityQNLInertiaFunction(const string& pluginPath, const string& functionName)
{
  IN("LagrangianDS::setComputeJacobianVelocityQNLInertiaFunction\n");
  this->computeJacobianVelocityQNLInertiaPtr = 0;
  cShared.setFunction(&computeJacobianVelocityQNLInertiaPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->jacobianVelocityQNLInertiaFunctionName = plugin + ":" + functionName;

  OUT("LagrangianDS::setComputeJacobianVelocityQNLInertiaFunction\n");
}


void LagrangianDS::saveDSToXML()
{
  IN("LagrangianDS::saveDSToXML\n");

  //--- general DS data---
  saveDSDataToXML();
  // --- other data ---
  if (this->dsxml != 0)
  {
    LagrangianDSXML* lgptr = static_cast <LagrangianDSXML*>(this->dsxml);
    lgptr->setNdof(this->ndof);
    lgptr->setMPlugin(this->massFunctionName);
    lgptr->setQ(this->q);
    lgptr->setQ0(this->q0);
    lgptr->setQMemory(&(this->qMemory));
    lgptr->setVelocity(this->velocity);
    lgptr->setVelocity0(this->velocity0);
    lgptr->setVelocityMemory(&(this->velocityMemory));

    // FExt
    if (lgptr->hasFext())
    {
      if (!lgptr->isFextPlugin())
      {
        lgptr->setFextVector(this->fExt);
      }
    }
    else
    {
      lgptr->setFextPlugin(this->fExtFunctionName);
    }

    // FInt
    if (lgptr->hasFint())
    {
      if (!lgptr->isFintPlugin())
      {
        if (this->fInt->size() > 0)
          lgptr->setFintVector(this->fInt);
        else cout << "Warning : Fint can't be saved, the Fint vector is not defined." << endl;
      }
    }
    else
    {
      lgptr->setFintPlugin(this->fIntFunctionName);
    }

    // JacobianQFInt
    if (lgptr->hasJacobianQFint())
    {
      if (!lgptr->isJacobianQFintPlugin())
      {
        lgptr->setJacobianQFintMatrix(this->jacobianQFInt);
      }
    }
    else
    {
      lgptr->setJacobianQFintPlugin(this->jacobianQFIntFunctionName);
    }

    // JacobianVelocityFInt
    if (lgptr->hasJacobianVelocityFint())
    {
      if (!lgptr->isJacobianVelocityFintPlugin())
      {
        lgptr->setJacobianVelocityFintMatrix(this->jacobianVelocityFInt);
      }
    }
    else
    {
      lgptr->setJacobianVelocityFintPlugin(this->jacobianVelocityFIntFunctionName);
    }

    // JacobianQQNLInertia
    if (lgptr->hasJacobianQQNLInertia())
    {
      if (!lgptr->isJacobianQQNLInertiaPlugin())
      {
        lgptr->setJacobianQQNLInertiaMatrix(this->jacobianQQNLInertia);
      }
    }
    else
    {
      lgptr->setJacobianQQNLInertiaPlugin(this->jacobianQQNLInertiaFunctionName);
    }

    // JacobianVelocityQNLInertiaFunction
    if (lgptr->hasJacobianVelocityQNLInertia())
    {
      if (!lgptr->isJacobianVelocityQNLInertiaPlugin())
      {
        lgptr->setJacobianVelocityQNLInertiaMatrix(this->jacobianVelocityQNLInertia);
      }
    }
    else
    {
      lgptr->setJacobianVelocityQNLInertiaPlugin(this->jacobianVelocityQNLInertiaFunctionName);
    }

    // QNLInertia
    if (lgptr->hasQNLInertia())
    {
      if (!lgptr->isQNLInertiaPlugin())
      {
        lgptr->setQNLInertiaVector(this->QNLInertia);
      }
    }
    else
    {
      lgptr->setQNLInertiaPlugin(this->QNLInertiaFunctionName);
    }
  }
  else RuntimeException::selfThrow("LagrangianDS::saveDSToXML - object DSXML does not exist");
  OUT("LagrangianDS::saveDSToXML\n");
}

void LagrangianDS::display() const
{
  IN("LagrangianDS::display\n");

  cout << "-----------------------------------------------------" << endl;
  DynamicalSystem::display();
  cout << "____ data of the LagrangianDS " << endl;
  cout << "| ndof : " << this->ndof << endl;
  cout << "| q " << endl;
  if (q != 0) this->q->display();
  else cout << "-> 0" << endl;
  cout << "| q0 " << endl;
  if (q0 != 0) this->q0->display();
  else cout << "-> 0" << endl;
  cout << "| qFree " << endl;
  if (qFree != 0) this->qFree->display();
  else cout << "-> 0" << endl;
  cout << "| velocity " << endl;
  if (velocity != 0) this->velocity->display();
  else cout << "-> 0" << endl;
  cout << "| velocity0 " << endl;
  if (velocity0 != 0) this->velocity0->display();
  else cout << "-> 0" << endl;
  cout << "| velocityFree " << endl;
  if (velocityFree != 0) this->velocityFree->display();
  else cout << "-> 0" << endl;
  cout << "| p " << endl;
  if (p != 0) this->p->display();
  else cout << "-> 0" << endl;
  cout << "-----------------------------------------------------" << endl << endl;

  OUT("LagrangianDS::display\n");
}

// --- Functions for memory handling ---
void LagrangianDS::initMemory(const int& steps)
{
  IN("LagrangianDS::initMemory\n");
  DynamicalSystem::initMemory(steps);

  qMemory = SiconosMemory::SiconosMemory(steps, this->qMemory.getSiconosMemoryXML());
  velocityMemory = SiconosMemory::SiconosMemory(steps, this->velocityMemory.getSiconosMemoryXML());

  OUT("LagrangianDS::initMemory\n");
}


void LagrangianDS::swapInMemory(void)
{
  IN("LagrangianDS::swapInMemory(void)\n");

  // This operation should be made only if necessary. See todo note.
  DynamicalSystem::swapInMemory();

  this->qMemory.swap(q);
  this->velocityMemory.swap(velocity);

  // initialization of the reaction force due to the non smooth law
  this->p->zero();

  OUT("LagrangianDS::swapInMemory(void)\n");
}

LagrangianDS* LagrangianDS::convert(DynamicalSystem* ds)
{
  cout << "LagrangianDS::convert (DynamicalSystem* ds)" << endl;
  LagrangianDS* lnlds = dynamic_cast<LagrangianDS*>(ds);
  return lnlds;
}

double LagrangianDS::dsConvergenceIndicator() const
{
  double dsCvgIndic;
  // Velocity is used to calculate the indicator.
  // Remember that free state contains the previous Newton step
  SimpleVector *diff = new SimpleVector(this->velocity->size());
  dsCvgIndic = diff->norm() / (this -> velocityFree)->norm();
  delete diff;
  cout << " DSSDSD cvg" << dsCvgIndic << endl;
  return (dsCvgIndic);
}

// -- Default constructor --
LagrangianDS::LagrangianDS(): DynamicalSystem(), ndof(0), q(0), q0(0), qFree(0), velocity(0), velocity0(0),
  velocityFree(0), p(0), mass(0), fInt(0), fExt(0), QNLInertia(0), jacobianQFInt(0),
  jacobianVelocityFInt(0), jacobianQQNLInertia(0), jacobianVelocityQNLInertia(0),
  computeMassPtr(0), computeFIntPtr(0), computeFExtPtr(0), computeQNLInertiaPtr(0), computeJacobianQFIntPtr(0),
  computeJacobianVelocityFIntPtr(0), computeJacobianQQNLInertiaPtr(0), computeJacobianVelocityQNLInertiaPtr(0)
{
  IN("LagrangianDS::LagrangianDS() - Default constructor\n");
  this->DSType = LNLDS;

  // --- plugins connected to BasicPlugin ---
  this->setComputeMassFunction("BasicPlugin.so", "computeMass");
  this->setComputeFIntFunction("BasicPlugin.so", "computeFInt");
  this->setComputeFExtFunction("BasicPlugin.so", "computeFExt");
  this->setComputeQNLInertiaFunction("BasicPlugin.so", "computeQNLInertia");
  this->setComputeJacobianQFIntFunction("BasicPlugin.so", "computeJacobianQFInt");
  this->setComputeJacobianVelocityFIntFunction("BasicPlugin.so", "computeJacobianVelocityFInt");
  this->setComputeJacobianQQNLInertiaFunction("BasicPlugin.so", "computeJacobianQQNLInertia");
  this->setComputeJacobianVelocityQNLInertiaFunction("BasicPlugin.so", "computeJacobianVelocityQNLInertia");
  OUT("LagrangianDS::LagrangianDS() - Default constructor\n");
}
