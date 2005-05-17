#include "LagrangianDS.h"
#include "check.h"
#include "LinearBC.h"
#include "NLinearBC.h"
#include "PeriodicBC.h"
#include "LinearDSIO.h"
#include "LagrangianDSIO.h"
#include "LagrangianLinearDSIO.h"


// --- Constructor from an xml file ---
LagrangianDS::LagrangianDS(DSXML * dsXML): DynamicalSystem(), ndof(0), q(NULL), q0(NULL), qFree(NULL), velocity(NULL), velocity0(NULL),
  velocityFree(NULL), p(NULL), mass(NULL), fInt(NULL), fExt(NULL), QNLInertia(NULL), jacobianQFInt(NULL),
  jacobianVelocityFInt(NULL), jacobianQQNLInertia(NULL), jacobianVelocityQNLInertia(NULL),
  computeMassPtr(NULL), computeFIntPtr(NULL), computeFExtPtr(NULL), computeQNLInertiaPtr(NULL), computeJacobianQFIntPtr(NULL),
  computeJacobianVelocityFIntPtr(NULL), computeJacobianQQNLInertiaPtr(NULL), computeJacobianVelocityQNLInertiaPtr(NULL)
{
  IN("LagrangianDS::LagrangianDS() - Xml constructor\n");
  if (dsXML != NULL)
  {
    // --- DS BASE-CLASS MEMBERS ---
    // --- Settings and xml load ---
    DSType = LNLDS;
    dsxml = dsXML;
    number = dsxml->getNumber();
    if (dsxml->hasId() == true) id = dsxml->getId();
    // -- Memory allocation for vector and matrix members --
    x = new CompositeVector();
    x0 = new CompositeVector();
    xDot = new CompositeVector();
    xFree = new CompositeVector();
    //
    r = new SimpleVector(n);
    // jacobianX = new SiconosMatrix(n,n);
    if (dsxml->hasStepsInMemory() == true) stepsInMemory = dsxml->getStepsInMemory();

    // -- plugins --
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
    // -- Boundary conditions --
    fillBoundaryConditionsFromXml();
    // -- DS input-output --
    fillDsioFromXml();

    // --- LAGRANGIAN INHERITED CLASS MEMBERS ---
    // --- Settings and xml load ---
    // -- Lagrangian  xml object --
    LagrangianDSXML* lgptr = static_cast <LagrangianDSXML*>(dsxml);

    // -- Size of the system and number of degrees of freedom --
    ndof = lgptr->getNdof();
    n = 2 * ndof;

    // -- Vector and matrix members memory allocation --
    q = new SimpleVector(ndof);
    q0 = new SimpleVector(ndof);
    qFree = new SimpleVector(ndof);
    velocity = new SimpleVector(ndof);
    velocity0 = new SimpleVector(ndof);
    velocityFree = new SimpleVector(ndof);
    p = new SimpleVector(ndof);
    mass = new SiconosMatrix(ndof, ndof);
    fInt = new SimpleVector(ndof);
    fExt = new SimpleVector(ndof);
    QNLInertia = new SimpleVector(ndof);
    jacobianQFInt = new SiconosMatrix(ndof, ndof);
    jacobianVelocityFInt = new SiconosMatrix(ndof, ndof);
    jacobianQQNLInertia = new SiconosMatrix(ndof, ndof);
    jacobianVelocityQNLInertia = new SiconosMatrix(ndof, ndof);

    // -- Plugins --
    // Mass
    if ((static_cast <LagrangianDSXML*>(dsxml))->isMPlugin())
    {
      plugin = lgptr->getMPlugin();
      setComputeMassFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    // \warning : VA:  It is a very good idea to take the constant Mass Matrix, but for the moment a constant
    //  Mass Matrix is only read by a LagrangianLinearTIDS
    else *mass = lgptr->getMMatrix();

    // q0, q and qMemory
    *q0 = lgptr->getQ0();
    if (lgptr->hasQ()) *q = lgptr->getQ();
    else *q = *q0;
    if (lgptr->hasQMemory()) qMemory = SiconosMemory::SiconosMemory(lgptr->getQMemoryXML());

    // velocity0, velocity and velocityMemory
    *velocity0 = lgptr->getVelocity0();
    if (lgptr->hasVelocity()) *velocity = lgptr->getVelocity();
    else *velocity = *velocity0;
    if (lgptr->hasVelocityMemory()) velocityMemory = SiconosMemory::SiconosMemory(lgptr->getVelocityMemoryXML());

    // fill in x, x0 and xFree with q and velocity
    static_cast<CompositeVector*>(x)->add(q);
    static_cast<CompositeVector*>(x)->add(velocity);
    static_cast<CompositeVector*>(x0)->add(q0);
    static_cast<CompositeVector*>(x0)->add(velocity0);
    static_cast<CompositeVector*>(xFree)->add(qFree);
    static_cast<CompositeVector*>(xFree)->add(velocityFree);

    // FExt
    if (lgptr->isFextPlugin())
    {
      plugin = lgptr->getFextPlugin();
      setComputeFExtFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else *fExt = lgptr->getFextVector();

    // -- FInt and its jacobian --
    if (lgptr->hasFint())
    {
      if (lgptr->isFintPlugin())
      {
        plugin = lgptr->getFintPlugin();
        setComputeFIntFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
      }
      else *fInt = lgptr->getFintVector();
    }
    // Jacobian Q FInt
    if (lgptr->isJacobianQFintPlugin())
    {
      plugin = lgptr->getJacobianQFintPlugin();
      setComputeJacobianQFIntFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else *jacobianQFInt = lgptr->getJacobianQFintMatrix();
    // Jacobian Velocity FInt
    if (lgptr->isJacobianVelocityFintPlugin())
    {
      plugin = lgptr->getJacobianVelocityFintPlugin();
      setComputeJacobianVelocityFIntFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else *jacobianVelocityFInt = lgptr->getJacobianVelocityFintMatrix();

    // -- QNL Inertia and its jacobian --
    if (lgptr->isQNLInertiaPlugin())
    {
      plugin = lgptr->getQNLInertiaPlugin();
      setComputeQNLInertiaFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else *QNLInertia = lgptr->getQNLInertiaVector();
    // Jacobian Q QNLInertia
    if (lgptr->isJacobianQQNLInertiaPlugin())
    {
      plugin = lgptr->getJacobianQQNLInertiaPlugin();
      setComputeJacobianQQNLInertiaFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else *jacobianQQNLInertia = lgptr->getJacobianQQNLInertiaMatrix();
    // Jacobian Velocity QNLInertia
    if (lgptr->isJacobianVelocityQNLInertiaPlugin())
    {
      plugin = lgptr->getJacobianVelocityQNLInertiaPlugin();
      setComputeJacobianVelocityQNLInertiaFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else *jacobianVelocityQNLInertia = lgptr->getJacobianVelocityQNLInertiaMatrix();

  }
  else RuntimeException::selfThrow("LagrangianDS::LagrangianDS - DSXML paramater must not be NULL");
  OUT("LagrangianDS::LagrangianDS() - Xml constructor\n");
}

// From a minimum set of data
LagrangianDS::LagrangianDS(int newNumber, int newNdof,
                           SiconosVector* newQ0, SiconosVector* newVelocity0, string massName,
                           string fIntName, string fExtName, string jacobianQFIntName, string jacobianVelocityFIntName,
                           string jacobianQQNLInertiaName, string jacobianVelocityQNLInertiaName, string QNLInertiaName):
  DynamicalSystem(), ndof(newNdof), q(NULL), q0(NULL), qFree(NULL), velocity(NULL), velocity0(NULL), velocityFree(NULL), p(NULL), mass(NULL),
  fInt(NULL), fExt(NULL), QNLInertia(NULL), jacobianQFInt(NULL), jacobianVelocityFInt(NULL), jacobianQQNLInertia(NULL), jacobianVelocityQNLInertia(NULL),
  computeMassPtr(NULL), computeFIntPtr(NULL), computeFExtPtr(NULL), computeQNLInertiaPtr(NULL), computeJacobianQFIntPtr(NULL),
  computeJacobianVelocityFIntPtr(NULL), computeJacobianQQNLInertiaPtr(NULL), computeJacobianVelocityQNLInertiaPtr(NULL)

{
  IN("LagrangianDS::LagrangianDS - From a minimum set of data\n");
  // --- DS BASE-CLASS MEMBERS ---
  // --- Settings and xml load ---
  DSType = LNLDS;
  number = newNumber;
  n = 2 * ndof;

  // -- Memory allocation for vector and matrix members --
  x = new CompositeVector();
  x0 = new CompositeVector();
  xDot = new CompositeVector();
  xFree = new CompositeVector();
  //
  // r = new SimpleVector(n);
  // jacobianX = new SiconosMatrix(n,n);
  // -- plugins --
  string plugin;
  // VectorField
  // JacobianX

  // --- LAGRANGIAN INHERITED CLASS MEMBERS ---
  // -- Memory allocation for vector and matrix members --
  mass = new SiconosMatrix(ndof, ndof);
  q = new SimpleVector(ndof);
  q0 = new SimpleVector(ndof);
  qFree = new SimpleVector(ndof);
  velocity = new SimpleVector(ndof);
  velocity0 = new SimpleVector(ndof);
  velocityFree = new SimpleVector(ndof);
  p = new SimpleVector(ndof);
  fInt = new SimpleVector(ndof);
  fExt = new SimpleVector(ndof);
  QNLInertia = new SimpleVector(ndof);
  jacobianQFInt = new SiconosMatrix(ndof, ndof);
  jacobianVelocityFInt = new SiconosMatrix(ndof, ndof);
  jacobianQQNLInertia = new SiconosMatrix(ndof, ndof);
  jacobianVelocityQNLInertia = new SiconosMatrix(ndof, ndof);

  // --- initial state filling ---
  *q0 = *newQ0;
  *q = *newQ0;
  *velocity0 = *newVelocity0;
  *velocity = *newVelocity0;

  // --- x, xDot and xFree update ---

  static_cast<CompositeVector*>(x)->add(q);
  static_cast<CompositeVector*>(x)->add(velocity);
  static_cast<CompositeVector*>(x0)->add(q0);
  static_cast<CompositeVector*>(x0)->add(velocity0);
  static_cast<CompositeVector*>(xFree)->add(qFree);
  static_cast<CompositeVector*>(xFree)->add(velocityFree);

  //   --- plugins ---
  setComputeMassFunction(cShared.getPluginName(massName), cShared.getPluginFunctionName(massName));
  setComputeFIntFunction(cShared.getPluginName(fIntName), cShared.getPluginFunctionName(fIntName));
  setComputeFExtFunction(cShared.getPluginName(fExtName), cShared.getPluginFunctionName(fExtName));
  setComputeJacobianQFIntFunction(cShared.getPluginName(jacobianQFIntName), cShared.getPluginFunctionName(jacobianQFIntName));
  setComputeJacobianVelocityFIntFunction(cShared.getPluginName(jacobianVelocityFIntName), cShared.getPluginFunctionName(jacobianQFIntName));
  setComputeQNLInertiaFunction(cShared.getPluginName(QNLInertiaName), cShared.getPluginFunctionName(QNLInertiaName));
  setComputeJacobianQQNLInertiaFunction(cShared.getPluginName(jacobianQQNLInertiaName), cShared.getPluginFunctionName(jacobianQQNLInertiaName));
  setComputeJacobianVelocityQNLInertiaFunction(cShared.getPluginName(jacobianVelocityQNLInertiaName), cShared.getPluginFunctionName(jacobianVelocityQNLInertiaName));
  OUT("LagrangianDS::LagrangianDS - From a minimum set of data\n");
}

LagrangianDS::~LagrangianDS()
{
  IN("LagrangianDS::~LagrangianDS()\n");
  delete q;
  q = NULL;
  delete q0;
  q0 = NULL;
  delete qFree;
  qFree = NULL;
  delete velocity ;
  velocity = NULL;
  delete velocity0 ;
  velocity0 = NULL;
  delete velocityFree ;
  velocityFree = NULL;
  delete p ;
  p = NULL;
  delete mass;
  mass = NULL;
  delete fInt ;
  fInt = NULL;
  delete fExt ;
  fExt = NULL;
  delete QNLInertia ;
  QNLInertia = NULL;
  delete jacobianQFInt  ;
  jacobianQFInt = NULL;
  delete  jacobianVelocityFInt ;
  jacobianVelocityFInt = NULL;
  delete jacobianQQNLInertia ;
  jacobianQQNLInertia = NULL;
  delete jacobianVelocityQNLInertia ;
  jacobianVelocityQNLInertia = NULL;
  OUT("LagrangianDS::~LagrangianDS()\n");
}


// --- Plugins related functions ---
void LagrangianDS::computeMass(double time)
{
  if (computeMassPtr == NULL)
    RuntimeException::selfThrow("computeMass() is not linked to a plugin function");
  int size = q->size();
  computeMassPtr(&size, &time, &(*q)(0), &(*mass)(0, 0));
}

void LagrangianDS::computeMass(double time, SimpleVector *q)
{
  if (computeMassPtr == NULL)
    RuntimeException::selfThrow("computeMass() is not linked to a plugin function");

  int size = q->size();
  computeMassPtr(&size, &time, &(*q)(0), &(*mass)(0, 0));
}

void LagrangianDS::computeFInt(double time)
{
  if (computeFIntPtr == NULL)
    RuntimeException::selfThrow("computeFInt() is not linked to a plugin function");

  int size = q->size();
  computeFIntPtr(&size, &time, &(*q)(0), &(*velocity)(0), &(*fInt)(0));
}
void LagrangianDS::computeFInt(double time, SimpleVector *q, SimpleVector *velocity)
{
  if (computeFIntPtr == NULL)
    RuntimeException::selfThrow("computeFInt() is not linked to a plugin function");

  int size = q->size();
  computeFIntPtr(&size, &time, &(*q)(0), &(*velocity)(0), &(*fInt)(0));
}

void LagrangianDS::computeFExt(double time)
{
  IN("LagrangianDS::computeFExt(double time)\n");
  if (computeFExtPtr == NULL)
    RuntimeException::selfThrow("computeFExt() is not linked to a plugin function");

  int size = q->size();

  computeFExtPtr(&size, &time, &(*fExt)(0));
  OUT("LagrangianDS::computeFExt(double time)\n");

}

void LagrangianDS::computeQNLInertia()
{
  if (computeQNLInertiaPtr == NULL)
    RuntimeException::selfThrow("computeQ() is not linked to a plugin function");

  int size = q->size();
  computeQNLInertiaPtr(&size, &(*q)(0), &(*velocity)(0), &(*QNLInertia)(0));
}

void LagrangianDS::computeQNLInertia(SimpleVector *q, SimpleVector *velocity)
{
  if (computeQNLInertiaPtr == NULL)
    RuntimeException::selfThrow("computeQ() is not linked to a plugin function");

  int size = q->size();
  computeQNLInertiaPtr(&size, &(*q)(0), &(*velocity)(0), &(*QNLInertia)(0));
}

void LagrangianDS::computeJacobianQFInt(double time)
{
  if (computeJacobianQFIntPtr == NULL)
    RuntimeException::selfThrow("computeJacobianQFInt() is not linked to a plugin function");

  int size = q->size();
  computeJacobianQFIntPtr(&size, &time, &(*q)(0), &(*velocity)(0), &(*jacobianQFInt)(0, 0));
}

void LagrangianDS::computeJacobianQFInt(double time, SimpleVector *q, SimpleVector *velocity)
{
  if (computeJacobianQFIntPtr == NULL)
    RuntimeException::selfThrow("computeJacobianQFInt() is not linked to a plugin function");

  int size = q->size();
  computeJacobianQFIntPtr(&size, &time, &(*q)(0), &(*velocity)(0), &(*jacobianQFInt)(0, 0));
}

void LagrangianDS::computeJacobianVelocityFInt(double time)
{
  if (computeJacobianVelocityFIntPtr == NULL)
    RuntimeException::selfThrow("computeJacobianVelocityFInt() is not linked to a plugin function");

  int size = q->size();
  computeJacobianVelocityFIntPtr(&size, &time, &(*q)(0), &(*velocity)(0), &(*jacobianVelocityFInt)(0, 0));
}
void LagrangianDS::computeJacobianVelocityFInt(double time, SimpleVector *q, SimpleVector *velocity)
{
  if (computeJacobianVelocityFIntPtr == NULL)
    RuntimeException::selfThrow("computeJacobianVelocityFInt() is not linked to a plugin function");

  int size = q->size();
  computeJacobianVelocityFIntPtr(&size, &time, &(*q)(0), &(*velocity)(0), &(*jacobianVelocityFInt)(0, 0));
}


void LagrangianDS::computeJacobianQQNLInertia(double time)
{
  if (computeJacobianQQNLInertiaPtr == NULL)
    RuntimeException::selfThrow("computeJacobianQQNLInertia() is not linked to a plugin function");

  int size = q->size();
  computeJacobianQQNLInertiaPtr(&size, &(*q)(0), &(*velocity)(0), &(*jacobianQQNLInertia)(0, 0));
}

void LagrangianDS::computeJacobianQQNLInertia(double time, SimpleVector *q, SimpleVector *velocity)
{
  if (computeJacobianQQNLInertiaPtr == NULL)
    RuntimeException::selfThrow("computeJacobianQQNLInertia() is not linked to a plugin function");

  int size = q->size();
  computeJacobianQQNLInertiaPtr(&size, &(*q)(0), &(*velocity)(0), &(*jacobianQQNLInertia)(0, 0));
}

void LagrangianDS::computeJacobianVelocityQNLInertia(double time)
{
  if (computeJacobianVelocityQNLInertiaPtr == NULL)
    RuntimeException::selfThrow("computeJacobianVelocityQNLInertia() is not linked to a plugin function");

  int size = q->size();
  computeJacobianVelocityQNLInertiaPtr(&size, &(*q)(0), &(*velocity)(0), &(*jacobianVelocityQNLInertia)(0, 0));
}
void LagrangianDS::computeJacobianVelocityQNLInertia(double time, SimpleVector *q, SimpleVector *velocity)
{
  if (computeJacobianVelocityQNLInertiaPtr == NULL)
    RuntimeException::selfThrow("computeJacobianVelocityQNLInertia() is not linked to a plugin function");

  int size = q->size();
  computeJacobianVelocityQNLInertiaPtr(&size, &(*q)(0), &(*velocity)(0), &(*jacobianVelocityQNLInertia)(0, 0));
}

void LagrangianDS::setComputeMassFunction(const string& pluginPath, const string& functionName)
{
  IN("LagrangianDS::setComputeMassFunction\n");
  computeMassPtr = NULL;
  cShared.setFunction(&computeMassPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  massFunctionName = plugin + ":" + functionName;

  OUT("LagrangianDS::setComputeMassFunction\n");

}

void LagrangianDS::setComputeFIntFunction(const string& pluginPath, const string& functionName)
{
  IN("LagrangianDS::setComputeFIntFunction\n");
  computeFIntPtr = NULL;
  cShared.setFunction(&computeFIntPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  fIntFunctionName = plugin + ":" + functionName;

  OUT("LagrangianDS::setComputeFIntFunction\n");
}

void LagrangianDS::setComputeFExtFunction(const string& pluginPath, const string& functionName)
{
  IN("LagrangianDS::setComputeFExtFunction\n");
  computeFExtPtr = NULL;
  cShared.setFunction(&computeFExtPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  fExtFunctionName = plugin + ":" + functionName;

  OUT("LagrangianDS::setComputeFExtFunction\n");
}

void LagrangianDS::setComputeQNLInertiaFunction(const string& pluginPath, const string& functionName)
{
  IN("LagrangianDS::setComputeQNLInertiaFunction\n");
  computeQNLInertiaPtr = NULL;
  cShared.setFunction(&computeQNLInertiaPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  QNLInertiaFunctionName = plugin + ":" + functionName;

  OUT("LagrangianDS::setComputeQNLInertiaFunction\n");
}

void LagrangianDS::setComputeJacobianQFIntFunction(const string& pluginPath, const string& functionName)
{
  IN("LagrangianDS::setComputeJacobianQFIntFunction\n");
  computeJacobianQFIntPtr = NULL;
  cShared.setFunction(&computeJacobianQFIntPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  jacobianQFIntFunctionName = plugin + ":" + functionName;

  OUT("LagrangianDS::setComputeJacobianQFIntFunction\n");
}

void LagrangianDS::setComputeJacobianVelocityFIntFunction(const string& pluginPath, const string& functionName)
{
  IN("LagrangianDS::setComputeJacobianVelocityFIntFunction\n");
  computeJacobianVelocityFIntPtr = NULL;
  cShared.setFunction(&computeJacobianVelocityFIntPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  jacobianVelocityFIntFunctionName = plugin + ":" + functionName;

  OUT("LagrangianDS::setComputeJacobianVelocityFIntFunction\n");
}

void LagrangianDS::setComputeJacobianQQNLInertiaFunction(const string& pluginPath, const string& functionName)
{
  IN("LagrangianDS::setComputeJacobianQQNLInertiaFunction\n");
  computeJacobianQQNLInertiaPtr = NULL;
  cShared.setFunction(&computeJacobianQQNLInertiaPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  jacobianQQNLInertiaFunctionName = plugin + ":" + functionName;

  OUT("LagrangianDS::setComputeJacobianQQNLInertiaFunction\n");
}

void LagrangianDS::setComputeJacobianVelocityQNLInertiaFunction(const string& pluginPath, const string& functionName)
{
  IN("LagrangianDS::setComputeJacobianVelocityQNLInertiaFunction\n");
  computeJacobianVelocityQNLInertiaPtr = NULL;
  cShared.setFunction(&computeJacobianVelocityQNLInertiaPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  jacobianVelocityQNLInertiaFunctionName = plugin + ":" + functionName;

  OUT("LagrangianDS::setComputeJacobianVelocityQNLInertiaFunction\n");
}


void LagrangianDS::saveDSToXML()
{
  IN("LagrangianDS::saveDSToXML\n");

  //--- general DS data---
  saveDSDataToXML();
  // --- other data ---
  if (dsxml != NULL)
  {
    LagrangianDSXML* lgptr = static_cast <LagrangianDSXML*>(dsxml);
    lgptr->setNdof(ndof);
    lgptr->setMPlugin(massFunctionName);
    lgptr->setQ(q);
    lgptr->setQ0(q0);
    lgptr->setQMemory(&(qMemory));
    lgptr->setVelocity(velocity);
    lgptr->setVelocity0(velocity0);
    lgptr->setVelocityMemory(&(velocityMemory));

    // FExt
    if (lgptr->hasFext())
    {
      if (!lgptr->isFextPlugin())
      {
        lgptr->setFextVector(fExt);
      }
    }
    else
    {
      lgptr->setFextPlugin(fExtFunctionName);
    }

    // FInt
    if (lgptr->hasFint())
    {
      if (!lgptr->isFintPlugin())
      {
        if (fInt->size() > 0)
          lgptr->setFintVector(fInt);
        else cout << "Warning : Fint can't be saved, the Fint vector is not defined." << endl;
      }
    }
    else
    {
      lgptr->setFintPlugin(fIntFunctionName);
    }

    // JacobianQFInt
    if (lgptr->hasJacobianQFint())
    {
      if (!lgptr->isJacobianQFintPlugin())
      {
        lgptr->setJacobianQFintMatrix(jacobianQFInt);
      }
    }
    else
    {
      lgptr->setJacobianQFintPlugin(jacobianQFIntFunctionName);
    }

    // JacobianVelocityFInt
    if (lgptr->hasJacobianVelocityFint())
    {
      if (!lgptr->isJacobianVelocityFintPlugin())
      {
        lgptr->setJacobianVelocityFintMatrix(jacobianVelocityFInt);
      }
    }
    else
    {
      lgptr->setJacobianVelocityFintPlugin(jacobianVelocityFIntFunctionName);
    }

    // JacobianQQNLInertia
    if (lgptr->hasJacobianQQNLInertia())
    {
      if (!lgptr->isJacobianQQNLInertiaPlugin())
      {
        lgptr->setJacobianQQNLInertiaMatrix(jacobianQQNLInertia);
      }
    }
    else
    {
      lgptr->setJacobianQQNLInertiaPlugin(jacobianQQNLInertiaFunctionName);
    }

    // JacobianVelocityQNLInertiaFunction
    if (lgptr->hasJacobianVelocityQNLInertia())
    {
      if (!lgptr->isJacobianVelocityQNLInertiaPlugin())
      {
        lgptr->setJacobianVelocityQNLInertiaMatrix(jacobianVelocityQNLInertia);
      }
    }
    else
    {
      lgptr->setJacobianVelocityQNLInertiaPlugin(jacobianVelocityQNLInertiaFunctionName);
    }

    // QNLInertia
    if (lgptr->hasQNLInertia())
    {
      if (!lgptr->isQNLInertiaPlugin())
      {
        lgptr->setQNLInertiaVector(QNLInertia);
      }
    }
    else
    {
      lgptr->setQNLInertiaPlugin(QNLInertiaFunctionName);
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
  cout << "| ndof : " << ndof << endl;
  cout << "| q " << endl;
  if (q != NULL) q->display();
  else cout << "-> NULL" << endl;
  cout << "| q0 " << endl;
  if (q0 != NULL) q0->display();
  else cout << "-> NULL" << endl;
  cout << "| qFree " << endl;
  if (qFree != NULL) qFree->display();
  else cout << "-> NULL" << endl;
  cout << "| velocity " << endl;
  if (velocity != NULL) velocity->display();
  else cout << "-> NULL" << endl;
  cout << "| velocity0 " << endl;
  if (velocity0 != NULL) velocity0->display();
  else cout << "-> NULL" << endl;
  cout << "| velocityFree " << endl;
  if (velocityFree != NULL) velocityFree->display();
  else cout << "-> NULL" << endl;
  cout << "| p " << endl;
  if (p != NULL) p->display();
  else cout << "-> NULL" << endl;
  cout << "-----------------------------------------------------" << endl << endl;

  OUT("LagrangianDS::display\n");
}

// --- Functions for memory handling ---
void LagrangianDS::initMemory(const int& steps)
{
  IN("LagrangianDS::initMemory\n");
  DynamicalSystem::initMemory(steps);

  qMemory = SiconosMemory::SiconosMemory(steps, qMemory.getSiconosMemoryXML());
  velocityMemory = SiconosMemory::SiconosMemory(steps, velocityMemory.getSiconosMemoryXML());

  OUT("LagrangianDS::initMemory\n");
}


void LagrangianDS::swapInMemory(void)
{
  IN("LagrangianDS::swapInMemory(void)\n");

  // This operation should be made only if necessary. See todo note.
  DynamicalSystem::swapInMemory();

  qMemory.swap(q);
  velocityMemory.swap(velocity);

  // initialization of the reaction force due to the non smooth law
  p->zero();

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
  SimpleVector *diff = new SimpleVector(velocity->size());
  // Compute difference between present and previous Newton steps
  *diff =  *(getVelocityPtr()) - *(getVelocityFreePtr());
  dsCvgIndic = diff->norm(); // / velocityFree->norm();
  delete diff;
  return (dsCvgIndic);
}

// -- Default constructor --
LagrangianDS::LagrangianDS(): DynamicalSystem(), ndof(0), q(NULL), q0(NULL), qFree(NULL), velocity(NULL), velocity0(NULL),
  velocityFree(NULL), p(NULL), mass(NULL), fInt(NULL), fExt(NULL), QNLInertia(NULL), jacobianQFInt(NULL),
  jacobianVelocityFInt(NULL), jacobianQQNLInertia(NULL), jacobianVelocityQNLInertia(NULL),
  computeMassPtr(NULL), computeFIntPtr(NULL), computeFExtPtr(NULL), computeQNLInertiaPtr(NULL), computeJacobianQFIntPtr(NULL),
  computeJacobianVelocityFIntPtr(NULL), computeJacobianQQNLInertiaPtr(NULL), computeJacobianVelocityQNLInertiaPtr(NULL)
{
  IN("LagrangianDS::LagrangianDS() - Default constructor\n");
  DSType = LNLDS;

  // --- plugins connected to BasicPlugin ---
  setComputeMassFunction("BasicPlugin.so", "computeMass");
  setComputeFIntFunction("BasicPlugin.so", "computeFInt");
  setComputeFExtFunction("BasicPlugin.so", "computeFExt");
  setComputeQNLInertiaFunction("BasicPlugin.so", "computeQNLInertia");
  setComputeJacobianQFIntFunction("BasicPlugin.so", "computeJacobianQFInt");
  setComputeJacobianVelocityFIntFunction("BasicPlugin.so", "computeJacobianVelocityFInt");
  setComputeJacobianQQNLInertiaFunction("BasicPlugin.so", "computeJacobianQQNLInertia");
  setComputeJacobianVelocityQNLInertiaFunction("BasicPlugin.so", "computeJacobianVelocityQNLInertia");
  OUT("LagrangianDS::LagrangianDS() - Default constructor\n");
}
