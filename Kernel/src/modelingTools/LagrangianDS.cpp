#include "LagrangianDS.h"
using namespace std;

// --- Constructor from an xml file ---
LagrangianDS::LagrangianDS(DSXML * dsXML):
  DynamicalSystem(), ndof(0), q(NULL), q0(NULL), qFree(NULL), velocity(NULL), velocity0(NULL),
  velocityFree(NULL), p(NULL), mass(NULL), fInt(NULL), fExt(NULL), QNLInertia(NULL), jacobianQFInt(NULL),
  jacobianVelocityFInt(NULL), jacobianQQNLInertia(NULL), jacobianVelocityQNLInertia(NULL),
  isPAllocatedIn(true), isMassAllocatedIn(true),
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
    isXAllocatedIn = true;
    x0 = new CompositeVector();
    isX0AllocatedIn = true;
    xDot = new CompositeVector();
    isXDotAllocatedIn = true;
    xFree = new CompositeVector();
    isXFreeAllocatedIn = true;
    //
    r = new SimpleVector(n);
    isRAllocatedIn = true;
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
    isQAllocatedIn.resize(4, true);
    isQAllocatedIn[3] = false; // qMemory
    velocity = new SimpleVector(ndof);
    velocity0 = new SimpleVector(ndof);
    velocityFree = new SimpleVector(ndof);
    isVelocityAllocatedIn.resize(4, true);
    isVelocityAllocatedIn[3] = false; // velocityMemory
    p = new SimpleVector(ndof);
    mass = new SiconosMatrix(ndof, ndof);
    fInt = new SimpleVector(ndof);
    fExt = new SimpleVector(ndof);
    QNLInertia = new SimpleVector(ndof);
    areForcesAllocatedIn.resize(3, true);
    jacobianQFInt = new SiconosMatrix(ndof, ndof);
    jacobianVelocityFInt = new SiconosMatrix(ndof, ndof);
    jacobianQQNLInertia = new SiconosMatrix(ndof, ndof);
    jacobianVelocityQNLInertia = new SiconosMatrix(ndof, ndof);
    isJacobianAllocatedIn.resize(4, true);

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
    if (lgptr->hasQMemory())
    {
      qMemory = new SiconosMemory(lgptr->getQMemoryXML());
      isQAllocatedIn[3] = true;
    }

    // velocity0, velocity and velocityMemory
    *velocity0 = lgptr->getVelocity0();
    if (lgptr->hasVelocity()) *velocity = lgptr->getVelocity();
    else *velocity = *velocity0;
    if (lgptr->hasVelocityMemory())
    {
      velocityMemory = new SiconosMemory(lgptr->getVelocityMemoryXML());
      isVelocityAllocatedIn[3] = true;
    }

    // fill in x, x0 and xFree with q and velocity
    static_cast<CompositeVector*>(x)->addPtr(q);
    static_cast<CompositeVector*>(x)->addPtr(velocity);
    static_cast<CompositeVector*>(x0)->addPtr(q0);
    static_cast<CompositeVector*>(x0)->addPtr(velocity0);
    static_cast<CompositeVector*>(xFree)->addPtr(qFree);
    static_cast<CompositeVector*>(xFree)->addPtr(velocityFree);

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
    if (lgptr ->hasJacobianQFint())
    {
      if (lgptr->isJacobianQFintPlugin())
      {
        plugin = lgptr->getJacobianQFintPlugin();
        setComputeJacobianQFIntFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
      }
      else *jacobianQFInt = lgptr->getJacobianQFintMatrix();
    }
    // Jacobian Velocity FInt

    if (lgptr ->hasJacobianVelocityFint())
    {
      if (lgptr->isJacobianVelocityFintPlugin())
      {
        plugin = lgptr->getJacobianVelocityFintPlugin();
        setComputeJacobianVelocityFIntFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
      }
      else *jacobianVelocityFInt = lgptr->getJacobianVelocityFintMatrix();
    }

    // -- QNL Inertia and its jacobian --
    if (lgptr ->hasQNLInertia())
    {
      if (lgptr->isQNLInertiaPlugin())
      {
        plugin = lgptr->getQNLInertiaPlugin();
        setComputeQNLInertiaFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
      }
      else *QNLInertia = lgptr->getQNLInertiaVector();
    }

    // Jacobian Q QNLInertia
    if (lgptr -> hasJacobianQQNLInertia())
    {
      if (lgptr->isJacobianQQNLInertiaPlugin())
      {
        plugin = lgptr->getJacobianQQNLInertiaPlugin();
        setComputeJacobianQQNLInertiaFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
      }
      else *jacobianQQNLInertia = lgptr->getJacobianQQNLInertiaMatrix();
    }

    // Jacobian Velocity QNLInertia
    if (lgptr ->  hasJacobianVelocityQNLInertia())
    {
      if (lgptr->isJacobianVelocityQNLInertiaPlugin())
      {
        plugin = lgptr->getJacobianVelocityQNLInertiaPlugin();
        setComputeJacobianVelocityQNLInertiaFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
      }
      else *jacobianVelocityQNLInertia = lgptr->getJacobianVelocityQNLInertiaMatrix();
    }
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
  isPAllocatedIn(true), isMassAllocatedIn(true),
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
  isXAllocatedIn = true;
  x0 = new CompositeVector();
  isX0AllocatedIn = true;
  xDot = new CompositeVector();
  isXDotAllocatedIn = true;
  xFree = new CompositeVector();
  isXFreeAllocatedIn = true;

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
  isQAllocatedIn.resize(4, true);
  isQAllocatedIn[3] = false; // qMemory
  velocity = new SimpleVector(ndof);
  velocity0 = new SimpleVector(ndof);
  velocityFree = new SimpleVector(ndof);
  isVelocityAllocatedIn.resize(4, true);
  isVelocityAllocatedIn[3] = false; // velocityMemory
  p = new SimpleVector(ndof);
  fInt = new SimpleVector(ndof);
  fExt = new SimpleVector(ndof);
  QNLInertia = new SimpleVector(ndof);
  areForcesAllocatedIn.resize(3, true);
  jacobianQFInt = new SiconosMatrix(ndof, ndof);
  jacobianVelocityFInt = new SiconosMatrix(ndof, ndof);
  jacobianQQNLInertia = new SiconosMatrix(ndof, ndof);
  jacobianVelocityQNLInertia = new SiconosMatrix(ndof, ndof);
  isJacobianAllocatedIn.resize(4, true);

  // --- initial state filling ---
  *q0 = *newQ0;
  *q = *newQ0;
  *velocity0 = *newVelocity0;
  *velocity = *newVelocity0;

  // --- x, xDot and xFree update ---

  static_cast<CompositeVector*>(x)->addPtr(q);
  static_cast<CompositeVector*>(x)->addPtr(velocity);
  static_cast<CompositeVector*>(x0)->addPtr(q0);
  static_cast<CompositeVector*>(x0)->addPtr(velocity0);
  static_cast<CompositeVector*>(xFree)->addPtr(qFree);
  static_cast<CompositeVector*>(xFree)->addPtr(velocityFree);

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
  if (isQAllocatedIn[0])
  {
    delete q;
    q = NULL;
  }
  if (isQAllocatedIn[1])
  {
    delete q0;
    q0 = NULL;
  }
  if (isQAllocatedIn[2])
  {
    delete qFree;
    qFree = NULL;
  }
  if (isQAllocatedIn[3])
  {
    delete qMemory;
    qMemory = NULL;
  }
  if (isVelocityAllocatedIn[0])
  {
    delete velocity ;
    velocity = NULL;
  }
  if (isVelocityAllocatedIn[1])
  {
    delete velocity0 ;
    velocity0 = NULL;
  }
  if (isVelocityAllocatedIn[2])
  {
    delete velocityFree ;
    velocityFree = NULL;
  }
  if (isVelocityAllocatedIn[3])
  {
    delete velocityMemory;
    velocityMemory = NULL;
  }
  if (isPAllocatedIn)
  {
    delete p ;
    p = NULL;
  }
  if (isMassAllocatedIn)
  {
    delete mass;
    mass = NULL;
  }
  if (areForcesAllocatedIn[0])
  {
    delete fInt ;
    fInt = NULL;
  }
  if (areForcesAllocatedIn[1])
  {
    delete fExt ;
    fExt = NULL;
  }
  if (areForcesAllocatedIn[2])
  {
    delete QNLInertia ;
    QNLInertia = NULL;
  }
  if (isJacobianAllocatedIn[0])
  {
    delete jacobianQFInt  ;
    jacobianQFInt = NULL;
  }
  if (isJacobianAllocatedIn[0])
  {
    delete  jacobianVelocityFInt ;
    jacobianVelocityFInt = NULL;
  }
  if (isJacobianAllocatedIn[0])
  {
    delete jacobianQQNLInertia ;
    jacobianQQNLInertia = NULL;
  }
  if (isJacobianAllocatedIn[0])
  {
    delete jacobianVelocityQNLInertia ;
    jacobianVelocityQNLInertia = NULL;
  }
  OUT("LagrangianDS::~LagrangianDS()\n");
}

// --- GETTERS/SETTERS ---

void LagrangianDS::setQPtr(SimpleVector *newPtr)
{
  if (isQAllocatedIn[0]) delete q;
  q = newPtr;
  isQAllocatedIn[0] = false;
}

void LagrangianDS::setQ0Ptr(SimpleVector *newPtr)
{
  if (isQAllocatedIn[1]) delete q0;
  q0 = newPtr;
  isQAllocatedIn[1] = false;
}

void LagrangianDS::setQFreePtr(SimpleVector *newPtr)
{
  if (isQAllocatedIn[2]) delete qFree;
  qFree = newPtr;
  isQAllocatedIn[2] = false;
}

void LagrangianDS::setQMemoryPtr(SiconosMemory * newPtr)
{
  if (isQAllocatedIn[3]) delete qMemory;
  qMemory = newPtr;
  isQAllocatedIn[3] = false;
}

void LagrangianDS::setVelocityPtr(SimpleVector *newPtr)
{
  if (isVelocityAllocatedIn[0]) delete velocity;
  velocity = newPtr;
  isVelocityAllocatedIn[0] = false;
}

void LagrangianDS::setVelocity0Ptr(SimpleVector *newPtr)
{
  if (isVelocityAllocatedIn[1]) delete velocity0;
  velocity0 = newPtr;
  isVelocityAllocatedIn[1] = false;
}

void LagrangianDS::setVelocityFreePtr(SimpleVector *newPtr)
{
  if (isVelocityAllocatedIn[2]) delete velocityFree;
  velocityFree = newPtr;
  isVelocityAllocatedIn[2] = false;
}

void LagrangianDS::setVelocityMemoryPtr(SiconosMemory * newPtr)
{
  if (isVelocityAllocatedIn[3]) delete velocityMemory;
  velocityMemory = newPtr;
  isVelocityAllocatedIn[3] = false;
}

void LagrangianDS::setPPtr(SimpleVector *newPtr)
{
  if (isPAllocatedIn) delete p;
  p = newPtr;
  isPAllocatedIn = false;
}

void LagrangianDS::setMassPtr(SiconosMatrix *newPtr)
{
  if (isMassAllocatedIn) delete mass;
  mass = newPtr;
  isMassAllocatedIn = false;
}

void LagrangianDS::setFIntPtr(SimpleVector *newPtr)
{
  if (areForcesAllocatedIn[0]) delete fInt;
  fInt = newPtr;
  areForcesAllocatedIn[0] = false;
}

void LagrangianDS::setFExtPtr(SimpleVector *newPtr)
{
  if (areForcesAllocatedIn[1]) delete fExt;
  fExt = newPtr;
  areForcesAllocatedIn[1] = false;
}

void LagrangianDS::setQNLInertiaPtr(SimpleVector *newPtr)
{
  if (areForcesAllocatedIn[2]) delete QNLInertia;
  QNLInertia = newPtr;
  areForcesAllocatedIn[2] = false;
}

void LagrangianDS::setJacobianQFIntPtr(SiconosMatrix *newPtr)
{
  if (isJacobianAllocatedIn[0]) delete jacobianQFInt;
  jacobianQFInt = newPtr;
  isJacobianAllocatedIn[0] = false;
}

void LagrangianDS::setJacobianVelocityFIntPtr(SiconosMatrix *newPtr)
{
  if (isJacobianAllocatedIn[1]) delete jacobianVelocityFInt;
  jacobianVelocityFInt = newPtr;
  isJacobianAllocatedIn[1] = false;
}

void LagrangianDS::setJacobianQQNLInertiaPtr(SiconosMatrix *newPtr)
{
  if (isJacobianAllocatedIn[2]) delete jacobianQQNLInertia;
  jacobianQQNLInertia = newPtr;
  isJacobianAllocatedIn[2] = false;
}

void LagrangianDS::setJacobianVelocityQNLInertiaPtr(SiconosMatrix *newPtr)
{
  if (isJacobianAllocatedIn[3]) delete jacobianVelocityQNLInertia;
  jacobianVelocityQNLInertia = newPtr;
  isJacobianAllocatedIn[3] = false;
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
    lgptr->setQMemory(qMemory);
    lgptr->setVelocity(velocity);
    lgptr->setVelocity0(velocity0);
    lgptr->setVelocityMemory(velocityMemory);

    // FExt
    if (lgptr->hasFext())
    {
      if (!lgptr->isFextPlugin())
        lgptr->setFextVector(fExt);
    }
    else
      lgptr->setFextPlugin(fExtFunctionName);

    // FInt
    if (lgptr->hasFint())
    {
      if (!lgptr->isFintPlugin())
        if (fInt->size() > 0)
          lgptr->setFintVector(fInt);
        else cout << "Warning : Fint can't be saved, the Fint vector is not defined." << endl;
    }
    else
      lgptr->setFintPlugin(fIntFunctionName);

    // JacobianQFInt
    if (lgptr->hasJacobianQFint())
    {
      if (!lgptr->isJacobianQFintPlugin())
        lgptr->setJacobianQFintMatrix(jacobianQFInt);
    }
    else
      lgptr->setJacobianQFintPlugin(jacobianQFIntFunctionName);

    // JacobianVelocityFInt
    if (lgptr->hasJacobianVelocityFint())
    {
      if (!lgptr->isJacobianVelocityFintPlugin())
        lgptr->setJacobianVelocityFintMatrix(jacobianVelocityFInt);
    }
    else
      lgptr->setJacobianVelocityFintPlugin(jacobianVelocityFIntFunctionName);

    // JacobianQQNLInertia
    if (lgptr->hasJacobianQQNLInertia())
    {
      if (!lgptr->isJacobianQQNLInertiaPlugin())
        lgptr->setJacobianQQNLInertiaMatrix(jacobianQQNLInertia);
    }
    else
      lgptr->setJacobianQQNLInertiaPlugin(jacobianQQNLInertiaFunctionName);

    // JacobianVelocityQNLInertiaFunction
    if (lgptr->hasJacobianVelocityQNLInertia())
    {
      if (!lgptr->isJacobianVelocityQNLInertiaPlugin())
        lgptr->setJacobianVelocityQNLInertiaMatrix(jacobianVelocityQNLInertia);
    }
    else
      lgptr->setJacobianVelocityQNLInertiaPlugin(jacobianVelocityQNLInertiaFunctionName);

    // QNLInertia
    if (lgptr->hasQNLInertia())
    {
      if (!lgptr->isQNLInertiaPlugin())
        lgptr->setQNLInertiaVector(QNLInertia);
    }
    else
      lgptr->setQNLInertiaPlugin(QNLInertiaFunctionName);
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
void LagrangianDS::initMemory(const unsigned int& steps)
{
  IN("LagrangianDS::initMemory\n");
  DynamicalSystem::initMemory(steps);
  // warning : depends on xml loading ??? To be reviewed
  //qMemory = new SiconosMemory( steps, qMemory.getSiconosMemoryXML() );
  //velocityMemory = new SiconosMemory( steps, velocityMemory.getSiconosMemoryXML() );
  if (isQAllocatedIn[3]) delete qMemory;
  qMemory = new SiconosMemory(steps);
  isQAllocatedIn[3] = true;
  if (isVelocityAllocatedIn[3]) delete velocityMemory;
  velocityMemory = new SiconosMemory(steps);
  isVelocityAllocatedIn[3] = true;
  OUT("LagrangianDS::initMemory\n");
}


void LagrangianDS::swapInMemory(void)
{
  IN("LagrangianDS::swapInMemory(void)\n");
  DynamicalSystem::swapInMemory();
  qMemory->swap(*q);
  velocityMemory->swap(*velocity);
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

double LagrangianDS::dsConvergenceIndicator()
{
  double dsCvgIndic;
  // Velocity is used to calculate the indicator.
  SimpleVector *diff = new SimpleVector(velocity->size());
  // Compute difference between present and previous Newton steps
  SimpleVector * valRef = tmpWorkVector["LagNLDSMoreau"];
  *diff =  *(getVelocityPtr()) - *valRef;
  if (valRef->norm() != 0)
    dsCvgIndic = diff->norm() / (valRef->norm());
  else
    dsCvgIndic = diff->norm();
  delete diff;
  return (dsCvgIndic);
}

// -- Default constructor --
LagrangianDS::LagrangianDS():
  DynamicalSystem(), ndof(0), q(NULL), q0(NULL), qFree(NULL), velocity(NULL), velocity0(NULL),
  velocityFree(NULL), p(NULL), mass(NULL), fInt(NULL), fExt(NULL), QNLInertia(NULL), jacobianQFInt(NULL),
  jacobianVelocityFInt(NULL), jacobianQQNLInertia(NULL), jacobianVelocityQNLInertia(NULL),
  isPAllocatedIn(false), isMassAllocatedIn(false),
  computeMassPtr(NULL), computeFIntPtr(NULL), computeFExtPtr(NULL), computeQNLInertiaPtr(NULL), computeJacobianQFIntPtr(NULL),
  computeJacobianVelocityFIntPtr(NULL), computeJacobianQQNLInertiaPtr(NULL), computeJacobianVelocityQNLInertiaPtr(NULL)
{
  IN("LagrangianDS::LagrangianDS() - Default constructor\n");
  DSType = LNLDS;
  isQAllocatedIn.resize(4, false);
  isVelocityAllocatedIn.resize(4, false);
  areForcesAllocatedIn.resize(3, false);
  isJacobianAllocatedIn.resize(4, false);

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
