#include "LinearDS.h"
using namespace std;

// --- Constructors ---

// From xml file
LinearDS::LinearDS(DSXML * dsXML):
  DynamicalSystem(dsXML), A(NULL), f(NULL), uSize(0), u(NULL), B(NULL),
  computeAPtr(NULL), computeFPtr(NULL), computeUPtr(NULL), computeBPtr(NULL),
  isAAllocatedIn(true), isFAllocatedIn(false),  isUAllocatedIn(false), isBAllocatedIn(false)
{
  IN("LinearDS::LinearDS - XML constructor\n");
  if (dsXML != NULL)
  {
    DSType = LDS;

    // pointer to xml
    LinearDSXML * ldsxml = (static_cast <LinearDSXML*>(dsxml));

    isPlugin.resize(4, false);

    // --- vector and matrix members memory allocation ---
    // (only those specific to LinearDS) and values loading
    string plugin;

    A = new SiconosMatrix(n, n);
    if (ldsxml->isAPlugin())
    {
      plugin = ldsxml->getAPlugin();
      setComputeFFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
      isPlugin[0] = true;
    }
    else
      *A = ldsxml->getA();

    if (ldsxml->hasF())
    {
      f = new SimpleVector(n);
      isFAllocatedIn = true;
      if (ldsxml->isFPlugin())
      {
        plugin = ldsxml->getFPlugin();
        setComputeFFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
        isPlugin[1] = true;
      }
      else
        *f = ldsxml->getFVector();
    }

    if (ldsxml->hasU())
    {
      if (ldsxml->isUPlugin())
      {
        if (!ldsxml->hasUSize())
          RuntimeException::selfThrow("LinearDS - xml constructor needs uSize in xml file!!");
        uSize = ldsxml->getUSize();
        u = new SimpleVector(uSize);
        isUAllocatedIn = true;
        plugin = ldsxml->getUPlugin();
        setComputeUFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
        isPlugin[2] = true;
      }
      else
      {
        uSize = ldsxml->getUVector().size();
        u = new SimpleVector(uSize);
        isUAllocatedIn = true;
        *u = ldsxml->getUVector();
      }
    }

    if (ldsxml->hasB())
    {
      if (ldsxml->isBPlugin())
      {
        if (!ldsxml->hasUSize())
          RuntimeException::selfThrow("LinearDS - xml constructor needs uSize in xml file!!");
        uSize = ldsxml->getUSize();
        B = new SiconosMatrix(n, uSize);
        isBAllocatedIn = true;
        plugin = ldsxml->getBPlugin();
        setComputeBFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
        isPlugin[3] = true;
      }
      else
      {
        uSize = ldsxml->getUVector().size();
        B = new SiconosMatrix(n, uSize);
        isBAllocatedIn = true;
        *B = ldsxml->getB();
      }
    }
  }
  else
    RuntimeException::selfThrow("LinearDS - xml constructor, xml file = NULL");

  OUT("LinearDS::LinearDS - XML constructor\n");

}

// For the following constructors, only A is required. If necessary f, u and B can be defined thanks
// to set or setCompute, depending on they are plugins or not.

// From a minimum set of data, A from a plugin
LinearDS::LinearDS(const int& newNumber, const unsigned int& newN, const SiconosVector& newX0,
                   const string& pluginPath, const string& functionName):
  DynamicalSystem(newNumber, newN, newX0, "BasicPlugin:vectorField"),
  A(NULL), f(NULL), uSize(0), u(NULL), B(NULL),
  computeAPtr(NULL), computeFPtr(NULL), computeUPtr(NULL), computeBPtr(NULL),
  isAAllocatedIn(true), isFAllocatedIn(false),  isUAllocatedIn(false), isBAllocatedIn(false)
{
  DSType = LDS;
  isPlugin.resize(4, false);

  A = new SiconosMatrix(n, n);
  setComputeAFunction(pluginPath, functionName);
  isPlugin[0] = true;
}

// From a minimum set of data, A from a given matrix
LinearDS::LinearDS(const int& newNumber, const SiconosVector& newX0,
                   const SiconosMatrix& newA):
  DynamicalSystem(newNumber, newA.size(0), newX0, "BasicPlugin:vectorField"),
  A(NULL), f(NULL), uSize(0), u(NULL), B(NULL),
  computeAPtr(NULL), computeFPtr(NULL), computeUPtr(NULL), computeBPtr(NULL),
  isAAllocatedIn(true), isFAllocatedIn(false), isUAllocatedIn(false), isBAllocatedIn(false)
{
  if (newA.size(0) != n || newA.size(1) != n)
    RuntimeException::selfThrow("LinearDS - constructor(3): inconsistent dimensions with problem size for input matrix A");

  DSType = LDS;
  A = new SiconosMatrix(n, n);
  *A = newA;
  isPlugin.resize(4, false);
}

LinearDS::~LinearDS()
{
  if (isAAllocatedIn)
  {
    delete A;
    A = NULL ;
  }
  if (isBAllocatedIn)
  {
    delete B;
    B = NULL;
  }
  if (isUAllocatedIn)
  {
    delete u;
    u = NULL;
  }
  if (isFAllocatedIn)
  {
    delete f;
    f = NULL;
  }
}

void LinearDS::setAPtr(SiconosMatrix *newPtr)
{
  if (isAAllocatedIn) delete A;
  A = newPtr;
  isAAllocatedIn = false;
  isPlugin[0] = false;
}

void LinearDS::setF(const SimpleVector& newValue)
{
  if (newValue.size() != n)
    RuntimeException::selfThrow("LinearDS - setF: inconsistent dimensions with problem size for input vector f");

  if (f == NULL)
  {
    f = new SimpleVector(n);
    isFAllocatedIn = true;
  }
  *f = newValue;
  isPlugin[1] = false;
}

void LinearDS::setFPtr(SimpleVector *newPtr)
{
  if (isFAllocatedIn) delete f;
  f = newPtr;
  isFAllocatedIn = false;
  isPlugin[1] = false;
}

void LinearDS::setUSize(const unsigned int& newUSize)
{
  uSize = newUSize;
  if (isUAllocatedIn) delete u;
  if (isBAllocatedIn) delete B;
  u = new SimpleVector(uSize);
  B = new SiconosMatrix(n, uSize);
}

void LinearDS::setU(const SimpleVector& newValue)
{
  // check if u has already been allocated and if so, if sises are consistent.

  if (u == NULL)
  {
    uSize = newValue.size();
    u = new SimpleVector(uSize);
    isUAllocatedIn = true;
  }
  else if (uSize != newValue.size())
    RuntimeException::selfThrow("LinearDS - setU: inconsistent dimensions with already allocated u for new input vector");

  *u = newValue;
  isPlugin[2] = false;
}

void LinearDS::setUPtr(SimpleVector *newPtr)
{
  if (isUAllocatedIn) delete u;
  u = newPtr;
  uSize = newPtr->size();
  isUAllocatedIn = false;
  isPlugin[2] = false;
}

void LinearDS::setB(const SiconosMatrix& newValue)
{
  if (newValue.size(0) != n)
    RuntimeException::selfThrow("LinearDS - setF: inconsistent dimensions with problem size for input matrix B");

  if (B == NULL)
  {
    uSize = newValue.size(1);
    B = new SiconosMatrix(n, uSize);
    isBAllocatedIn = true;
  }
  else if (uSize != newValue.size(1))
    RuntimeException::selfThrow("LinearDS - setB: inconsistent dimensions with already allocated B for new input matrix");

  *B = newValue;
  isPlugin[3] = false;
}

void LinearDS::setBPtr(SiconosMatrix *newPtr)
{
  if (isBAllocatedIn) delete B;
  B = newPtr;
  uSize = newPtr->size(1);
  isBAllocatedIn = false;
  isPlugin[3] = false;
}

void LinearDS::setComputeAFunction(const string& pluginPath, const string& functionName)
{
  IN("LinearDS::setComputeAFunction\n");
  cShared.setFunction(&computeAPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  AFunctionName = plugin + ":" + functionName;
  isPlugin[0] = true;
  OUT("LinearDS::setComputeAFunction\n");

}
void LinearDS::setComputeFFunction(const string& pluginPath, const string& functionName)
{
  cShared.setFunction(&computeFPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  fFunctionName = plugin + ":" + functionName;
  isPlugin[1] = true;
}

void LinearDS::setComputeUFunction(const string& pluginPath, const string& functionName)
{
  cShared.setFunction(&computeUPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  uFunctionName = plugin + ":" + functionName;
  isPlugin[2] = true;
}
void LinearDS::setComputeBFunction(const string& pluginPath, const string& functionName)
{
  IN("LinearDS::setComputeBFunction\n");
  cShared.setFunction(&computeBPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  BFunctionName = plugin + ":" + functionName;
  isPlugin[3] = true;
  OUT("LinearDS::setComputeBFunction\n");

}

void LinearDS::computeA(const double& time)
{
  if (computeAPtr == NULL)
    RuntimeException::selfThrow("computeA() is not linked to a plugin function");
  computeAPtr(&n, &(*A)(0, 0), &time);
}

void LinearDS::computeF(const double& time)
{
  if (computeFPtr == NULL)
    RuntimeException::selfThrow("computeF() is not linked to a plugin function");

  unsigned int size = f->size();
  computeFPtr(&size, &(*f)(0), &time);
}

void LinearDS::computeU(const double& time)
{
  if (computeUPtr == NULL)
    RuntimeException::selfThrow("computeU() is not linked to a plugin function");

  unsigned int size = u->size();
  computeUPtr(&size, &(*u)(0), &time);
}

void LinearDS::computeB(const double& time)
{
  if (computeBPtr == NULL)
    RuntimeException::selfThrow("computeB() is not linked to a plugin function");
  computeBPtr(&n, &uSize, &(*B)(0, 0), &time);
}

void LinearDS::display() const
{
  DynamicalSystem::display();
  cout << "=== Linear system display ===" << endl;
  cout << "- A " << endl;
  if (A != NULL) A->display();
  else cout << "-> NULL" << endl;
  cout << "- f " << endl;
  if (f != NULL) f->display();
  else cout << "-> NULL" << endl;
  cout << "- u " << endl;
  if (u != NULL) u->display();
  else cout << "-> NULL" << endl;
  cout << "- B " << endl;
  if (B != NULL) B->display();
  else cout << "-> NULL" << endl;
  cout << "=============================" << endl;
}

void LinearDS::saveDSToXML()
{
  IN("LinearDS::saveDSToXML\n");

  //--- Common data ---
  saveDSDataToXML();
  // --- other data ---
  if (dsxml != NULL)
  {
    dsxml->setN(n);
    static_cast<LinearDSXML*>(dsxml)->setA(*A);
    static_cast<LinearDSXML*>(dsxml)->setB(*B);

    // u
    if (!(static_cast <LinearDSXML*>(dsxml))->isUPlugin())
    {
      static_cast<LinearDSXML*>(dsxml)->setUVector(*u);
    }

    // f
    if (!(static_cast <LinearDSXML*>(dsxml))->isFPlugin())
    {
      static_cast<LinearDSXML*>(dsxml)->setFVector(*f);
    }
  }
  else RuntimeException::selfThrow("LinearDS::saveDSToXML - The DSXML object doesn't exists");

  OUT("LinearDS::saveDSToXML\n");
}

LinearDS* LinearDS::convert(DynamicalSystem* ds)
{
  cout << "LinearDS::convert (DynamicalSystem* ds)" << endl;
  LinearDS* lsds = dynamic_cast<LinearDS*>(ds);
  return lsds;
}

// Default constructor
LinearDS::LinearDS():
  DynamicalSystem(NULL), A(NULL), f(NULL), uSize(0), u(NULL), B(NULL),
  computeAPtr(NULL), computeFPtr(NULL), computeUPtr(NULL), computeBPtr(NULL),
  isAAllocatedIn(true), isFAllocatedIn(false), isUAllocatedIn(false), isBAllocatedIn(false)
{
  IN("LinearDS::LinearDS - Default constructor\n");
  DSType = LDS;
  isPlugin.resize(4, false);
  OUT("LinearDS::LinearDS - Default constructor\n");
}

// \todo link to vectorField function
/*void LinearDS::computeVectorField(const double& time)
{
  if(isPlugin[0])
    computeA(time);
  // if necessary, compute the required functions from plugins
  // Warning: B != NULL means u!=NULL
  if(isPlugin[0])
    computeA(time)

  if(f!=NULL)
    {
      if(isPlugin[1])
  computeF(time);
    }
  if(B!=NULL)
    {
      if(isPlugin[2])
  computeU(time);
      if(isPlugin[3])
  computeB(time);
    }

  // compute vectorField
  if(f!=NULL && B==NULL)
    vectorField = *A * *x + *f;
  else if(f!=NULL && B!=NULL)
    vectorField = *A * *x + *f + *B * *u;
  else if(f==NULL && B!= NULL)
    vectorField = *A * *x + *B * *u;
  else
    vectorField = *A * *x;

}
*/
