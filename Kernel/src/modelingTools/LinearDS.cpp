#include "LinearDS.h"
using namespace std;

// --- Constructors ---

// From xml file (newNsds is optional)
LinearDS::LinearDS(DSXML * dsXML, NonSmoothDynamicalSystem* newNsds):
  DynamicalSystem(dsXML, newNsds), A(NULL), b(NULL), uSize(0), u(NULL), E(NULL),
  AFunctionName("none"), bFunctionName("none"), uFunctionName("none"), EFunctionName("none"),
  computeAPtr(NULL), computeBPtr(NULL), computeUPtr(NULL), computeEPtr(NULL),
  isAAllocatedIn(true), isBAllocatedIn(false),  isUAllocatedIn(false), isEAllocatedIn(false)
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
      setComputeAFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
      isPlugin[0] = true;
    }
    else
      *A = ldsxml->getA();

    if (ldsxml->hasB())
    {
      b = new SimpleVector(n);
      isBAllocatedIn = true;
      if (ldsxml->isBPlugin())
      {
        plugin = ldsxml->getBPlugin();
        setComputeBFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
        isPlugin[1] = true;
      }
      else
        *b = ldsxml->getBVector();
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

    if (ldsxml->hasE())
    {
      if (ldsxml->isEPlugin())
      {
        if (!ldsxml->hasUSize())
          RuntimeException::selfThrow("LinearDS - xml constructor needs uSize in xml file!!");
        uSize = ldsxml->getUSize();
        E = new SiconosMatrix(n, uSize);
        isEAllocatedIn = true;
        plugin = ldsxml->getEPlugin();
        setComputeEFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
        isPlugin[3] = true;
      }
      else
      {
        uSize = ldsxml->getE().size(1);
        E = new SiconosMatrix(n, uSize);
        isEAllocatedIn = true;
        *E = ldsxml->getE();
      }
    }
  }
  else
    RuntimeException::selfThrow("LinearDS - xml constructor, xml file = NULL");

  OUT("LinearDS::LinearDS - XML constructor\n");

}

// For the following constructors, only A is required. If necessary b, u and E can be defined thanks
// to set or setCompute, depending on they are plugins or not.

// From a minimum set of data, A from a plugin
LinearDS::LinearDS(const int& newNumber, const unsigned int& newN, const SiconosVector& newX0,
                   const string& pluginPath, const string& functionName):
  DynamicalSystem(newNumber, newN, newX0, "BasicPlugin:vectorField"),
  A(NULL), b(NULL), uSize(0), u(NULL), E(NULL),
  AFunctionName("none"), bFunctionName("none"), uFunctionName("none"), EFunctionName("none"),
  computeAPtr(NULL), computeBPtr(NULL), computeUPtr(NULL), computeEPtr(NULL),
  isAAllocatedIn(true), isBAllocatedIn(false),  isUAllocatedIn(false), isEAllocatedIn(false)
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
  A(NULL), b(NULL), uSize(0), u(NULL), E(NULL),
  AFunctionName("none"), bFunctionName("none"), uFunctionName("none"), EFunctionName("none"),
  computeAPtr(NULL), computeBPtr(NULL), computeUPtr(NULL), computeEPtr(NULL),
  isAAllocatedIn(true), isBAllocatedIn(false), isUAllocatedIn(false), isEAllocatedIn(false)
{
  if (newA.size(0) != n || newA.size(1) != n)
    RuntimeException::selfThrow("LinearDS - constructor(3): inconsistent dimensions with problem size for input matrix A");

  DSType = LDS;
  A = new SiconosMatrix(n, n);
  *A = newA;
  isPlugin.resize(4, false);
}

// Copy constructor
LinearDS::LinearDS(const DynamicalSystem & newDS):
  DynamicalSystem(newDS), A(NULL), b(NULL), uSize(0), u(NULL), E(NULL),
  AFunctionName("none"), bFunctionName("none"), uFunctionName("none"), EFunctionName("none"),
  computeAPtr(NULL), computeBPtr(NULL), computeUPtr(NULL), computeEPtr(NULL),
  isAAllocatedIn(true), isBAllocatedIn(false), isUAllocatedIn(false), isEAllocatedIn(false)

{
  if (newDS.getType() != LDS)
    RuntimeException::selfThrow("LinearDS - copy constructor: try to copy into a LinearDS a DS of type: " + newDS.getType());

  DSType = LDS;

  // convert newDS to linearDS by keeping const options
  const LinearDS * lds = static_cast<const LinearDS*>(&newDS);

  A = new SiconosMatrix(lds->getA());
  isAAllocatedIn = true;
  if (lds->getBPtr() != NULL)
  {
    b = new SimpleVector(lds->getB());
    isBAllocatedIn = true;
  }
  if (lds->getUPtr() != NULL)
  {
    u = new SimpleVector(lds->getU());
    uSize = u->size();
    isUAllocatedIn = true;
  }
  if (lds->getEPtr() != NULL)
  {
    E = new SiconosMatrix(lds->getE());
    isEAllocatedIn = true;
  }

  isPlugin = lds->getIsPlugin();
  string pluginPath, functionName;
  if (isPlugin[0])
  {
    AFunctionName = lds -> getAFunctionName();
    functionName = cShared.getPluginFunctionName(AFunctionName);
    pluginPath  = cShared.getPluginName(AFunctionName);
    setComputeAFunction(pluginPath, functionName);
  }
  if (isPlugin[1])
  {
    bFunctionName = lds -> getBFunctionName();
    functionName = cShared.getPluginFunctionName(bFunctionName);
    pluginPath  = cShared.getPluginName(bFunctionName);
    setComputeBFunction(pluginPath, functionName);
  }
  if (isPlugin[2])
  {
    uFunctionName = lds -> getUFunctionName();
    functionName = cShared.getPluginFunctionName(uFunctionName);
    pluginPath  = cShared.getPluginName(uFunctionName);
    setComputeUFunction(pluginPath, functionName);
  }

  if (isPlugin[3])
  {
    EFunctionName = lds -> getEFunctionName();
    functionName = cShared.getPluginFunctionName(EFunctionName);
    pluginPath  = cShared.getPluginName(EFunctionName);
    setComputeEFunction(pluginPath, functionName);
  }
}

LinearDS::~LinearDS()
{
  if (isAAllocatedIn)
  {
    delete A;
    A = NULL ;
  }
  if (isEAllocatedIn)
  {
    delete E;
    E = NULL;
  }
  if (isUAllocatedIn)
  {
    delete u;
    u = NULL;
  }
  if (isBAllocatedIn)
  {
    delete b;
    b = NULL;
  }
}

void LinearDS::setAPtr(SiconosMatrix *newPtr)
{
  if (isAAllocatedIn) delete A;
  A = newPtr;
  isAAllocatedIn = false;
  isPlugin[0] = false;
}

void LinearDS::setB(const SimpleVector& newValue)
{
  if (newValue.size() != n)
    RuntimeException::selfThrow("LinearDS - setB: inconsistent dimensions with problem size for input vector b");

  if (b == NULL)
  {
    b = new SimpleVector(n);
    isBAllocatedIn = true;
  }
  *b = newValue;
  isPlugin[1] = false;
}

void LinearDS::setBPtr(SimpleVector *newPtr)
{
  if (isBAllocatedIn) delete b;
  b = newPtr;
  isBAllocatedIn = false;
  isPlugin[1] = false;
}

void LinearDS::setUSize(const unsigned int& newUSize)
{
  uSize = newUSize;
  if (isUAllocatedIn) delete u;
  if (isEAllocatedIn) delete E;
  u = new SimpleVector(uSize);
  E = new SiconosMatrix(n, uSize);
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

void LinearDS::setE(const SiconosMatrix& newValue)
{
  if (newValue.size(0) != n)
    RuntimeException::selfThrow("LinearDS - setE: inconsistent dimensions with problem size for input matrix E");

  if (E == NULL)
  {
    uSize = newValue.size(1);
    E = new SiconosMatrix(n, uSize);
    isEAllocatedIn = true;
  }
  else if (uSize != newValue.size(1))
    RuntimeException::selfThrow("LinearDS - setE: inconsistent dimensions with already allocated E for new input matrix");

  *E = newValue;
  isPlugin[3] = false;
}

void LinearDS::setEPtr(SiconosMatrix *newPtr)
{
  if (isEAllocatedIn) delete E;
  E = newPtr;
  uSize = newPtr->size(1);
  isEAllocatedIn = false;
  isPlugin[3] = false;
}

void LinearDS::setComputeAFunction(const string& pluginPath, const string& functionName)
{
  IN("LinearDS::setComputeAFunction\n");
  cShared.setFunction(&computeAPtr, pluginPath, functionName);
  cout << " COMPUTEA PTR " << computeAPtr << endl;
  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  AFunctionName = plugin + ":" + functionName;
  isPlugin[0] = true;
  OUT("LinearDS::setComputeAFunction\n");

}
void LinearDS::setComputeBFunction(const string& pluginPath, const string& functionName)
{
  cShared.setFunction(&computeBPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  bFunctionName = plugin + ":" + functionName;
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
void LinearDS::setComputeEFunction(const string& pluginPath, const string& functionName)
{
  IN("LinearDS::setComputeEFunction\n");
  cShared.setFunction(&computeEPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  EFunctionName = plugin + ":" + functionName;
  isPlugin[3] = true;
  OUT("LinearDS::setComputeEFunction\n");

}

void LinearDS::computeA(const double& time)
{
  if (computeAPtr == NULL)
    RuntimeException::selfThrow("computeA() is not linked to a plugin function");
  computeAPtr(&n, &(*A)(0, 0), &time);
}

void LinearDS::computeB(const double& time)
{
  if (computeBPtr == NULL)
    RuntimeException::selfThrow("computeB() is not linked to a plugin function");

  unsigned int size = b->size();
  computeBPtr(&size, &(*b)(0), &time);
}

void LinearDS::computeU(const double& time)
{
  if (computeUPtr == NULL)
    RuntimeException::selfThrow("computeU() is not linked to a plugin function");

  unsigned int size = u->size();
  computeUPtr(&size, &(*u)(0), &time);
}

void LinearDS::computeE(const double& time)
{
  if (computeEPtr == NULL)
    RuntimeException::selfThrow("computeE() is not linked to a plugin function");
  computeEPtr(&n, &uSize, &(*E)(0, 0), &time);
}

void LinearDS::display() const
{
  DynamicalSystem::display();
  cout << "=== Linear system display ===" << endl;
  cout << "- A " << endl;
  if (A != NULL) A->display();
  else cout << "-> NULL" << endl;
  cout << "- b " << endl;
  if (b != NULL) b->display();
  else cout << "-> NULL" << endl;
  cout << "- u " << endl;
  if (u != NULL) u->display();
  else cout << "-> NULL" << endl;
  cout << "- E " << endl;
  if (E != NULL) E->display();
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
    static_cast<LinearDSXML*>(dsxml)->setE(*E);

    // u
    if (!(static_cast <LinearDSXML*>(dsxml))->isUPlugin())
    {
      static_cast<LinearDSXML*>(dsxml)->setUVector(*u);
    }

    // b
    if (!(static_cast <LinearDSXML*>(dsxml))->isBPlugin())
    {
      static_cast<LinearDSXML*>(dsxml)->setBVector(*b);
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
  DynamicalSystem(NULL), A(NULL), b(NULL), uSize(0), u(NULL), E(NULL),
  computeAPtr(NULL), computeBPtr(NULL), computeUPtr(NULL), computeEPtr(NULL),
  isAAllocatedIn(true), isBAllocatedIn(false), isUAllocatedIn(false), isEAllocatedIn(false)
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
  // Warning: E != NULL means u!=NULL
  if(isPlugin[0])
    computeA(time)

  if(b!=NULL)
    {
      if(isPlugin[1])
  computeB(time);
    }
  if(E!=NULL)
    {
      if(isPlugin[2])
  computeU(time);
      if(isPlugin[3])
  computeE(time);
    }

  // compute vectorField
  if(b!=NULL && E==NULL)
    vectorField = *A * *x + *b;
  else if(b!=NULL && E!=NULL)
    vectorField = *A * *x + *b + *E * *u;
  else if(b==NULL && E!= NULL)
    vectorField = *A * *x + *E * *u;
  else
    vectorField = *A * *x;

}
*/
