#include "LinearDS.h"
using namespace std;

// --- Constructors ---

// From xml file (newNsds is optional)
LinearDS::LinearDS(DynamicalSystemXML * dsXML, NonSmoothDynamicalSystem* newNsds):
  DynamicalSystem(dsXML, newNsds), A(NULL), b(NULL),
  AFunctionName("none"), bFunctionName("none"),
  computeAPtr(NULL), computeBPtr(NULL),
  isAAllocatedIn(true), isBAllocatedIn(false)
{
  IN("LinearDS::LinearDS - XML constructor\n");
  if (dsXML != NULL)
  {
    DSType = LDS;

    // pointer to xml
    LinearDSXML * ldsxml = (static_cast <LinearDSXML*>(dsxml));

    isLDSPlugin.resize(2, false);

    // --- vector and matrix members memory allocation ---
    // (only those specific to LinearDS) and values loading
    string plugin;

    A = new SiconosMatrix(n, n);
    if (ldsxml->isAPlugin())
    {
      plugin = ldsxml->getAPlugin();
      setComputeAFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
      isLDSPlugin[0] = true;
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
        isLDSPlugin[1] = true;
      }
      else
        *b = ldsxml->getBVector();
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
  A(NULL), b(NULL),
  AFunctionName("none"), bFunctionName("none"),
  computeAPtr(NULL), computeBPtr(NULL),
  isAAllocatedIn(true), isBAllocatedIn(false)
{
  DSType = LDS;
  isLDSPlugin.resize(2, false);

  A = new SiconosMatrix(n, n);
  setComputeAFunction(pluginPath, functionName);
  isLDSPlugin[0] = true;
}

// From a minimum set of data, A from a given matrix
LinearDS::LinearDS(const int& newNumber, const SiconosVector& newX0,
                   const SiconosMatrix& newA):
  DynamicalSystem(newNumber, newA.size(0), newX0),
  A(NULL), b(NULL),
  AFunctionName("none"), bFunctionName("none"),
  computeAPtr(NULL), computeBPtr(NULL),
  isAAllocatedIn(true), isBAllocatedIn(false)
{
  if (newA.size(0) != n || newA.size(1) != n)
    RuntimeException::selfThrow("LinearDS - constructor(3): inconsistent dimensions with problem size for input matrix A");

  DSType = LDS;
  A = new SiconosMatrix(n, n);
  *A = newA;
  isLDSPlugin.resize(2, false);
}

// Copy constructor
LinearDS::LinearDS(const DynamicalSystem & newDS):
  DynamicalSystem(newDS), A(NULL), b(NULL),
  AFunctionName("none"), bFunctionName("none"),
  computeAPtr(NULL), computeBPtr(NULL),
  isAAllocatedIn(true), isBAllocatedIn(false)

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

  isLDSPlugin = lds->getIsLDSPlugin();
  string pluginPath, functionName;
  if (isLDSPlugin[0])
  {
    AFunctionName = lds -> getAFunctionName();
    functionName = cShared.getPluginFunctionName(AFunctionName);
    pluginPath  = cShared.getPluginName(AFunctionName);
    setComputeAFunction(pluginPath, functionName);
  }
  if (isLDSPlugin[1])
  {
    bFunctionName = lds -> getBFunctionName();
    functionName = cShared.getPluginFunctionName(bFunctionName);
    pluginPath  = cShared.getPluginName(bFunctionName);
    setComputeBFunction(pluginPath, functionName);
  }
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
    delete b;
    b = NULL;
  }
}

void LinearDS::setAPtr(SiconosMatrix *newPtr)
{
  if (isAAllocatedIn) delete A;
  A = newPtr;
  isAAllocatedIn = false;
  isLDSPlugin[0] = false;
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
  isLDSPlugin[1] = false;
}

void LinearDS::setBPtr(SimpleVector *newPtr)
{
  if (isBAllocatedIn) delete b;
  b = newPtr;
  isBAllocatedIn = false;
  isLDSPlugin[1] = false;
}

void LinearDS::setComputeAFunction(const string& pluginPath, const string& functionName)
{
  IN("LinearDS::setComputeAFunction\n");
  cShared.setFunction(&computeAPtr, pluginPath, functionName);
  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  AFunctionName = plugin + ":" + functionName;
  isLDSPlugin[0] = true;
  OUT("LinearDS::setComputeAFunction\n");

}
void LinearDS::setComputeBFunction(const string& pluginPath, const string& functionName)
{
  cShared.setFunction(&computeBPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  bFunctionName = plugin + ":" + functionName;
  isLDSPlugin[1] = true;
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

    // b
    if (!(static_cast <LinearDSXML*>(dsxml))->isBPlugin())
    {
      static_cast<LinearDSXML*>(dsxml)->setBVector(*b);
    }
  }
  else RuntimeException::selfThrow("LinearDS::saveDSToXML - The DynamicalSystemXML object doesn't exists");

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
  DynamicalSystem(NULL), A(NULL), b(NULL),
  computeAPtr(NULL), computeBPtr(NULL),
  isAAllocatedIn(true), isBAllocatedIn(false)
{
  IN("LinearDS::LinearDS - Default constructor\n");
  DSType = LDS;
  isLDSPlugin.resize(2, false);
  OUT("LinearDS::LinearDS - Default constructor\n");
}

// \todo link to vectorField function
/*void LinearDS::computeVectorField(const double& time)
{
  if(isLDSPlugin[0])
    computeA(time);
  // if necessary, compute the required functions from plugins
  // Warning: T != NULL means u!=NULL
  if(isLDSPlugin[0])
    computeA(time)

  if(b!=NULL)
    {
      if(isLDSPlugin[1])
  computeB(time);
    }
  if(T!=NULL)
    {
      if(isPlugin[0])
  computeU(time);
      if(isPlugin[1])
  computeT(time);
    }

  // compute vectorField
  if(b!=NULL && T==NULL)
    vectorField = *A * *x + *b;
  else if(b!=NULL && T!=NULL)
    vectorField = *A * *x + *b + *T * *u;
  else if(b==NULL && T!= NULL)
    vectorField = *A * *x + *T * *u;
  else
    vectorField = *A * *x;

}
*/
