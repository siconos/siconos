#include "LinearSystemDS.h"
using namespace std;

// --- Constructors ---

// From xml file
LinearSystemDS::LinearSystemDS(DSXML * dsXML): DynamicalSystem(dsXML), A(NULL), B(NULL), uSize(1), u(NULL), f(NULL)
{
  IN("LinearSystemDS::LinearSystemDS - XML constructor\n");
  if (dsXML != NULL)
  {
    DSType = LDS;
    // size of u
    // \warning: nothing available in LinearSystemDSXML to load uSize
    //if ( dsxml->hasUSize()==true) uSize = dsxml->getUSize();
    uSize = n;
    // --- vector and matrix members memory allocation ---
    // only those specific to LinearDS
    A = new SiconosMatrix(n, n);
    B = new SiconosMatrix(uSize, uSize);
    u = new SimpleVector(uSize);
    f = new SimpleVector(n);

    // --- xml loading and plugins ---
    string plugin;
    *A = (static_cast <LinearSystemDSXML*>(dsxml))->getA();
    *B = (static_cast <LinearSystemDSXML*>(dsxml))->getB();
    // u
    if ((static_cast <LinearSystemDSXML*>(dsxml))->isUPlugin())
    {
      plugin = (static_cast <LinearSystemDSXML*>(dsxml))->getUPlugin();
      setComputeUFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else
    {
      setComputeUFunction("BasicPlugin.so", "computeU");
      *u = (static_cast <LinearSystemDSXML*>(dsxml))->getUVector();
    }
    // f
    if ((static_cast <LinearSystemDSXML*>(dsxml))->isFPlugin())
    {
      plugin = (static_cast <LinearSystemDSXML*>(dsxml))->getFPlugin();
      setComputeFFunction(cShared.getPluginName(plugin), cShared.getPluginFunctionName(plugin));
    }
    else
    {
      setComputeFFunction("BasicPlugin.so", "computeF");
      *f = (static_cast <LinearSystemDSXML*>(dsxml))->getFVector();
    }
  }
  else
  {
    cout << "LinearSystemDS::LinearSystemDS - DSXML paramater musn't be NULL" << endl;
  }
  OUT("LinearSystemDS::LinearSystemDS - XML constructor\n");

}

// From a minimum set of data
LinearSystemDS::LinearSystemDS(int number, int n, SiconosVector* x0): DynamicalSystem(number, n, x0, "BasicPlugin:vectorField")
{
  DSType = LDS;
  A = new SiconosMatrix(n, n);
  B = new SiconosMatrix(n, n);
  u = new SimpleVector(n);
  f = new SimpleVector(n);

  // --- plugins ---
  setComputeFFunction("BasicPlugin.so", "computeF");
  setComputeUFunction("BasicPlugin.so", "computeU");
}

LinearSystemDS::~LinearSystemDS()
{
  delete A;
  A = NULL ;
  delete B;
  B = NULL;
  delete u;
  u = NULL;
  delete f;
  f = NULL;
}

void LinearSystemDS::setComputeFFunction(const string& pluginPath, const string& functionName)
{
  cShared.setFunction(&computeFPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  fFunctionName = plugin + ":" + functionName;
}

void LinearSystemDS::setComputeUFunction(const string& pluginPath, const string& functionName)
{
  cShared.setFunction(&computeUPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  uFunctionName = plugin + ":" + functionName;
}

void LinearSystemDS::computeF(const double& time)
{
  if (computeFPtr == NULL)
    RuntimeException::selfThrow("computeF() is not linked to a plugin function");

  int size = f->size();
  // const_cast to be deleted: problem const in C function signature? To see ...
  computeFPtr(&size, &(*f)(0), const_cast<double*>(&time));
}

void LinearSystemDS::computeU(const double& time)
{
  if (computeUPtr == NULL)
    RuntimeException::selfThrow("computeU() is not linked to a plugin function");

  int size = u->size();
  // const_cast to be deleted: problem const in C function signature? To see ...
  computeUPtr(&size, &(*u)(0), const_cast<double*>(&time));
}

void LinearSystemDS::display() const
{
  cout << "-----------------------------------------------------" << endl;
  cout << "____ data of the LinearSystemDS " << endl;
  DynamicalSystem::display();
  cout << "| A " << endl;
  if (A != NULL) A->display();
  else cout << "-> NULL" << endl;
  cout << "| B " << endl;
  if (B != NULL) B->display();
  else cout << "-> NULL" << endl;
  cout << "| u " << endl;
  if (u != NULL) u->display();
  else cout << "-> NULL" << endl;
  cout << "| f " << endl;
  if (f != NULL) f->display();
  else cout << "-> NULL" << endl;
  cout << "-----------------------------------------------------" << endl << endl;
}

void LinearSystemDS::saveDSToXML()
{
  IN("LinearSystemDS::saveDSToXML\n");

  //--- Common data ---
  saveDSDataToXML();
  // --- other data ---
  if (dsxml != NULL)
  {
    dsxml->setN(n);
    static_cast<LinearSystemDSXML*>(dsxml)->setA(A);
    static_cast<LinearSystemDSXML*>(dsxml)->setB(B);

    // u
    if (!(static_cast <LinearSystemDSXML*>(dsxml))->isUPlugin())
    {
      static_cast<LinearSystemDSXML*>(dsxml)->setUVector(u);
    }

    // f
    if (!(static_cast <LinearSystemDSXML*>(dsxml))->isFPlugin())
    {
      static_cast<LinearSystemDSXML*>(dsxml)->setFVector(f);
    }
  }
  else RuntimeException::selfThrow("LinearSystemDS::saveDSToXML - The DSXML object doesn't exists");

  OUT("LinearSystemDS::saveDSToXML\n");
}

LinearSystemDS* LinearSystemDS::convert(DynamicalSystem* ds)
{
  cout << "LinearSystemDS::convert (DynamicalSystem* ds)" << endl;
  LinearSystemDS* lsds = dynamic_cast<LinearSystemDS*>(ds);
  return lsds;
}

// Default constructor
LinearSystemDS::LinearSystemDS(): DynamicalSystem(), A(NULL), B(NULL), uSize(0), u(NULL), f(NULL)
{
  IN("LinearSystemDS::LinearSystemDS - Default constructor\n");
  DSType = LDS;
  setComputeFFunction("BasicPlugin.so", "computeF");
  setComputeUFunction("BasicPlugin.so", "computeU");
  OUT("LinearSystemDS::LinearSystemDS - Default constructor\n");
}
