#include "LinearSystemDS.h"
#include "check.h"

// --- Constructors ---

// From xml file
LinearSystemDS::LinearSystemDS(DSXML * dsXML): DynamicalSystem(dsXML), A(0), B(0), uSize(1), u(0), f(0)
{
  IN("LinearSystemDS::LinearSystemDS - XML constructor\n");
  if (dsXML != 0)
  {
    this->DSType = LDS;
    // size of u
    // \warning: nothing available in LinearSystemDSXML to load uSize
    //if ( this->dsxml->hasUSize()==true) this->uSize = this->dsxml->getUSize();
    uSize = this->n;
    // --- vector and matrix members memory allocation ---
    // only those specific to LinearDS
    this->A = new SiconosMatrix(this->n, this->n);
    this->B = new SiconosMatrix(this->uSize, this->uSize);
    this->u = new SimpleVector(this->uSize);
    this->f = new SimpleVector(this->n);

    // --- xml loading and plugins ---
    string plugin;
    *this->A = (static_cast <LinearSystemDSXML*>(this->dsxml))->getA();
    *this->B = (static_cast <LinearSystemDSXML*>(this->dsxml))->getB();
    // u
    if ((static_cast <LinearSystemDSXML*>(this->dsxml))->isUPlugin())
    {
      plugin = (static_cast <LinearSystemDSXML*>(this->dsxml))->getUPlugin();
      this->setComputeUFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
    }
    else
    {
      this->setComputeUFunction("BasicPlugin.so", "computeU");
      *this->u = (static_cast <LinearSystemDSXML*>(this->dsxml))->getUVector();
    }
    // f
    if ((static_cast <LinearSystemDSXML*>(this->dsxml))->isFPlugin())
    {
      plugin = (static_cast <LinearSystemDSXML*>(this->dsxml))->getFPlugin();
      this->setComputeFFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
    }
    else
    {
      this->setComputeFFunction("BasicPlugin.so", "computeF");
      *this->f = (static_cast <LinearSystemDSXML*>(this->dsxml))->getFVector();
    }
  }
  else
  {
    cout << "LinearSystemDS::LinearSystemDS - DSXML paramater musn't be 0" << endl;
  }
  OUT("LinearSystemDS::LinearSystemDS - XML constructor\n");

}

// From a minimum set of data
LinearSystemDS::LinearSystemDS(int number, int n, SiconosVector* x0): DynamicalSystem(number, n, x0, "BasicPlugin:vectorField")
{
  this->DSType = LDS;
  this->A = new SiconosMatrix(this->n, this->n);
  this->B = new SiconosMatrix(this->n, this->n);
  this->u = new SimpleVector(this->n);
  this->f = new SimpleVector(this->n);

  // --- plugins ---
  this->setComputeFFunction("BasicPlugin.so", "computeF");
  this->setComputeUFunction("BasicPlugin.so", "computeU");
}

LinearSystemDS::~LinearSystemDS()
{
  delete A;
  A = 0 ;
  delete B;
  B = 0;
  delete u;
  u = 0;
  delete f;
  f = 0;
}

void LinearSystemDS::setComputeFFunction(const string& pluginPath, const string& functionName)
{
  cShared.setFunction(&computeFPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->fFunctionName = plugin + ":" + functionName;
}

void LinearSystemDS::setComputeUFunction(const string& pluginPath, const string& functionName)
{
  cShared.setFunction(&computeUPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->uFunctionName = plugin + ":" + functionName;
}

void LinearSystemDS::computeF(const double& time)
{
  if (computeFPtr == 0)
    RuntimeException::selfThrow("computeF() is not linked to a plugin function");

  int size = f->size();
  // const_cast to be deleted: problem const in C function signature? To see ...
  this->computeFPtr(&size, &(*f)(0), const_cast<double*>(&time));
}

void LinearSystemDS::computeU(const double& time)
{
  if (computeUPtr == 0)
    RuntimeException::selfThrow("computeU() is not linked to a plugin function");

  int size = u->size();
  // const_cast to be deleted: problem const in C function signature? To see ...
  this->computeUPtr(&size, &(*u)(0), const_cast<double*>(&time));
}

void LinearSystemDS::display() const
{
  cout << "-----------------------------------------------------" << endl;
  cout << "____ data of the LinearSystemDS " << endl;
  DynamicalSystem::display();
  cout << "| A " << endl;
  if (A != 0) this->A->display();
  else cout << "-> 0" << endl;
  cout << "| B " << endl;
  if (B != 0) this->B->display();
  else cout << "-> 0" << endl;
  cout << "| u " << endl;
  if (u != 0) this->u->display();
  else cout << "-> 0" << endl;
  cout << "| f " << endl;
  if (f != 0) this->f->display();
  else cout << "-> 0" << endl;
  cout << "-----------------------------------------------------" << endl << endl;
}

void LinearSystemDS::saveDSToXML()
{
  IN("LinearSystemDS::saveDSToXML\n");

  //--- Common data ---
  saveDSDataToXML();
  // --- other data ---
  if (this->dsxml != 0)
  {
    this->dsxml->setN(this->n);
    static_cast<LinearSystemDSXML*>(this->dsxml)->setA(this->A);
    static_cast<LinearSystemDSXML*>(this->dsxml)->setB(this->B);

    // u
    if (!(static_cast <LinearSystemDSXML*>(this->dsxml))->isUPlugin())
    {
      static_cast<LinearSystemDSXML*>(this->dsxml)->setUVector(this->u);
    }

    // f
    if (!(static_cast <LinearSystemDSXML*>(this->dsxml))->isFPlugin())
    {
      static_cast<LinearSystemDSXML*>(this->dsxml)->setFVector(this->f);
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
LinearSystemDS::LinearSystemDS(): DynamicalSystem(), A(0), B(0), uSize(0), u(0), f(0)
{
  IN("LinearSystemDS::LinearSystemDS - Default constructor\n");
  this->DSType = LDS;
  this->setComputeFFunction("BasicPlugin.so", "computeF");
  this->setComputeUFunction("BasicPlugin.so", "computeU");
  OUT("LinearSystemDS::LinearSystemDS - Default constructor\n");
}
