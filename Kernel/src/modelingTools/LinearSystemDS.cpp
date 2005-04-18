#include "LinearSystemDS.h"
#include "check.h"


LinearSystemDS::LinearSystemDS(): DynamicalSystem()
{
  this->DSType = LSDS;
  this->init();
}

LinearSystemDS::LinearSystemDS(DSXML * dsXML)
{
  if (dsXML != NULL)
  {
    this->DSType = LSDS;
    DynamicalSystem::init();
    this->dsxml = dsXML;

    this->fillDSWithDSXML();
    this->linkDSXML();
  }
  else
  {
    cout << "LinearSystemDS::LinearSystemDS - DSXML paramater musn't be NULL" << endl;
  }
}

LinearSystemDS::LinearSystemDS(int number, int n, SiconosVector* x0)
{
  this->DSType = LSDS;
  this->number = number;
  this->n = n;
  *(this->x0) = *x0;
}

LinearSystemDS::~LinearSystemDS()
{}

SiconosMatrix* LinearSystemDS::getAPtr(void)
{
  return &this->A;
}

SiconosMatrix* LinearSystemDS::getBPtr(void)
{
  return &this->B;
}

SimpleVector* LinearSystemDS::getUPtr(void)
{
  return &this->u;
}

SimpleVector* LinearSystemDS::getFPtr(void)
{
  return &this->f;
}


/////////////////////////////


void LinearSystemDS::setComputeFFunction(string pluginPath, string functionName)
{
  cShared.setFunction(&computeFPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->fFunctionName = plugin + ":" + functionName;
}

void LinearSystemDS::setComputeUFunction(string pluginPath, string functionName)
{
  cShared.setFunction(&computeUPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->uFunctionName = plugin + ":" + functionName;
}

void LinearSystemDS::computeF(double time)
{
  if (computeFPtr == NULL)
    RuntimeException::selfThrow("computeF() is not linked to a plugin function");

  int size = f.size();
  this->computeFPtr(&size, &f(0), &time);
}

void LinearSystemDS::computeU(double time)
{
  if (computeUPtr == NULL)
    RuntimeException::selfThrow("computeU() is not linked to a plugin function");

  int size = u.size();
  this->computeUPtr(&size, &u(0), &time);
}

///////////////////////////////

void LinearSystemDS::fillDSWithDSXML()
{
  IN("LinearSystemDS::fillDSWithDSXML\n");

  string plugin;

  DynamicalSystem::fillDSWithDSXML();
  if (this->dsxml != NULL)
  {
    this->A = (static_cast <LinearSystemDSXML*>(this->dsxml))->getA();
    this->B = (static_cast <LinearSystemDSXML*>(this->dsxml))->getB();

    // u
    if ((static_cast <LinearSystemDSXML*>(this->dsxml))->isUPlugin())
    {
      plugin = (static_cast <LinearSystemDSXML*>(this->dsxml))->getUPlugin();
      this->setComputeUFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
    }
    else
    {
      this->u = (static_cast <LinearSystemDSXML*>(this->dsxml))->getUVector();
    }

    // f
    if ((static_cast <LinearSystemDSXML*>(this->dsxml))->isFPlugin())
    {
      plugin = (static_cast <LinearSystemDSXML*>(this->dsxml))->getFPlugin();
      this->setComputeFFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
    }
    else
    {
      this->f = (static_cast <LinearSystemDSXML*>(this->dsxml))->getFVector();
    }

  }
  else RuntimeException::selfThrow("LinearSystemDS::fillDSWithDSXML - The DSXML object doesn't exists");

  OUT("LinearSystemDS::fillDSWithDSXML\n");
}

void LinearSystemDS::display() const
{
  cout << "-----------------------------------------------------" << endl;
  cout << "____ data of the LinearSystemDS " << endl;
  DynamicalSystem::display();
  cout << "| A " << endl;
  this->A.display();
  cout << "| B " << endl;
  this->B.display();
  cout << "| u " << endl;
  this->u.display();
  cout << "| f " << endl;
  this->f.display();
  cout << "-----------------------------------------------------" << endl << endl;
}

void LinearSystemDS::init()
{
  IN("LinearSystemDS::init\n");
  this->setComputeFFunction("BasicPlugin.so", "computeF");
  this->setComputeUFunction("BasicPlugin.so", "computeU");
  OUT("LinearSystemDS::init\n");
}

void LinearSystemDS::saveDSToXML()
{
  IN("LinearSystemDS::saveDSToXML\n");
  DynamicalSystem::saveDSToXML();

  if (this->dsxml != NULL)
  {
    static_cast<LinearSystemDSXML*>(this->dsxml)->setA(&(this->A));
    static_cast<LinearSystemDSXML*>(this->dsxml)->setB(&(this->B));

    // u
    if (!(static_cast <LinearSystemDSXML*>(this->dsxml))->isUPlugin())
    {
      static_cast<LinearSystemDSXML*>(this->dsxml)->setUVector(&(this->u));
    }

    // f
    if (!(static_cast <LinearSystemDSXML*>(this->dsxml))->isFPlugin())
    {
      static_cast<LinearSystemDSXML*>(this->dsxml)->setFVector(&(this->f));
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

