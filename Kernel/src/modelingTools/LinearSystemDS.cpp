#include "LinearSystemDS.h"
#include "check.h"


LinearSystemDS::LinearSystemDS(): DynamicalSystem()
{
  this->DSType = LSDS;
  this->init();
}

LinearSystemDS::LinearSystemDS(DSXML* dsxml): DynamicalSystem(dsxml)
{
  this->DSType = LSDS;
  this->init();
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

/*SiconosVector*/SimpleVector* LinearSystemDS::getUPtr(void)
{
  return &this->u;
}

/*SiconosVector*/SimpleVector* LinearSystemDS::getFPtr(void)
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

void LinearSystemDS::createDynamicalSystem(DSXML * dsXML, int number, int n,
    SiconosVector* x0)//, NSDS * nsds)
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
    this->DSType = LSDS;
    this->number = number;
    this->n = n;
    *(this->x0) = *x0;
  }
}

LinearSystemDS* LinearSystemDS::convert(DynamicalSystem* ds)
{
  cout << "LinearSystemDS::convert (DynamicalSystem* ds)" << endl;
  LinearSystemDS* lsds = dynamic_cast<LinearSystemDS*>(ds);
  return lsds;
}

//$Log: LinearSystemDS.cpp,v $
//Revision 1.30  2005/02/11 17:36:01  charlety
//
//_ little "inspection of code"
//_ basic getters and setters passed inline
//_ getters functions passed const
//
//Revision 1.29  2005/02/10 10:35:19  jbarbier
//- new file regrouping all the const values of the model, modelingTools and numericalStrategy
//
//- new function in the LagrangianLinearR to get the H matrix corresponding to one of the 2 dynamical systems linked to the relation
//
//- new atribute of the OneStepNSProblem. A visibility table of the Interaction.
//
//Revision 1.28  2005/01/31 16:26:21  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.27  2004/09/21 11:49:09  jbarbier
//- correction in the XML save for a manual construction of the platform :
//    DS_Concerned of the Interaction
//    DS_Concerned of the Integrator
//
//- test updated for these changes
//
//Revision 1.26  2004/09/10 11:26:14  charlety
//
//_ Integration of the new version of the SiconosVector in the platform. the class simpleVector is used mostly to replace old SiconosVector. When a vector can be composite or simple, like the state of a dynamical system, a pointer on SiconosVector is used, and the vector is initialized simple or composite when the system is initialized.
//
//_ All the tests which worked with the previous version of the vector are OK with the new version.
//
//_ Example SICONOS and bouncingBall are OK
//
//_ some comments have still to be adapted to NewSiconosVector .
//
//_ This version of NewSiconosVector could be called 0.9. some details have to be fixed, it will be done before the end of September.
//
//Revision 1.25  2004/08/23 14:30:01  jbarbier
//- All the dynamical systems can be created in a comand program and added to a
//NSDS. The save is OK, but the creation of the boundary conditions is not yet
//finished.
//
//Revision 1.24  2004/08/13 11:26:58  jbarbier
//- function createNSDS complete
//
//- function createDynamicalSystem and createLinearSystemDS complete
//
//- function  createLagrangianDS in progress
//
//Revision 1.23  2004/08/12 11:55:14  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
//Revision 1.22  2004/07/29 14:25:37  jbarbier
