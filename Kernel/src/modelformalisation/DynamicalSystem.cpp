//$Id: DynamicalSystem.cpp,v 1.77 2005/03/14 16:05:26 jbarbier Exp $
#include "DynamicalSystem.h"

#include "LinearBC.h"
#include "NLinearBC.h"
#include "PeriodicBC.h"

#include "LinearDSIO.h"
#include "LagrangianDSIO.h"


DynamicalSystem::DynamicalSystem()
{
  this->init();
  this->DSType = NLSDS;
  this->dsxml = NULL;
}

DynamicalSystem::DynamicalSystem(DSXML *dsxml)
{
  IN("DynamicalSystem::DynamicalSystem(DSXML *dsxml)\n");
  this->init();
  this->DSType = NLSDS;
  this->dsxml = dsxml;
  OUT("DynamicalSystem::DynamicalSystem(DSXML *dsxml)\n");
}

DynamicalSystem::~DynamicalSystem()
{

  IN("DynamicalSystem::~DynamicalSystem()\n");

  if (this->x != NULL) delete x;
  if (this->x0 != NULL) delete x0;
  if (this->xFree != NULL) delete xFree;
  for (int i = 0; i < this->dsioVector.size(); i++)
    if (this->dsioVector[i] != NULL) delete this->dsioVector[i];

  OUT("DynamicalSystem::~DynamicalSystem()\n");
}
////////////////////////////


vector<DSInputOutput*> DynamicalSystem::getDSInputOutputs(void)
{
  return dsioVector;
}

DSInputOutput* DynamicalSystem::getDSInputOutput(int i)
{
  if (i < this->dsioVector.size())
  {
    return this->dsioVector[i];
  }
  //cout<<"# i = "<<i<<endl;
  RuntimeException::selfThrow("EqualityConstraint - getDSInputOutput : \'i\' is out of range");
}

void DynamicalSystem::setDSInputOutputs(vector<DSInputOutput*> dsioVect)
{
  this->dsioVector = dsioVect;
}

void DynamicalSystem::addDSInputOutput(DSInputOutput* dsio)
{
  //  DSInputOutput* dsioTmp;
  //  dsioTmp = new DSInputOutput();
  //  *dsioTmp = *dsio;
  //  this->dsioVector.push_back( dsioTmp );
  this->dsioVector.push_back(dsio);
}

////////////////////////////

void DynamicalSystem::setVectorFieldFunction(string pluginPath, string functionName)
{
  this->vectorFieldPtr = NULL;
  cShared.setFunction(&vectorFieldPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->vectorFieldFunctionName = plugin + ":" + functionName;
}

void DynamicalSystem::setComputeJacobianXFunction(string pluginPath, string functionName)
{
  this->computeJacobianXPtr = NULL;
  cShared.setFunction(&computeJacobianXPtr, pluginPath, functionName);

  string plugin;
  plugin = pluginPath.substr(0, pluginPath.length() - 3);
  this->computeJacobianXFunctionName = plugin + ":" + functionName;
}

void DynamicalSystem::vectorField(double time)
{
  if (vectorFieldPtr == NULL)
    RuntimeException::selfThrow("vectorField() is not linked to a plugin function");

  int size = x->size();
  this->vectorFieldPtr(&size, &time, &(*x)(0) , &xDot(0));
}

void DynamicalSystem::computeJacobianX(double time)
{
  if (computeJacobianXPtr == NULL)
    RuntimeException::selfThrow("computeJacobianX() is not linked to a plugin function");

  int size = x->size();
  this->computeJacobianXPtr(&size, &time, &(*x)(0), &jacobianX(0, 0));
}


void DynamicalSystem::swapInMemory(void)
{
  IN("DynamicalSystem::swapInMemory\n ");
  xMemory.swap(this->x);
  xDotMemory.swap(&this->xDot);
  rMemory.swap(&r);

  OUT("DynamicalSystem::swapInMemory\n ");
}


void DynamicalSystem::fillDSWithDSXML()
{
  IN("DynamicalSystem::fillDSWithDSXML\n");
  if (this->dsxml != NULL)
  {
    this->number = this->dsxml->getNumber();

    if (this->dsxml->hasId() == true) this->id = this->dsxml->getId();
    else cout << "Warning : Id is not defined in the XML " << endl;

    if (this->dsxml->hasN() == true) this->n = this->dsxml->getN();
    else cout << "Warning : n is not defined in the XML " << endl;

    if (this->dsxml->hasX0() == true)
      *(this->x0) = this->dsxml->getX0();
    else cout << "Warning : x0 is not defined in the XML " << endl;

    if (this->dsxml->hasX() == true)
    {
      *(this->x) = this->dsxml->getX();
    }
    else cout << "Warning : x is not defined in the XML " << endl;

    if (this->dsxml->hasXDot() == true)(this->xDot) = this->dsxml->getXDot();
    else cout << "Warning : xDot is not defined in the XML " << endl;

    if (this->dsxml->hasXMemory() == true) this->xMemory = SiconosMemory::SiconosMemory(this->dsxml->getXMemoryXML()); //this->dsxml->getXMemory();
    else cout << "Warning : xMemory is not defined in the XML " << endl;

    if (this->dsxml->hasXDotMemory() == true) this->xDotMemory = SiconosMemory::SiconosMemory(this->dsxml->getXDotMemoryXML()); //this->dsxml->getXDotMemory();
    else cout << "Warning : xDotMemory is not defined in the XML " << endl;

    if (this->dsxml->hasStepsInMemory() == true) this->stepsInMemory = this->dsxml->getStepsInMemory();
    else cout << "Warning : stepsInMemory is not defined in the XML " << endl;

    string plugin;
    // vectorField
    if (this->dsxml->hasVectorFieldPlugin() == true)
    {
      plugin = this->dsxml->getVectorFieldPlugin();
      this->setVectorFieldFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
    }
    else cout << "Warning : vectorFieldPlugin is not defined in the XML " << endl;


    // computeJacobianX
    if (this->dsxml->hasComputeJacobianXPlugin() == true)
    {
      plugin = this->dsxml->getComputeJacobianXPlugin();
      this->setComputeJacobianXFunction(this->cShared.getPluginName(plugin), this->cShared.getPluginFunctionName(plugin));
    }
    else cout << "Warning : computeJacobianXPlugin is not defined in the XML " << endl;

    this->r = SimpleVector(this->n);
  }
  else RuntimeException::selfThrow("DynamicalSystem::fillDSWithDSXML - DSXML object not exists");
  OUT("DynamicalSystem::fillDSWithDSXML\n");
}

void DynamicalSystem::display() const
{
  IN("DynamicalSystem::display\n");
  cout << "____ data of the Dynamical System " << endl;
  cout << "| number : " << this->number << endl;
  cout << "| id : " << this->id << endl;
  cout << "| n : " << this->n << endl;
  cout << "| x " << endl;
  this->x->display();
  cout << "| x0 " << endl;
  this->x0->display();
  cout << "| xFree " << endl;
  this->xFree->display();
  cout << "| xDot " << endl;
  this->xDot.display();
  cout << "| stepsInMemory : " << this->stepsInMemory << endl;
  cout << "| r " << endl;
  this->r.display();
  OUT("DynamicalSystem::display\n");

}

void DynamicalSystem::linkDSXML()
{
  IN("DynamicalSystem::linkDSXML\n");
  if (this->dsxml->getBoundaryConditionXML() != NULL)
  {
    //cout<<"#DynamicalSystem::linkDSXML - BC type == "<< this->dsxml->getBoundaryConditionXML()->getType() <<endl;
    if (this->dsxml->getBoundaryConditionXML()->getType() == LINEARBC_TAG)
    {
      // creation of the LinearBC with this constructor and call of a method to fill
      this->BC = new LinearBC();
      static_cast<LinearBC*>(this->BC)->createBoundaryCondition(this->dsxml->getBoundaryConditionXML());
    }

    else if (this->dsxml->getBoundaryConditionXML()->getType() == NON_LINEARBC_TAG)
    {
      // creation of the NLinearBC with this constructor and call of a method to fill
      this->BC = new NLinearBC();
      static_cast<NLinearBC*>(this->BC)->createBoundaryCondition(this->dsxml->getBoundaryConditionXML());
    }

    else if (this->dsxml->getBoundaryConditionXML()->getType() == PERIODICBC_TAG)
    {
      // creation of the PeriodicBC with this constructor and call of a method to fill
      this->BC = new PeriodicBC();
      static_cast<PeriodicBC*>(this->BC)->createBoundaryCondition(this->dsxml->getBoundaryConditionXML());
    }
    else RuntimeException::selfThrow("DynamicalSystem::linkDSXML - bad kind of BoundaryCondition : " + this->dsxml->getBoundaryConditionXML()->getType());
  }
  else this->BC = NULL;

  vector<int> nbDSIOtab = this->dsxml->getDSInputOutputNumbers();
  for (int i = 0; i < nbDSIOtab.size(); i++)
  {
    DSInputOutput *dsio;
    if (this->dsxml->getDSInputOutputXML(nbDSIOtab[i])->getType() == LINEAR_DSIO_TAG)
    {
      dsio = new LinearDSIO();
      this->dsioVector.push_back(dsio);
      static_cast<LinearDSIO*>(dsio)->createDSInputOutput(this->dsxml->getDSInputOutputXML(nbDSIOtab[i]));
    }

    else if (this->dsxml->getDSInputOutputXML(nbDSIOtab[i])->getType() == NON_LINEAR_DSIO_TAG)
    {
      dsio = new DSInputOutput();
      this->dsioVector.push_back(dsio);
      static_cast<DSInputOutput*>(dsio)->createDSInputOutput(this->dsxml->getDSInputOutputXML(nbDSIOtab[i]));
    }

    else if (this->dsxml->getDSInputOutputXML(nbDSIOtab[i])->getType() == LAGRANGIAN_DSIO_TAG)
    {
      dsio = new LagrangianDSIO();
      this->dsioVector.push_back(dsio);
      static_cast<LagrangianDSIO*>(dsio)->createDSInputOutput(this->dsxml->getDSInputOutputXML(nbDSIOtab[i]));
    }
    else RuntimeException::selfThrow("DynamicalSystem::linkDSXML - bad kind of BoundaryCondition : " + this->dsxml->getBoundaryConditionXML()->getType());
  }
  OUT("DynamicalSystem::linkDSXML\n");
}


void DynamicalSystem::init()
{
  IN("DynamicalSystem::init\n");
  //this->nsds = NULL;
  this->number = 0;
  this->id = "none";
  this->n = 0;

  this->x0 = new SimpleVector();
  this->x = new SimpleVector();
  this->xDot = SimpleVector::SimpleVector();
  this->xFree = new SimpleVector();

  this->r = SimpleVector::SimpleVector();
  this->BC = NULL;

  this->jacobianX = SiconosMatrix::SiconosMatrix();

  this->stepsInMemory = 1;
  this->setVectorFieldFunction("BasicPlugin.so", "vectorField");
  this->setComputeJacobianXFunction("BasicPlugin.so", "computeJacobianX");
  this->dsxml = NULL;

  OUT("DynamicalSystem::init\n");
}


void DynamicalSystem::initMemory(int steps)
{
  IN("DynamicalSystem::initMemory\n");
  if (steps < 0)
    RuntimeException::selfThrow("DynamicalSystem::initMemory(int steps) - steps < 0");
  else
  {
    this->stepsInMemory = steps;

    /*
     ** we made the initialization of the memories
     *
     * for rMemory, we don't need to load data for the DOM tree because there are no data saved in the XML for r
     *
     * the other memories are resized with the first parameter 'steps', and data are reloaded from the DOM tree
     * only if there are data in the DOM tree
     */

    this->rMemory = SiconosMemory::SiconosMemory(steps);
    this->xMemory = SiconosMemory::SiconosMemory(steps, this->xMemory.getSiconosMemoryXML());
    this->xDotMemory = SiconosMemory::SiconosMemory(steps, this->xDotMemory.getSiconosMemoryXML());
  }

  OUT("DynamicalSystem::initMemory\n");
}


void DynamicalSystem::saveDSToXML()
{
  IN("DynamicalSystem::saveDSToXML\n");

  /*
   * save of the BoundariesConditions
   */
  if (this->BC != NULL)
  {
    if (this->BC->getType() == LINEARBC)
      (static_cast<LinearBC*>(this->BC))->saveBCToXML();
    else if (this->BC->getType() == NLINEARBC)
      (static_cast<NLinearBC*>(this->BC))->saveBCToXML();
    else if (this->BC->getType() == PERIODICBC)
      (static_cast<PeriodicBC*>(this->BC))->saveBCToXML();
    else RuntimeException::selfThrow("DynamicalSystem::saveDSToXML - bad kind of BoundaryCondition");
  }

  if (this->dsioVector.size() != 0)
  {
    for (int i = 0; i < this->dsioVector.size(); i++)
    {
      if (this->dsioVector[i]->getType() == LINEARDSIO)
        (static_cast<LinearDSIO*>(this->dsioVector[i]))->saveDSInputOutputToXML();
      else if (this->dsioVector[i]->getType() == NLINEARDSIO)
        (static_cast<DSInputOutput*>(this->dsioVector[i]))->saveDSInputOutputToXML();
      else if (this->dsioVector[i]->getType() == LAGRANGIANDSIO)
        (static_cast<LagrangianDSIO*>(this->dsioVector[i]))->saveDSInputOutputToXML();
      else RuntimeException::selfThrow("DynamicalSystem::saveDSToXML - bad kind of DSInputOuput");
    }
  }

  if (this->dsxml != NULL)
  {
    this->dsxml->setId(this->id);
    if ((this->dsxml->getType() != LNLDS)
        && (this->dsxml->getType() != LTIDS))
      this->dsxml->setN(this->n);

    this->dsxml->setX0(this->x0);

    this->dsxml->setX(this->x);
    this->dsxml->setXMemory(&(this->xMemory));

    this->dsxml->setXDot(&this->xDot);
    this->dsxml->setXDotMemory(&(this->xDotMemory));

    this->dsxml->setStepsInMemory(this->stepsInMemory);

    this->dsxml->setR(&(this->r));

    /*
     * vectorField and computeJacobianX function must be saved only for NonLinearSystemDS
     */
    if (this->DSType == NLSDS)
    {
      this->dsxml->setVectorFieldPlugin(this->vectorFieldFunctionName);
      this->dsxml->setComputeJacobianXPlugin(this->computeJacobianXFunctionName);
    }
  }
  else RuntimeException::selfThrow("DynamicalSystem::saveDSToXML - The DSXML object doesn't exists");
  OUT("DynamicalSystem::saveDSToXML\n");
}

void DynamicalSystem::createDynamicalSystem(DSXML * dsXML, int number, int n,
    SiconosVector* x0, string vectorFieldPlugin)//, NonSmoothDynamicalSystem * nsds, BoundaryCondition* bc)
{
  IN("DynamicalSystem::createDynamicalSystem\n");
  if (dsXML != NULL)
  {
    this->DSType = NLSDS;
    //this->init();
    this->dsxml = dsXML;

    this->fillDSWithDSXML();
    this->linkDSXML();
  }
  else
  {
    this->DSType = NLSDS;
    this->number = number;
    this->n = n;

    this->x0 = new SimpleVector(n);
    this->x = new SimpleVector(n);
    this->xDot = /*new*/ SimpleVector::SimpleVector(n);
    this->xFree = new SimpleVector(n);

    *(this->x0) = *x0;
    this->setVectorFieldFunction(this->cShared.getPluginName(vectorFieldPlugin), this->cShared.getPluginFunctionName(vectorFieldPlugin));

  }
  OUT("DynamicalSystem::createDynamicalSystem\n");
}

BoundaryCondition* DynamicalSystem::createPeriodicBC()
{
  this->BC = new PeriodicBC();
  static_cast<PeriodicBC*>(this->BC)->createBoundaryCondition(NULL);
  return this->BC;
}

BoundaryCondition* DynamicalSystem::createLinearBC(SiconosVector* omega, SiconosMatrix* omega0, SiconosMatrix* omegaT)
{
  this->BC = new LinearBC();
  static_cast<LinearBC*>(this->BC)->createBoundaryCondition(NULL, omega, omega0, omegaT);
  return this->BC;
}

BoundaryCondition* DynamicalSystem::createNLinearBC()
{
  this->BC = new NLinearBC();
  static_cast<NLinearBC*>(this->BC)->createBoundaryCondition(NULL);
  return this->BC;
}

//$Log: DynamicalSystem.cpp,v $
//Revision 1.77  2005/03/14 16:05:26  jbarbier
//- manual creation of DSInputOutput saving OK
//
//- in progress for EqualityConstraint
//
//Revision 1.76  2005/03/11 15:06:20  jbarbier
//- save to XML methods of EqualityConstraint and DSInputOutput added
//
//- XML loading process modified : Model loads NSDS, then NSDS loads the DynamicalSystems, EqualityConstraints, Interactions; Modle loads Strategy, then Strategy loads TimeDiscretisation, then the Integrators, then the OneStepNSProblem
//
//Revision 1.75  2005/03/10 12:55:19  jbarbier
//- implmentation of the EqualityConstraint and DSInputOutput classes in progress
//    attributes H (DSIO) et G (EC) added in XML and managed in XML objects
//
//Revision 1.74  2005/03/09 15:30:24  jbarbier
//- add of LagrangianEC class
//
//- in progress : implementation of the EqualityConstraint and DSInputOutput - create methods
//
//Revision 1.73  2005/03/08 12:41:36  jbarbier
//- constant variables files modified :
//Some constants added in SiconosConst
//
//all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//
//Revision 1.72  2005/03/07 13:17:19  jbarbier
//- new test : Ball2D, with a ball moving in a 2D system
//
//- another constant variables moved/refactored in XMLTagsName
//- making uniform the name of the constant variables
//
//Revision 1.71  2005/02/11 17:35:54  charlety
//
//_ little "inspection of code"
//_ basic getters and setters passed inline
//_ getters functions passed const
//
//Revision 1.70  2005/02/01 11:08:41  charlety
//
//_ some displays of values during computations suppressed.
//
//Revision 1.69  2005/01/25 14:51:46  jbarbier
//- attributes id, type and XML object added to EqualityConstraint
//
//Revision 1.68  2005/01/25 13:56:02  jbarbier
//- link DynamicalSystem-DSInputOutput, NonSmoothDynamicalSystem-EqualityConstraint, EquaityConstraint-DSInputOutput and Relation-DSInputOutput available
//
//Revision 1.67  2005/01/20 14:44:48  jbarbier
//- NSDS class renamed NonSmoothDynamicalSystem
//
//- code reduce, some comments remove
//
//Revision 1.66  2005/01/11 17:08:30  jbarbier
//- last modification about the BoundaryCondition
//<BoundaryCondition>
//  <[type]>
//  <[type]/>
//<BoundaryCondition/>
//
//- modification of the xml files for this modification
//
//- version 1.2 of the xml schema
//
//Revision 1.65  2004/09/21 11:49:09  jbarbier
//- correction in the XML save for a manual construction of the platform :
//    DS_Concerned of the Interaction
//    DS_Concerned of the Integrator
//
//- test updated for these changes
//
//Revision 1.64  2004/09/16 11:35:24  jbarbier
//- save of the TimeDiscretisation in a XML file in manual creation of the
//platform which was forgotten is now available.
//
//- the save of the platform's data can be done when the platform is created with
//an XML input file and completed with dynmical systems, interactions, one-step
//non smooth problem and one-step integrator.
//
//Revision 1.63  2004/09/14 13:24:53  charlety
//
//_ changes in the interface of SiconosVector
//
//Revision 1.62  2004/09/10 11:26:07  charlety
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
//Revision 1.61  2004/09/03 14:41:41  jbarbier
//- new functions to create the boundary condition of the dynamical systems
//- new functions to add an interaction to a NonSmoothDynamicalSystem
//- new functions to create the relation and the non-smooth law of an interaciton
//
//Revision 1.60  2004/08/23 14:30:01  jbarbier
//- All the dynamical systems can be created in a comand program and added to a
//NonSmoothDynamicalSystem. The save is OK, but the creation of the boundary conditions is not yet
//finished.
//
//Revision 1.59  2004/08/20 15:26:43  jbarbier
//- creation of a Model and save in the XML is ok
//- creation of a NonSmoothDynamicalSystem and save in the XML is ok
//- creation of a NonLinearSystemDS and save in the XML is OK
//
//Revision 1.58  2004/08/18 14:37:16  jbarbier
//- creation of Model, NonSmoothDynamicalSystem, Strategy(TimeStepping and EventDriven) and
//DynamicalSystem available when the creation is in a command program
//
//Revision 1.57  2004/08/13 11:26:58  jbarbier
//- function createNSDS complete
//
//- function createDynamicalSystem and createLinearSystemDS complete
//
//- function  createLagrangianNLDS in progress
//
//Revision 1.56  2004/08/10 12:04:27  jbarbier
//- save of the plugin's name for fInt
//
//Revision 1.55  2004/08/09 15:00:52  jbarbier
//- changes in the cardinality of some attributes of the DynamicalSystem,
//OneStepIntegrator
//
//- modifications in classes Moreau, Lsodar, Adams for these new cardinalities
//
//- corrections in the test xml files
//
//Revision 1.54  2004/08/06 11:27:53  jbarbier
//- new tests with the XML and the optional attributes
//
//- tests on the save of the XML data
//
//Revision 1.53  2004/08/05 12:44:40  jbarbier
//- loading XML file with no OneStepNSProblem succesfull
//
//- NonLinearSystemDS is now available
//
//Revision 1.52  2004/08/04 11:03:22  jbarbier
//- about the SiconosMemory : when a SiconosMemory has a maxSize greater than the
//number of steps in memory required by an integrator, the oldest SiconosVector
//are deleted
//
//- the way to initialize the SiconosMemory by the integrator has been updated to
//match with these changes
//
//Revision 1.51  2004/08/03 12:07:11  jbarbier
//- all test on th eModel are successfull
//
//- new tests on the Model with the opening of XML file
//
//- link TimeDiscretisation -> Strategy
//
//- attribute T of the Model is now optional
//
//Revision 1.50  2004/08/02 09:26:25  jbarbier
//- xml save for SiconosMemory corrected
//- temporary operation in Moreau::integrate because of the current version of
//SiconosVector
//
//Revision 1.49  2004/07/30 14:37:14  jbarbier
//- saving methods for DynamicalSystemXML and LagrangianNLDSXML
//
//Revision 1.48  2004/07/28 14:13:48  charlety
//
//_ add of function to get SiconosMemory objects from XML (xMemory, xDotMemory in DynamicalSystem, etc.)
//
//Revision 1.47  2004/07/23 14:39:25  jbarbier
//- createModel, createNSDS, createDynamicalSystem, createBoundaryCondition OK
//it's the first step, it do the same thing that before, but the process is
//unified and it must simply add the traitment for the creation of the nodes in
//the DOM tree
//
//Revision 1.46  2004/07/22 15:16:22  jbarbier
//- xml_test.xml corrected with because of n which is equals to 2*ndof
//- manage the save of n (optional attribute) for Dynamical Systems
//
//Revision 1.45  2004/07/12 14:50:54  jbarbier
//- integration of the Memory in the XML (vol.1)
//
//Revision 1.44  2004/07/12 13:04:33  jbarbier
//- $id$, $log$, $date$, ... added in the XML management files
//- n id calculated with ndof for Lagrangian dynamical systems
//
//Revision 1.43  2004/07/12 11:22:53  charlety
//
//_ the state vector x of the dynamical system is now plugged to q and the velocity when this system is a Lagrangian one.
//
//_ A function to force a SiconosVector to be composite has been written. This is a temporary solution. We should use the operator = and the constructor by copy instead. This problem will be fixed later in the summer.
//
//Revision 1.41  2004/07/09 14:18:55  jbarbier
//-management of the optional attributes of the DynamicalSystem
//-new node t for current time in the time of the NonSmoothDynamicalSystem
//
//Revision 1.40  2004/07/09 11:14:53  charlety
//
//_ Added a constructor by copy and an operator = in class SiconosMemory
//_ the getters on memory in DynamicalSystems return now some pointers
//
//Revision 1.39  2004/07/05 12:38:08  charlety
//
//try of false plugin developed in LagrangianTIDS. The Moreau integrator considers it like a LagrangianNLDS, but this is not the plugin which is used to compute the external strength, but a function of the LagrangianTIDS.
//
//Revision 1.38  2004/07/02 14:48:28  acary
//Added MACRO IN and OUT
//
//Revision 1.37  2004/07/02 14:40:20  jbarbier
//some IN and OUT added
//BoundaryConditon saveToXML ok
//problem after, during the call of the destructors
//
//Revision 1.36  2004/06/30 13:08:21  jbarbier
//DoxygenSiconos.cfg moved into config/
//ICDLLSharedLibrary renamed SiconosSharedLibrary
//
//Revision 1.35  2004/06/29 15:05:53  acary
//Change  in the Display method of the Dynamical system
//
//Revision 1.34  2004/06/29 08:49:57  acary
//Ajout des commentaires Doxygen et des tages CVS
//
