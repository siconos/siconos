//$Id: NonSmoothDynamicalSystem.cpp,v 1.9 2005/03/14 16:05:27 jbarbier Exp $
#include "NonSmoothDynamicalSystem.h"

#include "LagrangianNLDS.h"
#include "LagrangianTIDS.h"
#include "LinearSystemDS.h"
//#include "LagrangianNLDSXML.h"
//#include "LagrangianTIDSXML.h"
//#include "LinearSystemDSXML.h"
#include "LinearEC.h"
#include "LinearTIEC.h"
#include "LagrangianEC.h"


#include "check.h"

NonSmoothDynamicalSystem::NonSmoothDynamicalSystem()
{
  this->nsdsxml = NULL;
  this->BVP = false;
}

NonSmoothDynamicalSystem::NonSmoothDynamicalSystem(bool bvp)
{
  this->nsdsxml = NULL;
  this->BVP = bvp;
}

NonSmoothDynamicalSystem::NonSmoothDynamicalSystem(NonSmoothDynamicalSystem* nsds)
{
  this->BVP = nsds->isBVP();
  this->DSVector = nsds->getDynamicalSystems();
  this->interactionVector = nsds->getInteractions();
  this->nsdsxml = nsds->getNSDSXML();
}

NonSmoothDynamicalSystem::NonSmoothDynamicalSystem(NSDSXML* nsdsxml)
{
  this->BVP = false;
  this->nsdsxml = nsdsxml;
}

NonSmoothDynamicalSystem::NonSmoothDynamicalSystem(string type)
{
  if (type == "BVP")
  {
    this->BVP = true;
  }
  else if (type == "IVP")
  {
    this->BVP = false;
  }
  else
  {
    // raise exception
    // or make IVP or BVP as default choice
  }
}

NonSmoothDynamicalSystem::~NonSmoothDynamicalSystem()
{
  IN("NonSmoothDynamicalSystem::~NonSmoothDynamicalSystem\n");
  int i;

  if (this->DSVector.size() > 0)
  {
    for (i = 0; i < this->DSVector.size(); i++)
    {
      delete this->DSVector[i];
    }
    this->DSVector.clear();
  }

  if (this->interactionVector.size() > 0)
  {
    for (i = 0; i < this->interactionVector.size(); i++)
    {
      delete this->interactionVector[i];
    }
    this->interactionVector.clear();
  }
  OUT("NonSmoothDynamicalSystem::~NonSmoothDynamicalSystem\n");
}


DynamicalSystem* NonSmoothDynamicalSystem::getDynamicalSystem(int nb)
{
  IN("DynamicalSystem* NonSmoothDynamicalSystem::getDynamicalSystem(int nb)\n");
  if (nb < this->DSVector.size())
  {
    OUT("DynamicalSystem* NonSmoothDynamicalSystem::getDynamicalSystem(int nb)\n");
    return this->DSVector[nb];
  }
  RuntimeException::selfThrow("NonSmoothDynamicalSystem - getDynamicalSystem : \'nb\' is out of range");
}

DynamicalSystem* NonSmoothDynamicalSystem::getDynamicalSystemOnNumber(int nb)
{
  for (unsigned int i = 0; i < this->DSVector.size(); i++)
  {
    if (DSVector[i]->getNumber() == nb)
    {
      return DSVector[i];
    }
  }
  RuntimeException::selfThrow("NonSmoothDynamicalSystem::getDynamicalSystemOnNumber : DSVector[i] == NULL");
}

bool NonSmoothDynamicalSystem::hasDynamicalSystemNumber(int nb)
{
  for (unsigned int i = 0; i < this->DSVector.size(); i++)
  {
    if (DSVector[i]->getNumber() == nb)
    {
      return true;
    }
  }
  return false;
}


Interaction* NonSmoothDynamicalSystem::getInteraction(int nb)
{
  if (nb < this->interactionVector.size())
  {
    return this->interactionVector[nb];
  }
  RuntimeException::selfThrow("NonSmoothDynamicalSystem - getInteraction : \'nb\' is out of range");
}

Interaction* NonSmoothDynamicalSystem::getInteractionOnNumber(int nb)
{
  for (unsigned int i = 0; i < this->interactionVector.size(); i++)
  {
    if (interactionVector[i] != NULL)
      if (interactionVector[i]->getNumber() == nb)
      {
        return interactionVector[i];
      }
  }
  RuntimeException::selfThrow("NonSmoothDynamicalSystem::getInteractionOnNumber : interactionVector[i] == NULL");
  return NULL;
}

bool NonSmoothDynamicalSystem::hasInteractionNumber(int nb)
{
  for (unsigned int i = 0; i < this->interactionVector.size(); i++)
  {
    if (interactionVector[i] != NULL)
      if (interactionVector[i]->getNumber() == nb)
      {
        return true;
      }
  }
  return false;
}


void NonSmoothDynamicalSystem::addDynamicalSystem(DynamicalSystem *ds)//, BoundaryCondition* bc)
{
  IN("NonSmoothDynamicalSystem::addDynamicalSystem\n");

  ds->setNSDS(this);
  this->DSVector.push_back(ds);

  //  DynamicalSystem* dsTmp;
  //
  //  if( ds->getType() == LNLDS )
  //  {
  //    dsTmp = new LagrangianNLDS();
  //    *dsTmp = *ds;
  //    this->DSVector.push_back( dsTmp );
  //  }
  //  else if( ds->getType() == LTIDS )
  //  {
  //    dsTmp = new LagrangianTIDS();
  //    *dsTmp = *ds;
  //    this->DSVector.push_back( dsTmp );
  //    this->DSVector.push_back( ds );
  //  }
  //  else if( ds->getType() == LSDS )
  //  {
  //    dsTmp = new LinearSystemDS();
  //    *dsTmp = *ds;
  //    this->DSVector.push_back( dsTmp );
  //  }
  //  else if( ds->getType() == NLSDS )
  //  {
  //    dsTmp = new DynamicalSystem();
  //    *dsTmp = *ds;
  //    this->DSVector.push_back( dsTmp );
  //  }
  //  else RuntimeException::selfThrow("NonSmoothDynamicalSystem::addDynamicalSystem - bad kind of DynamicalSystem");

  OUT("NonSmoothDynamicalSystem::addDynamicalSystem\n");
}


void NonSmoothDynamicalSystem::addInteraction(Interaction *inter)
{
  Interaction* interTmp;
  interTmp = new Interaction();
  *interTmp = *inter;
  this->interactionVector.push_back(interTmp);
}


EqualityConstraint* NonSmoothDynamicalSystem::getEqualityConstraint(int i)
{
  if (i < this->ecVector.size())
  {
    return this->ecVector[i];
  }
  RuntimeException::selfThrow("NonSmoothDynamicalSystem - getEqualityConstraint : \'i\' is out of range");
}


void NonSmoothDynamicalSystem::addEqualityConstraint(EqualityConstraint* ec)
{
  EqualityConstraint* ecTmp;
  ecTmp = new EqualityConstraint();
  *ecTmp = *ec;
  this->ecVector.push_back(ecTmp);
}

////////////////////////

void NonSmoothDynamicalSystem::linkNSDSXML()
{
  IN("NonSmoothDynamicalSystem::linkNSDSXML\n");

  int i = 0;

  // get all the DSXML objects then create the DS for this DSXML and add this DS to the vector of DS of the NonSmoothDynamicalSystem
  vector<int> nbDStab = this->nsdsxml->getDSNumbers();
  for (i = 0; i < nbDStab.size(); i++)
  {
    DynamicalSystem *ds;
    // with the data of the XML object, we know the type of DynamicalSystem, so we can instanciate the right type of DynamicalSystem

    // DynamicalSystem - LagrangianNLDS
    if ((this->nsdsxml->getDSXML(nbDStab[i]))->getType() == LAGRANGIAN_NON_LINEARDS_TAG)
    {
      // creation of the LagrangianNLDS with this constructor and call of a method to fill
      ds = new LagrangianNLDS();
      this->DSVector.push_back(ds);
      (static_cast<LagrangianNLDS*>(ds))->createDynamicalSystem(this->nsdsxml->getDSXML(nbDStab[i]));
      ds->setNSDS(this);
    }
    else if ((this->nsdsxml->getDSXML(nbDStab[i]))->getType() == LAGRANGIAN_TIME_INVARIANTDS_TAG)
    {
      ds = new LagrangianTIDS();
      this->DSVector.push_back(ds);
      (static_cast<LagrangianTIDS*>(ds))->createDynamicalSystem(this->nsdsxml->getDSXML(nbDStab[i]));
      ds->setNSDS(this);
    }
    else if ((this->nsdsxml->getDSXML(nbDStab[i]))->getType() == LINEAR_SYSTEMDS_TAG)
    {
      ds = new LinearSystemDS();
      this->DSVector.push_back(ds);
      (static_cast<LinearSystemDS*>(ds))->createDynamicalSystem(this->nsdsxml->getDSXML(nbDStab[i]));
      ds->setNSDS(this);
    }
    else if ((this->nsdsxml->getDSXML(nbDStab[i]))->getType() == NON_LINEAR_SYSTEMDS_TAG)
    {
      ds = new DynamicalSystem();
      this->DSVector.push_back(ds);
      ds->createDynamicalSystem(this->nsdsxml->getDSXML(nbDStab[i]));
      ds->setNSDS(this);
    }
    else RuntimeException::selfThrow("NonSmoothDynamicalSystem::LinkNSDSXML - bad kind of DynamicalSystem");

    // for the other DS, we must have the xxxDSXML.h .cpp objects needed, and the objects in the DynamicalSystems inherited of the platform
    /*
     *  other "if" to create LagrangianTIDS and LinearSystemDS
     *
     */
  }

  // get all the InteractionXML objects then create the Interaction for this InteractionXML and add this Interaction to the vector of Interaction of the NonSmoothDynamicalSystem
  vector<int> nbInteractionTab = this->nsdsxml->getInteractionNumbers();

  int ds1, ds2;
  for (i = 0; i < nbInteractionTab.size(); i++)
  {
    // the creation of the Interaction with this constructor call a method to fill
    Interaction * inter;
    inter = new Interaction();
    inter->createInteraction(this->nsdsxml->getInteractionXML(nbInteractionTab[i]));
    cout << "#@# createInteraction of the NonSmoothDynamicalSystem::linkNSDSXML ....... in progress" << endl;
    cout << " size of the vector of DS of the NonSmoothDynamicalSystem == " << this->DSVector.size() << endl;
    ds1 = (this->nsdsxml->getInteractionXML(nbInteractionTab[i])->getDSConcerned()[0])[0];
    ds2 = (this->nsdsxml->getInteractionXML(nbInteractionTab[i])->getDSConcerned()[0])[1];
    cout << " id of the ds1 == " << ds1 << endl;
    cout << " id of the ds2 == " << ds2 << endl;
    inter->setDynamicalSystems(this->getDynamicalSystemOnNumber(ds1), this->getDynamicalSystemOnNumber(ds2));
    this->interactionVector.push_back(inter);
  }
  if (nbInteractionTab.size() == 0) cout << "Warning : There are no Interaction defined." << endl;


  // get all the EqualityConstraintXML objects then create the EqualityConstraint for this EqualityConstraintXML
  vector<int> nbECtab = this->nsdsxml->getEqualityConstraintNumbers();
  for (i = 0; i < nbECtab.size(); i++)
  {
    EqualityConstraint *ec;
    if ((this->nsdsxml->getEqualityConstraintXML(nbECtab[i]))->getType() == LINEAR_EC_TAG)
    {
      ec = new LinearEC();
      this->ecVector.push_back(ec);
      (static_cast<LinearEC*>(ec))->createEqualityConstraint(this->nsdsxml->getEqualityConstraintXML(nbECtab[i]));
    }
    else if ((this->nsdsxml->getEqualityConstraintXML(nbECtab[i]))->getType() == NON_LINEAR_EC_TAG)
    {
      ec = new EqualityConstraint();
      this->ecVector.push_back(ec);
      (static_cast<EqualityConstraint*>(ec))->createEqualityConstraint(this->nsdsxml->getEqualityConstraintXML(nbECtab[i]));
    }
    else if ((this->nsdsxml->getEqualityConstraintXML(nbECtab[i]))->getType() == LINEAR_TIME_INVARIANT_EC_TAG)
    {
      ec = new LinearTIEC();
      this->ecVector.push_back(ec);
      (static_cast<LinearTIEC*>(ec))->createEqualityConstraint(this->nsdsxml->getEqualityConstraintXML(nbECtab[i]));
    }
    else if ((this->nsdsxml->getEqualityConstraintXML(nbECtab[i]))->getType() == LAGRANGIAN_EC_TAG)
    {
      ec = new LagrangianEC();
      this->ecVector.push_back(ec);
      (static_cast<LagrangianEC*>(ec))->createEqualityConstraint(this->nsdsxml->getEqualityConstraintXML(nbECtab[i]));
    }
    else RuntimeException::selfThrow("NonSmoothDynamicalSystem::LinkNSDSXML - bad kind of EqualityConstraint");
  }

  OUT("NonSmoothDynamicalSystem::linkNSDSXML\n");

}

void NonSmoothDynamicalSystem::fillNSDSWithNSDSXML()
{
  IN("NonSmoothDynamicalSystem::fillNSDSWithNSDSXML\n");

  if (this->nsdsxml != NULL)
  {
    //this->BVP = this->nsdsxml->isBVP();
  }
  else RuntimeException::selfThrow("NonSmoothDynamicalSystem::fillNSDSWithNSDSXML - The NSDSXML object doesn't exists");

  OUT("NonSmoothDynamicalSystem::fillNSDSWithNSDSXML\n");

}

void NonSmoothDynamicalSystem::saveNSDSToXML()
{
  IN("NonSmoothDynamicalSystem::saveNSDSToXML\n");
  int size, i;

  if (this->nsdsxml != NULL)
  {
    this->nsdsxml->setBVP(this->BVP);// no need to change the value of BVP, it mustn't change anyway

    size = this->DSVector.size();
    for (i = 0; i < size; i++)
    {
      if (this->DSVector[i]->getType() == LNLDS)
        (static_cast<LagrangianNLDS*>(this->DSVector[i]))->saveDSToXML();
      else if (this->DSVector[i]->getType() == LTIDS)
        (static_cast<LagrangianTIDS*>(this->DSVector[i]))->saveDSToXML();
      else if (this->DSVector[i]->getType() == LSDS)
        (static_cast<LinearSystemDS*>(this->DSVector[i]))->saveDSToXML();
      else if (this->DSVector[i]->getType() == NLSDS)
        this->DSVector[i]->saveDSToXML();
      else RuntimeException::selfThrow("NonSmoothDynamicalSystem::saveToXML - bad kind of DynamicalSystem");
    }

    size = this->ecVector.size();
    for (i = 0; i < size; i++)
    {
      if (this->ecVector[i]->getType() == LINEAREC)
        (static_cast<LinearEC*>(this->ecVector[i]))->saveEqualityConstraintToXML();
      else if (this->ecVector[i]->getType() == LINEARTIEC)
        (static_cast<LinearTIEC*>(this->ecVector[i]))->saveEqualityConstraintToXML();
      else if (this->ecVector[i]->getType() == NLINEAREC)
        this->ecVector[i]->saveEqualityConstraintToXML();
      else if (this->ecVector[i]->getType() == LAGRANGIANEC)
        (static_cast<LagrangianEC*>(this->ecVector[i]))->saveEqualityConstraintToXML();
      else RuntimeException::selfThrow("NonSmoothDynamicalSystem::saveToXML - bad kind of EqualityConstraint");
    }

    size = this->interactionVector.size();
    for (i = 0; i < size; i++)
      this->interactionVector[i]->saveInteractionToXML();

  }
  else RuntimeException::selfThrow("NonSmoothDynamicalSystem::saveNSDSToXML - The NSDSXML object doesn't exists");

  OUT("NonSmoothDynamicalSystem::saveNSDSToXML\n");
}

void NonSmoothDynamicalSystem::createNonSmoothDynamicalSystem(NSDSXML * nsdsXML, bool bvp)//, Model * model)//, xmlNode * rootDOMTreeNode)
{
  if (nsdsXML != NULL)
  {
    /*
     * if bvp is true, it is not used, data of the DOM tree are more important
     */
    this->nsdsxml = nsdsXML;
    this->fillNSDSWithNSDSXML();
    cout << "NonSmoothDynamicalSystem filled" << endl;
    this->linkNSDSXML();
    cout << "NonSmoothDynamicalSystem linked" << endl;
  }
  else
  {
    this->BVP = bvp;
  }
}

void NonSmoothDynamicalSystem::display() const
{
  IN("NonSmoothDynamicalSystem::display\n");

  cout << "| this = " << this << endl;
  cout << "| BVP = " << this->BVP << endl;
  cout << "| &nsdsxml = " << this->nsdsxml << endl;
  cout << "|===========================" << endl;

  OUT("NonSmoothDynamicalSystem::display\n");
}

DynamicalSystem* NonSmoothDynamicalSystem::addNonLinearSystemDS(int number, int n,
    SiconosVector* x0, string vectorFieldPlugin)
{
  if (!this->hasDynamicalSystemNumber(number))
  {
    DynamicalSystem* dsTmp;

    dsTmp = new DynamicalSystem();
    dsTmp->createDynamicalSystem(NULL, number, n, x0, vectorFieldPlugin);
    this->DSVector.push_back(dsTmp);
    dsTmp->setNSDS(this);

    return dsTmp;
  }
  else
  {
    char n[50];
    sprintf(n, "%i", number);
    string num = n;
    string msg = "NonSmoothDynamicalSystem::addNonLinearSystemDS : ERROR - The DynamicalSystem number " + num + " is already declared!";
    RuntimeException::selfThrow(msg);
  }
}

DynamicalSystem* NonSmoothDynamicalSystem::addLinearSystemDS(int number, int n,
    SiconosVector* x0)
{
  if (!this->hasDynamicalSystemNumber(number))
  {
    DynamicalSystem* dsTmp;

    dsTmp = new LinearSystemDS();
    static_cast<LinearSystemDS*>(dsTmp)->createDynamicalSystem(NULL, number, n, x0);
    this->DSVector.push_back(dsTmp);
    dsTmp->setNSDS(this);

    return dsTmp;
  }
  else
  {
    char n[50];
    sprintf(n, "%i", number);
    string num = n;
    string msg = "NonSmoothDynamicalSystem::addLinearSystemDS : ERROR - The DynamicalSystem number " + num + " is already declared!";
    RuntimeException::selfThrow(msg);
  }
}

DynamicalSystem* NonSmoothDynamicalSystem::addLagrangianNLDS(int number, int ndof,
    SiconosVector* q0, SiconosVector* velocity0,
    string mass, string fInt, string fExt,
    string jacobianQFInt, string jacobianVelocityFInt,
    string jacobianQQNLInertia, string jacobianVelocityQNLInertia,
    string QNLInertia)
{
  if (!this->hasDynamicalSystemNumber(number))
  {
    DynamicalSystem* dsTmp;

    dsTmp = new LagrangianNLDS();
    static_cast<LagrangianNLDS*>(dsTmp)->createDynamicalSystem(NULL, number, ndof,
        q0, velocity0, mass, fInt, fExt,
        jacobianQFInt, jacobianVelocityFInt,
        jacobianQQNLInertia, jacobianVelocityQNLInertia,
        QNLInertia);
    this->DSVector.push_back(dsTmp);
    dsTmp->setNSDS(this);

    return dsTmp;
  }
  else
  {
    char n[50];
    sprintf(n, "%i", number);
    string num = n;
    string msg = "NonSmoothDynamicalSystem::addLagrangianNLDS : ERROR - The DynamicalSystem number " + num + " is already declared!";
    RuntimeException::selfThrow(msg);
  }
}

DynamicalSystem* NonSmoothDynamicalSystem::addLagrangianTIDS(int number, int ndof,
    SiconosVector* q0, SiconosVector* velocity0,
    SiconosMatrix* mass, string fExt,
    SiconosMatrix* K, SiconosMatrix* C)
{
  if (!this->hasDynamicalSystemNumber(number))
  {
    DynamicalSystem* dsTmp;

    dsTmp = new LagrangianTIDS();
    static_cast<LagrangianTIDS*>(dsTmp)->createDynamicalSystem(NULL, number, ndof,
        q0, velocity0, mass, fExt, K, C);
    this->DSVector.push_back(dsTmp);
    dsTmp->setNSDS(this);

    return dsTmp;
  }
  else
  {
    char n[50];
    sprintf(n, "%i", number);
    string num = n;
    string msg = "NonSmoothDynamicalSystem::addLagrangianTIDS : ERROR - The DynamicalSystem number " + num + " is already declared!";
    RuntimeException::selfThrow(msg);
  }
}

Interaction* NonSmoothDynamicalSystem::addInteraction(int number, int nInter, vector<int>* status, vector<DynamicalSystem*>* dsConcerned)
{
  if (!this->hasInteractionNumber(number))
  {
    Interaction* interTmp;
    interTmp = new Interaction();
    interTmp->createInteraction(NULL, number, nInter, status, dsConcerned);

    this->interactionVector.push_back(interTmp);
    return interTmp;
  }
  else
  {
    char n[50];
    sprintf(n, "%i", number);
    string num = n;
    string msg = "NonSmoothDynamicalSystem::addInteraction : ERROR - The Interaction number " + num + " is already declared!";
    RuntimeException::selfThrow(msg);
  }
}

//$Log: NonSmoothDynamicalSystem.cpp,v $
//Revision 1.9  2005/03/14 16:05:27  jbarbier
//- manual creation of DSInputOutput saving OK
//
//- in progress for EqualityConstraint
//
//Revision 1.8  2005/03/11 15:06:20  jbarbier
//- save to XML methods of EqualityConstraint and DSInputOutput added
//
//- XML loading process modified : Model loads NSDS, then NSDS loads the DynamicalSystems, EqualityConstraints, Interactions; Modle loads Strategy, then Strategy loads TimeDiscretisation, then the Integrators, then the OneStepNSProblem
//
//Revision 1.7  2005/03/10 12:55:20  jbarbier
//- implmentation of the EqualityConstraint and DSInputOutput classes in progress
//    attributes H (DSIO) et G (EC) added in XML and managed in XML objects
//
//Revision 1.6  2005/03/09 15:30:32  jbarbier
//- add of LagrangianEC class
//
//- in progress : implementation of the EqualityConstraint and DSInputOutput - create methods
//
//Revision 1.5  2005/03/04 15:35:26  jbarbier
//- README files added for some samples
//
//- beginning of the refactoring of XML module constants
//
//Revision 1.4  2005/02/11 17:36:02  charlety
//
//_ little "inspection of code"
//_ basic getters and setters passed inline
//_ getters functions passed const
//
//Revision 1.3  2005/02/01 11:08:42  charlety
//
//_ some displays of values during computations suppressed.
//
//Revision 1.2  2005/01/25 13:56:03  jbarbier
//- link DynamicalSystem-DSInputOutput, NonSmoothDynamicalSystem-EqualityConstraint, EquaityConstraint-DSInputOutput and Relation-DSInputOutput available
//
//Revision 1.1  2005/01/20 14:44:49  jbarbier
//- NSDS class renamed NonSmoothDynamicalSystem
//
//- code reduce, some comments remove
//
//Revision 1.47  2005/01/10 17:06:37  jbarbier
//- attribute "size" is now unused in the code
//
//- xml schema v1.2 is in progress
//
//Revision 1.46  2004/09/22 10:54:44  jbarbier
//- light modification according to the attribute mass of the lagrangian dynamical
//systems. The lagrangianNLDS take always an function from a plugin to compute the
//mass, whereas the lagrangianTIDS needs only a matrix.
//
//- xml input files have been modified in consequence
//
//Revision 1.45  2004/09/21 11:49:09  jbarbier
//- correction in the XML save for a manual construction of the platform :
//    DS_Concerned of the Interaction
//    DS_Concerned of the Integrator
//
//- test updated for these changes
//
//Revision 1.44  2004/09/16 11:35:24  jbarbier
//- save of the TimeDiscretisation in a XML file in manual creation of the
//platform which was forgotten is now available.
//
//- the save of the platform's data can be done when the platform is created with
//an XML input file and completed with dynmical systems, interactions, one-step
//non smooth problem and one-step integrator.
//
//Revision 1.43  2004/09/10 11:26:15  charlety
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
//Revision 1.42  2004/09/10 08:04:47  jbarbier
//- XML save available for BoundaryCondition and Interaction
//
//Revision 1.41  2004/09/03 14:41:42  jbarbier
//- new functions to create the boundary condition of the dynamical systems
//- new functions to add an interaction to a NonSmoothDynamicalSystem
//- new functions to create the relation and the non-smooth law of an interaciton
//
//Revision 1.40  2004/08/23 14:30:02  jbarbier
//- All the dynamical systems can be created in a comand program and added to a
//NonSmoothDynamicalSystem. The save is OK, but the creation of the boundary conditions is not yet
//finished.
//
//Revision 1.39  2004/08/20 15:26:45  jbarbier
//- creation of a Model and save in the XML is ok
//- creation of a NonSmoothDynamicalSystem and save in the XML is ok
//- creation of a NonLinearSystemDS and save in the XML is OK
//
//Revision 1.38  2004/08/18 14:37:18  jbarbier
//- creation of Model, NonSmoothDynamicalSystem, Strategy(TimeStepping and EventDriven) and
//DynamicalSystem available when the creation is in a command program
//
//Revision 1.37  2004/08/13 11:26:58  jbarbier
//- function createNSDS complete
//
//- function createDynamicalSystem and createLinearSystemDS complete
//
//- function  createLagrangianNLDS in progress
//
//Revision 1.36  2004/08/12 11:55:17  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
//Revision 1.35  2004/08/06 10:46:30  charlety
//
//_ example Oscillator in progress
//_ corrected a bug : theXML  save of a system without OneStepIntegrator was not OK.
//
//Revision 1.34  2004/08/05 12:44:41  jbarbier
//- loading XML file with no OneStepNSProblem succesfull
//
//- NonLinearSystemDS is now available
//
//Revision 1.33  2004/08/04 14:51:02  jbarbier
//- new test using xml_uncomplete7.xml, test with no interaction defined in the
//XML input file
//
//- for the TimeDiscretisation, the triplet (t0,T,h), (t0,T,N) or (t0,h,N) ids
//required, the missing element is now computed
//
//Revision 1.32  2004/07/28 14:40:18  jbarbier
//- tests on the platform
//
//Revision 1.31  2004/07/26 13:28:05  jbarbier
//- add of the methods : createInteraction and createRelation
//
//Revision 1.30  2004/07/23 14:39:27  jbarbier
//- createModel, createNSDS, createDynamicalSystem, createBoundaryCondition OK
//it's the first step, it do the same thing that before, but the process is
//unified and it must simply add the traitment for the creation of the nodes in
//the DOM tree
//
//Revision 1.29  2004/07/05 12:45:13  jbarbier
//optionnal Matrix and Vector attributes save ... ok
//XML save  ok even without boundarycondition
//
//Revision 1.28  2004/07/02 14:48:29  acary
//Added MACRO IN and OUT
//
