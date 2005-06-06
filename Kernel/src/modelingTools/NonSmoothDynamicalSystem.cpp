#include "NonSmoothDynamicalSystem.h"

// includes to be deleted thanks to factories
#include "LagrangianDS.h"
#include "LagrangianLinearTIDS.h"
#include "LinearSystemDS.h"
#include "LinearEC.h"
#include "LinearTIEC.h"
#include "LagrangianEC.h"
#include "LagrangianLinearEC.h"

using namespace std;

// --- CONSTRUCTORS/DESTRUCTOR ---

// Default constructor
NonSmoothDynamicalSystem::NonSmoothDynamicalSystem():
  BVP(false), nsdsxml(NULL)
{}

// copy constructor
NonSmoothDynamicalSystem::NonSmoothDynamicalSystem(const NonSmoothDynamicalSystem& nsds):
  BVP(false), nsdsxml(NULL)
{
  BVP = nsds.isBVP();
  DSVector = nsds.getDynamicalSystems();
  interactionVector = nsds.getInteractions();
  *nsdsxml = *nsds.getNSDSXMLPtr();
}

// xml constuctor
NonSmoothDynamicalSystem::NonSmoothDynamicalSystem(NSDSXML* newNsdsxml):
  BVP(false), nsdsxml(newNsdsxml)
{
  if (this->nsdsxml != NULL)
  {
    int i = 0;
    // DS Vector fill-in
    vector<int> nbDStab = nsdsxml->getDSNumbers();
    for (i = 0; i < nbDStab.size(); i++)
    {
      DynamicalSystem *ds;
      // DynamicalSystem - LagrangianDS
      if ((nsdsxml->getDSXML(nbDStab[i]))->getType() == LAGRANGIAN_NON_LINEARDS_TAG)
      {
        ds = new LagrangianDS(nsdsxml->getDSXML(nbDStab[i]));
        DSVector.push_back(ds);
        isDSVectorAllocatedIn.push_back(true);
        ds->setNSDSPtr(this);
      }
      else if ((nsdsxml->getDSXML(nbDStab[i]))->getType() == LAGRANGIAN_TIME_INVARIANTDS_TAG)
      {
        ds = new LagrangianLinearTIDS(nsdsxml->getDSXML(nbDStab[i]));
        DSVector.push_back(ds);
        isDSVectorAllocatedIn.push_back(true);
        ds->setNSDSPtr(this);
      }
      else if ((nsdsxml->getDSXML(nbDStab[i]))->getType() == LINEAR_SYSTEMDS_TAG)
      {
        ds = new LinearSystemDS(nsdsxml->getDSXML(nbDStab[i]));
        DSVector.push_back(ds);
        isDSVectorAllocatedIn.push_back(true);
        ds->setNSDSPtr(this);
      }
      else if ((nsdsxml->getDSXML(nbDStab[i]))->getType() == NON_LINEAR_SYSTEMDS_TAG)
      {
        ds = new DynamicalSystem(nsdsxml->getDSXML(nbDStab[i]));
        DSVector.push_back(ds);
        isDSVectorAllocatedIn.push_back(true);
        ds->setNSDSPtr(this);
      }
      else RuntimeException::selfThrow("NonSmoothDynamicalSystem::xml constructor, wrong Dynamical System type" + ((nsdsxml->getDSXML(nbDStab[i]))->getType()));
    }
    // Interaction vector fill-in
    vector<int> nbInteractionTab = nsdsxml->getInteractionNumbers();
    int ds1, ds2;
    for (i = 0; i < nbInteractionTab.size(); i++)
    {
      // the creation of the Interaction with this constructor call a method to fill
      Interaction * inter = new Interaction(nsdsxml->getInteractionXML(nbInteractionTab[i]));
      ds1 = (nsdsxml->getInteractionXML(nbInteractionTab[i])->getDSConcerned()[0])[0];
      ds2 = (nsdsxml->getInteractionXML(nbInteractionTab[i])->getDSConcerned()[0])[1];
      inter->setDynamicalSystems(getDynamicalSystemPtrNumber(ds1), getDynamicalSystemPtrNumber(ds2));
      interactionVector.push_back(inter);
      isInteractionVectorAllocatedIn.push_back(true);
    }
    if (nbInteractionTab.size() == 0) cout << "Warning : No Interaction defined." << endl;

    // Algebraic constraints fill-in
    // get all the EqualityConstraintXML objects then create the EqualityConstraint for this EqualityConstraintXML
    vector<int> nbECtab = nsdsxml->getEqualityConstraintNumbers();
    for (i = 0; i < nbECtab.size(); i++)
    {
      EqualityConstraint *ec;
      if ((nsdsxml->getEqualityConstraintXML(nbECtab[i]))->getType() == LINEAR_EC_TAG)
      {
        ec = new LinearEC();
        ecVector.push_back(ec);
        isEcVectorAllocatedIn.push_back(true);
        (static_cast<LinearEC*>(ec))->createEqualityConstraint(nsdsxml->getEqualityConstraintXML(nbECtab[i]));
      }
      else if ((nsdsxml->getEqualityConstraintXML(nbECtab[i]))->getType() == NON_LINEAR_EC_TAG)
      {
        ec = new EqualityConstraint();
        ecVector.push_back(ec);
        isEcVectorAllocatedIn.push_back(true);
        (static_cast<EqualityConstraint*>(ec))->createEqualityConstraint(nsdsxml->getEqualityConstraintXML(nbECtab[i]));
      }
      else if ((nsdsxml->getEqualityConstraintXML(nbECtab[i]))->getType() == LINEAR_TIME_INVARIANT_EC_TAG)
      {
        ec = new LinearTIEC();
        ecVector.push_back(ec);
        isEcVectorAllocatedIn.push_back(true);
        (static_cast<LinearTIEC*>(ec))->createEqualityConstraint(nsdsxml->getEqualityConstraintXML(nbECtab[i]));
      }
      else if ((nsdsxml->getEqualityConstraintXML(nbECtab[i]))->getType() == LAGRANGIAN_EC_TAG)
      {
        ec = new LagrangianEC();
        ecVector.push_back(ec);
        isEcVectorAllocatedIn.push_back(true);
        (static_cast<LagrangianEC*>(ec))->createEqualityConstraint(nsdsxml->getEqualityConstraintXML(nbECtab[i]));
      }
      else if ((nsdsxml->getEqualityConstraintXML(nbECtab[i]))->getType() == LAGRANGIAN_LINEAR_EC_TAG)
      {
        ec = new LagrangianLinearEC();
        ecVector.push_back(ec);
        isEcVectorAllocatedIn.push_back(true);
        (static_cast<LagrangianLinearEC*>(ec))->createEqualityConstraint(nsdsxml->getEqualityConstraintXML(nbECtab[i]));
      }
      else RuntimeException::selfThrow("NonSmoothDynamicalSystem:: xml constructor, wrong kind of algebraic constraints");
    }
  }
  else RuntimeException::selfThrow("NonSmoothDynamicalSystem:: xml constructor, xml file=NULL");
}

// constructor from data
// TODO

// Other constructors
NonSmoothDynamicalSystem::NonSmoothDynamicalSystem(const bool& bvp):
  BVP(bvp), nsdsxml(NULL)
{}

NonSmoothDynamicalSystem::NonSmoothDynamicalSystem(const string& type):
  BVP(false), nsdsxml(NULL)
{
  if (type == "BVP")
    BVP = true; // else default value = IVP
}

NonSmoothDynamicalSystem::~NonSmoothDynamicalSystem()
{
  IN("NonSmoothDynamicalSystem::~NonSmoothDynamicalSystem\n");
  int i;

  if (DSVector.size() > 0)
  {
    for (i = 0; i < DSVector.size(); i++)
    {
      if (isDSVectorAllocatedIn[i])
      {
        delete DSVector[i];
        DSVector[i] = NULL;
      }
    }
    DSVector.clear();
  }

  if (interactionVector.size() > 0)
  {
    for (i = 0; i < interactionVector.size(); i++)
    {
      if (isInteractionVectorAllocatedIn[i])
      {
        delete interactionVector[i];
        interactionVector[i] = NULL;
      }
    }
    interactionVector.clear();
  }
  OUT("NonSmoothDynamicalSystem::~NonSmoothDynamicalSystem\n");
}

// --- GETTERS/SETTERS ---

DynamicalSystem* NonSmoothDynamicalSystem::getDynamicalSystemPtr(const int& nb) const
{
  IN("DynamicalSystem* NonSmoothDynamicalSystem::getDynamicalSystem(int nb)\n");
  if (nb < DSVector.size())
  {
    OUT("DynamicalSystem* NonSmoothDynamicalSystem::getDynamicalSystem(int nb)\n");
    return DSVector[nb];
  }
  RuntimeException::selfThrow("NonSmoothDynamicalSystem - getDynamicalSystem : \'nb\' is out of range");
}

DynamicalSystem* NonSmoothDynamicalSystem::getDynamicalSystemPtrNumber(const int& nb) const
{
  for (unsigned int i = 0; i < DSVector.size(); i++)
  {
    if (DSVector[i]->getNumber() == nb)
    {
      return DSVector[i];
    }
  }
  RuntimeException::selfThrow("NonSmoothDynamicalSystem::getDynamicalSystemOnNumber : DSVector[i] == NULL");
}

void NonSmoothDynamicalSystem::setDynamicalSystems(const std::vector<DynamicalSystem*>& newVect)
{
  for (int i = 0; i < DSVector.size(); i++)
  {
    if (isDSVectorAllocatedIn[i]) delete DSVector[i];
  }
  DSVector = newVect;
  isDSVectorAllocatedIn.clear();
  isDSVectorAllocatedIn.resize(newVect.size(), false);
}

Interaction* NonSmoothDynamicalSystem::getInteractionPtr(const int& nb) const
{
  if (nb < interactionVector.size())
  {
    return interactionVector[nb];
  }
  RuntimeException::selfThrow("NonSmoothDynamicalSystem - getInteraction : \'nb\' is out of range");
}

Interaction* NonSmoothDynamicalSystem::getInteractionPtrNumber(const int& nb) const
{
  for (unsigned int i = 0; i < interactionVector.size(); i++)
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

void NonSmoothDynamicalSystem::setInteractions(const std::vector<Interaction*>& newVect)
{
  for (int i = 0; i < interactionVector.size(); i++)
  {
    if (isInteractionVectorAllocatedIn[i]) delete interactionVector[i];
  }
  interactionVector = newVect;
  isInteractionVectorAllocatedIn.clear();
  isInteractionVectorAllocatedIn.resize(newVect.size(), false);
}

EqualityConstraint* NonSmoothDynamicalSystem::getEqualityConstraintPtr(const int& i) const
{
  if (i < ecVector.size())
  {
    return ecVector[i];
  }
  RuntimeException::selfThrow("NonSmoothDynamicalSystem - getEqualityConstraint : \'i\' is out of range");
}

void NonSmoothDynamicalSystem::setEqualityConstraints(const std::vector<EqualityConstraint*>& newVect)
{
  for (int i = 0; i < ecVector.size(); i++)
  {
    if (isEcVectorAllocatedIn[i]) delete ecVector[i];
  }
  ecVector = newVect;
  isEcVectorAllocatedIn.clear();
  isEcVectorAllocatedIn.resize(newVect.size(), false);
}

// --- OTHER FUNCTIONS ---

bool NonSmoothDynamicalSystem::hasDynamicalSystemNumber(const int& nb) const
{
  for (unsigned int i = 0; i < DSVector.size(); i++)
  {
    if (DSVector[i]->getNumber() == nb)
    {
      return true;
    }
  }
  return false;
}
bool NonSmoothDynamicalSystem::hasInteractionNumber(const int& nb) const
{
  for (unsigned int i = 0; i < interactionVector.size(); i++)
  {
    if (interactionVector[i] != NULL)
      if (interactionVector[i]->getNumber() == nb)
      {
        return true;
      }
  }
  return false;
}
// Xml management functions

void NonSmoothDynamicalSystem::saveNSDSToXML()
{
  IN("NonSmoothDynamicalSystem::saveNSDSToXML\n");

  if (nsdsxml != NULL)
  {
    int size, i;

    nsdsxml->setBVP(BVP);// no need to change the value of BVP, it mustn't change anyway

    size = DSVector.size();
    for (i = 0; i < size; i++)
    {
      if (DSVector[i]->getType() == LNLDS)
        (static_cast<LagrangianDS*>(DSVector[i]))->saveDSToXML();
      else if (DSVector[i]->getType() == LTIDS)
        (static_cast<LagrangianLinearTIDS*>(DSVector[i]))->saveDSToXML();
      else if (DSVector[i]->getType() == LDS)
        (static_cast<LinearSystemDS*>(DSVector[i]))->saveDSToXML();
      else if (DSVector[i]->getType() == NLDS)
        DSVector[i]->saveDSToXML();
      else RuntimeException::selfThrow("NonSmoothDynamicalSystem::saveToXML - bad kind of DynamicalSystem");
    }

    size = ecVector.size();
    for (i = 0; i < size; i++)
    {
      if (ecVector[i]->getType() == LINEAREC)
        (static_cast<LinearEC*>(ecVector[i]))->saveEqualityConstraintToXML();
      else if (ecVector[i]->getType() == LINEARTIEC)
        (static_cast<LinearTIEC*>(ecVector[i]))->saveEqualityConstraintToXML();
      else if (ecVector[i]->getType() == NLINEAREC)
        ecVector[i]->saveEqualityConstraintToXML();
      else if (ecVector[i]->getType() == LAGRANGIANEC)
        (static_cast<LagrangianEC*>(ecVector[i]))->saveEqualityConstraintToXML();
      else if (ecVector[i]->getType() == LAGRANGIANLINEAREC)
        (static_cast<LagrangianLinearEC*>(ecVector[i]))->saveEqualityConstraintToXML();
      else RuntimeException::selfThrow("NonSmoothDynamicalSystem::saveToXML - bad kind of EqualityConstraint");
    }

    size = interactionVector.size();
    for (i = 0; i < size; i++)
      interactionVector[i]->saveInteractionToXML();
  }
  else RuntimeException::selfThrow("NonSmoothDynamicalSystem::saveNSDSToXML - The NSDSXML object doesn't exists");

  OUT("NonSmoothDynamicalSystem::saveNSDSToXML\n");
}

void NonSmoothDynamicalSystem::display() const
{
  IN("NonSmoothDynamicalSystem::display\n");

  cout << "| this = " << this << endl;
  cout << "| BVP = " << BVP << endl;
  cout << "| &nsdsxml = " << nsdsxml << endl;
  cout << "|===========================" << endl;

  OUT("NonSmoothDynamicalSystem::display\n");
}

void NonSmoothDynamicalSystem::addDynamicalSystem(DynamicalSystem *ds)//, BoundaryCondition* bc)
{
  IN("NonSmoothDynamicalSystem::addDynamicalSystem\n");
  ds->setNSDSPtr(this);
  DSVector.push_back(ds);
  isDSVectorAllocatedIn.push_back(true);
  OUT("NonSmoothDynamicalSystem::addDynamicalSystem\n");
}


void NonSmoothDynamicalSystem::addInteraction(Interaction *inter)
{
  Interaction* interTmp;
  interTmp = new Interaction(*inter);
  interactionVector.push_back(interTmp);
  isInteractionVectorAllocatedIn.push_back(true);
}

Interaction* NonSmoothDynamicalSystem::addInteraction(const int& number, const int& nInter, vector<int>* status, vector<DynamicalSystem*>* dsConcerned)
{
  if (!hasInteractionNumber(number))
  {
    Interaction *interTmp = new Interaction("noId", number, nInter, status, dsConcerned);
    interactionVector.push_back(interTmp);
    isInteractionVectorAllocatedIn.push_back(true);
    return interTmp;
  }
  else
    RuntimeException::selfThrow("NonSmoothDynamicalSystem::addInteraction : Interaction already declared; number: " + number);
}

void NonSmoothDynamicalSystem::addEqualityConstraint(EqualityConstraint* ec)
{
  ecVector.push_back(ec);
  isEcVectorAllocatedIn.push_back(true);
}


double NonSmoothDynamicalSystem::nsdsConvergenceIndicator()
{
  // calculate the max value of all DS convergence indicators
  double convergenceIndicator = (* DSVector.begin()) -> dsConvergenceIndicator();
  vector <DynamicalSystem*>::iterator iter;
  for (iter = DSVector.begin(); iter != DSVector.end(); ++iter)
  {
    double dsIndic = (*iter)->dsConvergenceIndicator();
    if (convergenceIndicator > dsIndic) convergenceIndicator == dsIndic;
  }
  return(convergenceIndicator);
}
