/* Siconos-Kernel version 1.1.4, Copyright INRIA 2005-2006.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
#include "NonSmoothDynamicalSystem.h"

// includes to be deleted thanks to factories
#include "LagrangianDS.h"
#include "LagrangianLinearTIDS.h"
#include "LinearDS.h"
#include "LinearTIDS.h"
#include "LinearEC.h"
#include "LinearTIEC.h"
#include "LagrangianEC.h"
#include "LagrangianLinearEC.h"

using namespace std;

// --- CONSTRUCTORS/DESTRUCTOR ---

// Default constructor (isBvp is optional, default = false)
NonSmoothDynamicalSystem::NonSmoothDynamicalSystem(const bool& isBvp):
  BVP(isBvp), topology(NULL), nsdsxml(NULL), isTopologyAllocatedIn(false)
{
  topology = new Topology(this); // \todo use a copy constructor for topology?
  isTopologyAllocatedIn = true;
}

// copy constructor
NonSmoothDynamicalSystem::NonSmoothDynamicalSystem(const NonSmoothDynamicalSystem& nsds):
  BVP(false), topology(NULL), nsdsxml(NULL), isTopologyAllocatedIn(false)
{
  BVP = nsds.isBVP();
  unsigned int i;
  DSVector.resize(nsds.getDynamicalSystems().size(), NULL);

  isDSVectorAllocatedIn.resize(DSVector.size(), false);
  DynamicalSystem * dsTmp;
  for (i = 0; i < DSVector.size(); i++)
  {
    dsTmp = nsds.getDynamicalSystemPtr(i);
    if (dsTmp != NULL)
    {
      // \todo use factories to improve this copy.
      string type = dsTmp ->getType();
      if (type ==  NLDS)
        DSVector[i] = new DynamicalSystem(*dsTmp);
      else if (type ==  LDS)
        DSVector[i] = new LinearDS(*dsTmp);
      else if (type ==  LITIDS)
        DSVector[i] = new LinearTIDS(*dsTmp);
      else if (type ==  LNLDS)
        DSVector[i] = new LagrangianDS(*dsTmp);
      else if (type ==  LTIDS)
        DSVector[i] = new LagrangianLinearTIDS(*dsTmp);
      else
        RuntimeException::selfThrow("NonSmoothDynamicalSystem::copy constructor, unknown Dynamical system type:" + type);
      isDSVectorAllocatedIn[i] = true;
    }
  }
  interactionVector.resize(nsds.getInteractions().size(), NULL);
  isInteractionVectorAllocatedIn.resize(interactionVector.size(), false);
  Interaction * interTmp;
  for (i = 0; i < interactionVector.size(); i++)
  {
    interTmp = nsds.getInteractionPtr(i);
    if (interTmp != NULL)
    {
      interactionVector[i] = new Interaction(*interTmp);
      isInteractionVectorAllocatedIn[i] = true;
    }
  }

  topology = new Topology(this); // \todo use a copy constructor for topology?
  isTopologyAllocatedIn = true;

  /*  \todo: add this part when EC will be well-implemented
  ecVector.resize(nsds.getEqualityConstraints.size(), NULL);
  isEcVectorAllocatedIn.resize(ecVector.size(), false);
  equalityConstraint * ecTmp;
  for(i=0; i< ecVector.size();i++)
    {
      ecTmp = nsds.getEqualityConstraintPtr(i);
      if(ecTmp != NULL)
  {
    string typeEc = ecTmp->getType() ;
    if( typeEc == LINEAREC )
      ecVector[i] = new LinearEC( *ecTmp) ;
    else if( typeEc == LINEARTIEC )
      ecVector[i] = new LinearTIEC( *ecTmp) ;
    else if( typeEc == NLINEAREC )
      ecVector[i] = new NonLinearEC( *ecTmp) ;
    else if( typeEc == LAGRANGIANEC )
      ecVector[i] = new LagrangianEC( *ecTmp) ;
    else if( typeEc == LAGRANGIANLINEAREC )
      ecVector[i] = new LagrangianLinearEC( *ecTmp) ;
    else
      RuntimeException::selfThrow("NonSmoothDynamicalSystem::copy constructor, unknown Equality constraint type:"+typeEc);
    isEcVectorAllocatedIn[i] = true;
  }
    }
  */

  // Warning: xml link is not copied.
}

// xml constuctor
NonSmoothDynamicalSystem::NonSmoothDynamicalSystem(NonSmoothDynamicalSystemXML* newNsdsxml):
  BVP(false), nsdsxml(newNsdsxml), isTopologyAllocatedIn(false)
{
  if (nsdsxml != NULL)
  {
    unsigned int i = 0;
    // DS Vector fill-in
    vector<int> nbDStab = nsdsxml->getDSNumbers();
    unsigned int size = nbDStab.size();
    DSVector.resize(size, NULL);
    isDSVectorAllocatedIn.resize(size, false);

    for (i = 0; i < size; i++)
    {
      string type = (nsdsxml->getDynamicalSystemXML(nbDStab[i]))->getType();
      if (type  == LAGRANGIAN_NON_LINEARDS_TAG)  // LagrangianDS
      {
        DSVector[i] = new LagrangianDS(nsdsxml->getDynamicalSystemXML(nbDStab[i]), this);
        isDSVectorAllocatedIn[i] = true;
      }
      else if (type == LAGRANGIAN_TIDS_TAG)  // Lagrangian Linear Time Invariant
      {
        DSVector[i] = new LagrangianLinearTIDS(nsdsxml->getDynamicalSystemXML(nbDStab[i]), this);
        isDSVectorAllocatedIn[i] = true;
      }
      else if (type == LINEAR_DS_TAG)  // Linear DS
      {
        DSVector[i] = new LinearDS(nsdsxml->getDynamicalSystemXML(nbDStab[i]), this);
        isDSVectorAllocatedIn[i] = true;
      }
      else if (type == LINEAR_TIDS_TAG)  // Linear DS
      {
        DSVector[i] = new LinearTIDS(nsdsxml->getDynamicalSystemXML(nbDStab[i]), this);
        isDSVectorAllocatedIn[i] = true;
      }
      else if (type == NON_LINEAR_DS_TAG)  // Non linear DS
      {
        DSVector[i] = new DynamicalSystem(nsdsxml->getDynamicalSystemXML(nbDStab[i]), this);
        isDSVectorAllocatedIn[i] = true;
      }
      else RuntimeException::selfThrow("NonSmoothDynamicalSystem::xml constructor, wrong Dynamical System type" + type);
    }

    // Interaction vector fill-in
    vector<int> nbInteractionTab = nsdsxml->getInteractionNumbers();
    //int ds1, ds2;
    unsigned int sizeInter = nbInteractionTab.size();
    interactionVector.resize(sizeInter, NULL);
    isInteractionVectorAllocatedIn.resize(sizeInter, false);
    for (i = 0; i < sizeInter; i++)
    {
      // Creation of the Interaction
      interactionVector[i] = new Interaction(nsdsxml->getInteractionXML(nbInteractionTab[i]), this);
      isInteractionVectorAllocatedIn[i] = true;
    }
    if (nbInteractionTab.size() == 0) cout << " /!\\ Warning: you do not defined any interaction in the NonSmoothDymamicalSystem (Constructor: xml). /!\\ " << endl;

    // built topology:
    topology = new Topology(this);
    isTopologyAllocatedIn = true;

    // Algebraic constraints fill-in
    // get all the EqualityConstraintXML objects then create the EqualityConstraint for this EqualityConstraintXML
    // \todo: uncomment this part when EC will be well-implemented
    /*      vector<int> nbECtab = nsdsxml->getEqualityConstraintNumbers();
    for( i=0; i<nbECtab.size(); i++ )
    {
    EqualityConstraint *ec;
    if((nsdsxml->getEqualityConstraintXML( nbECtab[i]) )->getType() == LINEAR_EC_TAG )
    {
      ec = new LinearEC();
      ecVector.push_back( ec );
      isEcVectorAllocatedIn.push_back(true);
      (static_cast<LinearEC*>(ec))->createEqualityConstraint( nsdsxml->getEqualityConstraintXML( nbECtab[i]) );
    }
    else if((nsdsxml->getEqualityConstraintXML( nbECtab[i]) )->getType() == NON_LINEAR_EC_TAG )
    {
      ec = new EqualityConstraint();
      ecVector.push_back( ec );
      isEcVectorAllocatedIn.push_back(true);
      (static_cast<EqualityConstraint*>(ec))->createEqualityConstraint( nsdsxml->getEqualityConstraintXML( nbECtab[i]) );
    }
    else if((nsdsxml->getEqualityConstraintXML( nbECtab[i]) )->getType() == LINEAR_TIME_INVARIANT_EC_TAG )
    {
      ec = new LinearTIEC();
      ecVector.push_back( ec );
      isEcVectorAllocatedIn.push_back(true);
      (static_cast<LinearTIEC*>(ec))->createEqualityConstraint( nsdsxml->getEqualityConstraintXML( nbECtab[i]) );
    }
    else if((nsdsxml->getEqualityConstraintXML( nbECtab[i]) )->getType() == LAGRANGIAN_EC_TAG )
    {
      ec = new LagrangianEC();
      ecVector.push_back( ec );
      isEcVectorAllocatedIn.push_back(true);
      (static_cast<LagrangianEC*>(ec))->createEqualityConstraint( nsdsxml->getEqualityConstraintXML( nbECtab[i]) );
    }
    else if((nsdsxml->getEqualityConstraintXML( nbECtab[i]) )->getType() == LAGRANGIAN_LINEAR_EC_TAG )
    {
      ec = new LagrangianLinearEC();
      ecVector.push_back( ec );
      isEcVectorAllocatedIn.push_back(true);
      (static_cast<LagrangianLinearEC*>(ec))->createEqualityConstraint( nsdsxml->getEqualityConstraintXML( nbECtab[i]) );
    }
    else RuntimeException::selfThrow("NonSmoothDynamicalSystem:: xml constructor, wrong kind of algebraic constraints");
    }
    */
  }
  else RuntimeException::selfThrow("NonSmoothDynamicalSystem:: xml constructor, xml file=NULL");
}

// \todo add a  constructor from data ?

NonSmoothDynamicalSystem::~NonSmoothDynamicalSystem()
{
  unsigned int i;

  if (DSVector.size() > 0)
  {
    for (i = 0; i < DSVector.size(); i++)
    {
      if (isDSVectorAllocatedIn[i]) delete DSVector[i];
      DSVector[i] = NULL;
    }
  }
  DSVector.clear();

  for (i = 0; i < interactionVector.size(); i++)
  {
    if (isInteractionVectorAllocatedIn[i])
      delete interactionVector[i];
    interactionVector[i] = NULL;
  }
  interactionVector.clear();

  if (ecVector.size() > 0)
  {
    for (i = 0; i < ecVector.size(); i++)
    {
      if (isEcVectorAllocatedIn[i]) delete ecVector[i];
      ecVector[i] = NULL;
    }
  }
  ecVector.clear();

  if (isTopologyAllocatedIn)
  {
    delete topology;
    topology = NULL;
  }
}

// --- GETTERS/SETTERS ---

DynamicalSystem* NonSmoothDynamicalSystem::getDynamicalSystemPtr(const int& nb) const
{
  if ((unsigned int)nb >= DSVector.size())
    RuntimeException::selfThrow("NonSmoothDynamicalSystem - getDynamicalSystem : \'nb\' is out of range");
  return DSVector[nb];
}

DynamicalSystem* NonSmoothDynamicalSystem::getDynamicalSystemPtrNumber(const int& nb) const
{
  unsigned int i = 0;
  while (DSVector[i]->getNumber() != nb && i < DSVector.size())
  {
    i++ ;
  }
  if (i > DSVector.size())
    RuntimeException::selfThrow("NonSmoothDynamicalSystem::getDynamicalSystemOnNumber, out of range or chosen DSVector == NULL");
  return DSVector[i];
}

void NonSmoothDynamicalSystem::setDynamicalSystems(const std::vector<DynamicalSystem*>& newVect)
{
  for (unsigned int i = 0; i < DSVector.size(); i++)
  {
    if (isDSVectorAllocatedIn[i]) delete DSVector[i];
  }
  DSVector = newVect;
  isDSVectorAllocatedIn.clear();
  isDSVectorAllocatedIn.resize(newVect.size(), false);
}

Interaction* NonSmoothDynamicalSystem::getInteractionPtr(const int& nb) const
{
  if ((unsigned int)nb >= interactionVector.size())
    RuntimeException::selfThrow("NonSmoothDynamicalSystem - getInteraction : \'nb\' is out of range");
  return interactionVector[nb];
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
  for (unsigned int i = 0; i < interactionVector.size(); i++)
  {
    if (isInteractionVectorAllocatedIn[i]) delete interactionVector[i];
  }
  interactionVector = newVect;
  isInteractionVectorAllocatedIn.clear();
  isInteractionVectorAllocatedIn.resize(newVect.size(), false);

  // the topology should be updated
  topology->setUpToDate(false);

}

EqualityConstraint* NonSmoothDynamicalSystem::getEqualityConstraintPtr(const int& i) const
{
  if ((unsigned int)i >= ecVector.size())
    RuntimeException::selfThrow("NonSmoothDynamicalSystem - getEqualityConstraint : \'i\' is out of range");
  return ecVector[i];
}

void NonSmoothDynamicalSystem::setEqualityConstraints(const std::vector<EqualityConstraint*>& newVect)
{
  for (unsigned int i = 0; i < ecVector.size(); i++)
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
  bool hasDS = false;
  for (unsigned int i = 0; i < DSVector.size(); i++)
  {
    if (DSVector[i]->getNumber() == nb)
      hasDS = true;
  }
  return hasDS;
}

bool NonSmoothDynamicalSystem::hasInteractionNumber(const int& nb) const
{
  bool hasInter = false;
  for (unsigned int i = 0; i < interactionVector.size(); i++)
  {
    if (interactionVector[i] != NULL)
      if (interactionVector[i]->getNumber() == nb)
        hasInter = true;
  }
  return hasInter;
}
// Xml management functions

void NonSmoothDynamicalSystem::saveNSDSToXML()
{
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
        (static_cast<LinearDS*>(DSVector[i]))->saveDSToXML();
      else if (DSVector[i]->getType() == LITIDS)
        (static_cast<LinearDS*>(DSVector[i]))->saveDSToXML();
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
  else RuntimeException::selfThrow("NonSmoothDynamicalSystem::saveNSDSToXML - The NonSmoothDynamicalSystemXML object doesn't exists");
}

void NonSmoothDynamicalSystem::display() const
{
  cout << " ===== Non Smooth Dynamical System display ===== " << endl;
  cout << "| Adress of this object = " << this << endl;
  cout << "| isBVP = " << BVP << endl;
  if (nsdsxml != NULL)
    cout << "| &nsdsxml (xml link) = " << nsdsxml << endl;
  else
    cout << "| &nsdsxml (xml link) -> NULL " << endl;
  cout << "Dynamical systems list:" << endl;
  vector<DynamicalSystem*>::const_iterator itDS;
  for (itDS = DSVector.begin(); itDS != DSVector.end(); itDS++)
  {
    if ((*itDS) != NULL)
      (*itDS)->display();
    else
      cout << " DS-> NULL " << endl;
  }
  cout << "Interactions of this system:" << endl;
  vector<Interaction*>::const_iterator it;
  for (it = interactionVector.begin(); it != interactionVector.end(); it++)
  {
    if ((*it) != NULL)
      (*it)->display();
    else
      cout << " interaction -> NULL " << endl;
  }
  cout << "===================================================" << endl;
}

void NonSmoothDynamicalSystem::addDynamicalSystem(DynamicalSystem *ds)
{
  ds->setNonSmoothDynamicalSystemPtr(this);
  DSVector.push_back(ds);
  isDSVectorAllocatedIn.push_back(false);
}

void NonSmoothDynamicalSystem::addInteraction(Interaction *inter)
{
  interactionVector.push_back(inter);
  isInteractionVectorAllocatedIn.push_back(false);
  // the topology should be updated
  topology->setUpToDate(false);
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
    if (convergenceIndicator > dsIndic) convergenceIndicator = dsIndic;
  }
  return(convergenceIndicator);
}
