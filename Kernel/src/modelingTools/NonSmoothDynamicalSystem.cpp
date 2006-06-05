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

// xml constuctor
NonSmoothDynamicalSystem::NonSmoothDynamicalSystem(NonSmoothDynamicalSystemXML* newNsdsxml):
  BVP(false), topology(NULL), nsdsxml(newNsdsxml), isTopologyAllocatedIn(false)
{
  if (nsdsxml == NULL)
    RuntimeException::selfThrow("NonSmoothDynamicalSystem:: xml constructor, xml file=NULL");

  // === DS Vector fill-in: we sweep the list of DSXML and for each of them, add a new DynamicalSystem in the set allDS. ===
  SetOfDSXML dsList = nsdsxml->getDynamicalSystemsXML();
  SetOfDSXMLIt it;
  CheckInsertDS checkDS;
  string type;
  for (it = dsList.begin(); it != dsList.end(); ++it)
  {
    type = (*it)->getType();
    if (type  == LAGRANGIAN_NON_LINEARDS_TAG)  // LagrangianDS
      checkDS = allDS.insert(new LagrangianDS(*it , this));
    else if (type == LAGRANGIAN_TIDS_TAG)  // Lagrangian Linear Time Invariant
      checkDS = allDS.insert(new LagrangianLinearTIDS(*it, this));
    else if (type == LINEAR_DS_TAG)  // Linear DS
      checkDS = allDS.insert(new LinearDS(*it, this));
    else if (type == LINEAR_TIDS_TAG)  // Linear DS
      checkDS = allDS.insert(new LinearTIDS(*it, this));
    else if (type == NON_LINEAR_DS_TAG)  // Non linear DS
      checkDS = allDS.insert(new DynamicalSystem(*it, this));
    else RuntimeException::selfThrow("NonSmoothDynamicalSystem::xml constructor, wrong Dynamical System type" + type);
    // checkDS.first is an iterator that points to the DS inserted into the set.
    isDSAllocatedIn[*(checkDS.first)] = true ;
  }

  // ===  The same process is applied for Interactions ===
  SetOfInteractionsXML  interactionsList = nsdsxml->getInteractionsXML();
  SetOfInteractionsXMLIt it2;
  CheckInsertInteraction checkInter;
  for (it2 = interactionsList.begin(); it2 != interactionsList.end(); ++it2)
  {
    checkInter = allInteractions.insert(new Interaction(*it2, this));
    // checkInter.first is an iterator that points to the Interaction inserted into the set.
    isInteractionAllocatedIn[*(checkInter.first)] = true;
  }

  if (allInteractions.size() == 0) cout << " /!\\ Warning: you do not defined any interaction in the NonSmoothDymamicalSystem (Constructor: xml). /!\\ " << endl;

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

// copy constructor
NonSmoothDynamicalSystem::NonSmoothDynamicalSystem(const NonSmoothDynamicalSystem& nsds):
  BVP(false), topology(NULL), nsdsxml(NULL), isTopologyAllocatedIn(false)
{
  BVP = nsds.isBVP();

  // === copy of DynamicalSystems ===
  allDS = nsds.getDynamicalSystems(); // Warning!! This is a false copy, since pointers links remains between DS of each set
  DSIterator it;
  for (it = allDS.begin(); it != allDS.end(); ++it)
    isDSAllocatedIn[*it] = false;

  // === copy of Interactions ===
  allInteractions = nsds.getInteractions(); // Warning!! This is a false copy, since pointers links remains between Interactions of each set
  InteractionsIterator it2;
  for (it2 = allInteractions.begin(); it2 != allInteractions.end(); ++it2)
    isInteractionAllocatedIn[*it2] = false;

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
// \todo add a  constructor from data ?

NonSmoothDynamicalSystem::~NonSmoothDynamicalSystem()
{
  // == delete DS ==
  DSIterator it;
  for (it = allDS.begin(); it != allDS.end(); ++it)
  {
    if (isDSAllocatedIn[*it]) delete *it;
  }

  allDS.clear();
  isDSAllocatedIn.clear();

  // == delete Interactions ==
  InteractionsIterator it2;
  for (it2 = allInteractions.begin(); it2 != allInteractions.end(); ++it2)
  {
    if (isInteractionAllocatedIn[*it2]) delete *it2;
  }

  allInteractions.clear();
  isInteractionAllocatedIn.clear();

  if (ecVector.size() > 0)
  {
    for (unsigned int i = 0; i < ecVector.size(); i++)
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

// === DynamicalSystems management ===

DynamicalSystem* NonSmoothDynamicalSystem::getDynamicalSystemPtr(const int& nb) const
{
  // Mind that DS are sorted in a growing order according to their id number.
  DSIterator it = allDS.begin();
  for (int i = 0; i < nb; ++i)
    it++;

  if (it == allDS.end())
    RuntimeException::selfThrow("NonSmoothDynamicalSystem - getDynamicalSystemPtr(nb) : nb is out of range");

  return *it;
}

DynamicalSystem* NonSmoothDynamicalSystem::getDynamicalSystemPtrNumber(const int& nb) const
{
  if (! allDS.isDSIn(nb)) // if ds number nb is not in the set ...
    RuntimeException::selfThrow("NonSmoothDynamicalSystem::getDynamicalSystemOnNumber(nb), DS number nb is not in the set.");

  return allDS.getDynamicalSystem(nb);
}

void NonSmoothDynamicalSystem::setDynamicalSystems(const DSSet& newVect)
{
  // clear old set
  DSIterator it;
  for (it = allDS.begin(); it != allDS.end(); ++it)
  {
    if (isDSAllocatedIn[*it]) delete *it;
  }

  allDS.clear();
  isDSAllocatedIn.clear();

  // copy the new one
  allDS = newVect;
  for (it = allDS.begin(); it != allDS.end(); ++it)
    isDSAllocatedIn[*it] = false;

  topology->setUpToDate(false);
}

const bool NonSmoothDynamicalSystem::hasDynamicalSystemNumber(const int& nb) const
{
  return allDS.isDSIn(nb);
}

const bool NonSmoothDynamicalSystem::hasDynamicalSystem(DynamicalSystem* ds) const
{
  return allDS.isDSIn(ds);
}

// === Interactions management ===

Interaction* NonSmoothDynamicalSystem::getInteractionPtr(const int& nb) const
{
  // Mind that Interactions are sorted in a growing order according to their id number.
  InteractionsIterator it = allInteractions.begin();
  for (int i = 0; i < nb; ++i)
    it++;

  if (it == allInteractions.end())
    RuntimeException::selfThrow("NonSmoothDynamicalSystem - getInteractionPtr(nb) : nb is out of range");

  return *it;
}

Interaction* NonSmoothDynamicalSystem::getInteractionPtrNumber(const int& nb) const
{
  if (! allInteractions.isInteractionIn(nb)) // if Interaction number nb is not in the set ...
    RuntimeException::selfThrow("NonSmoothDynamicalSystem::getInteractionOnNumber(nb), Interaction number nb is not in the set.");

  return allInteractions.getInteraction(nb);
}

void NonSmoothDynamicalSystem::setInteractions(const InteractionsSet& newVect)
{
  // clear old set
  InteractionsIterator it;
  for (it = allInteractions.begin(); it != allInteractions.end(); ++it)
  {
    if (isInteractionAllocatedIn[*it]) delete *it;
  }

  allInteractions.clear();
  isInteractionAllocatedIn.clear();

  // copy the new one
  allInteractions = newVect;
  for (it = allInteractions.begin(); it != allInteractions.end(); ++it)
    isInteractionAllocatedIn[*it] = false;

  topology->setUpToDate(false);
}

const bool NonSmoothDynamicalSystem::hasInteractionNumber(const int& nb) const
{
  return allInteractions.isInteractionIn(nb);
}

const bool NonSmoothDynamicalSystem::hasInteraction(Interaction* inter) const
{
  return allInteractions.isInteractionIn(inter);
}

// === EqualityConstraints management ===

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

// === Xml management functions ===

void NonSmoothDynamicalSystem::saveNSDSToXML()
{
  if (nsdsxml != NULL)
  {
    int size, i;

    nsdsxml->setBVP(BVP);// no need to change the value of BVP, it mustn't change anyway

    DSIterator it;
    for (it = allDS.begin(); it != allDS.end(); ++it)
    {
      if ((*it)->getType() == LNLDS)
        (static_cast<LagrangianDS*>((*it)))->saveDSToXML();
      else if ((*it)->getType() == LTIDS)
        (static_cast<LagrangianLinearTIDS*>((*it)))->saveDSToXML();
      else if ((*it)->getType() == LDS)
        (static_cast<LinearDS*>((*it)))->saveDSToXML();
      else if ((*it)->getType() == LITIDS)
        (static_cast<LinearDS*>((*it)))->saveDSToXML();
      else if ((*it)->getType() == NLDS)
        (*it)->saveDSToXML();
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

    InteractionsIterator it2;
    for (it2 = allInteractions.begin(); it2 != allInteractions.end(); ++it2)
      (*it2)->saveInteractionToXML();
  }
  else RuntimeException::selfThrow("NonSmoothDynamicalSystem::saveNSDSToXML - The NonSmoothDynamicalSystemXML object doesn't exists");
}

void NonSmoothDynamicalSystem::display() const
{
  cout << " ===== Non Smooth Dynamical System display ===== " << endl;
  cout << "---> isBVP = " << BVP << endl;
  allDS.display();
  allInteractions.display();
  cout << "===================================================" << endl;
}

void NonSmoothDynamicalSystem::addDynamicalSystemPtr(DynamicalSystem *ds)
{
  ds->setNonSmoothDynamicalSystemPtr(this);
  allDS.insert(ds);
  isDSAllocatedIn[ds] = false;
  topology->setUpToDate(false);
}

void NonSmoothDynamicalSystem::addInteractionPtr(Interaction *inter)
{
  inter->setNonSmoothDynamicalSystemPtr(this);
  allInteractions.insert(inter);
  isInteractionAllocatedIn[inter] = false;
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
  double convergenceIndicator = (* allDS.begin()) -> dsConvergenceIndicator();
  DSIterator iter;
  for (iter = allDS.begin(); iter != allDS.end(); ++iter)
  {
    double dsIndic = (*iter)->dsConvergenceIndicator();
    if (convergenceIndicator > dsIndic) convergenceIndicator = dsIndic;
  }
  return(convergenceIndicator);
}
