/* Siconos-Kernel version 1.2.0, Copyright INRIA 2005-2006.
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
#include "Interaction.h"

// includes to be deleted thanks to factories
#include "LinearTIR.h"
#include "LagrangianLinearR.h"
#include "LagrangianR.h"
#include "ComplementarityConditionNSL.h"
#include "RelayNSL.h"
#include "NewtonImpactNSL.h"
#include "NewtonImpactFrictionNSL.h"

using namespace std;

// --- CONSTRUCTORS ---

// Copy constructor
Interaction::Interaction(const Interaction& newI):
  id("interactionCopy"), number(-1), interactionSize(newI.getInteractionSize()),
  numberOfRelations(newI.getNumberOfRelations()), sizeOfDS(newI.getSizeOfDS()), involvedDS(),
  nslaw(NULL), relation(NULL), NSDS(newI.getNonSmoothDynamicalSystemPtr()), interactionxml(NULL), isRelationAllocatedIn(false), isNsLawAllocatedIn(false)
{
  // Memory allocation and copy for y and lambda

  cout << "Warning, Interaction copy constructor: you need to set an id and a number to the new interaction (default={interactionCopy,-1})." << endl;

  unsigned int numberOfDerivatives = (newI.getY()).size();

  y.resize(numberOfDerivatives, NULL);
  yOld.resize(numberOfDerivatives, NULL);
  lambda.resize(numberOfDerivatives, NULL);
  lambdaOld.resize(numberOfDerivatives, NULL);

  for (unsigned int i = 0; i < numberOfDerivatives; i++)
  {
    y[i] = new SimpleVector(*(newI.getYPtr(i)));
    yOld[i] = new SimpleVector(*(newI.getYOldPtr(i)));
    lambda[i] = new SimpleVector(*(newI.getLambdaPtr(i)));
    lambdaOld[i] = new SimpleVector(*(newI.getLambdaOldPtr(i)));
  }

  isYAllocatedIn.resize(numberOfDerivatives, true);
  isYOldAllocatedIn.resize(numberOfDerivatives, true);
  isLambdaAllocatedIn.resize(numberOfDerivatives, true);
  isLambdaOldAllocatedIn.resize(numberOfDerivatives, true);

  involvedDS = newI.getDynamicalSystems(); // Warning: this keeps pointers links between Dynamical Systems of new I and this.

  // Nslaw (warning! nslaw is an abstract class)
  string NslawType = newI.getNonSmoothLawPtr()->getType();
  if (NslawType ==  COMPLEMENTARITYCONDITIONNSLAW)
    nslaw = new ComplementarityConditionNSL();
  else if (NslawType == NEWTONIMPACTNSLAW)
    nslaw = new NewtonImpactNSL();
  else if (NslawType == NEWTONIMPACTFRICTIONNSLAW)
    nslaw = new NewtonImpactFrictionNSL();
  else if (NslawType == RELAYNSLAW)
    nslaw = new  RelayNSL();
  else RuntimeException::selfThrow("Interaction::copy constructor, unknown NSLAW type :" + nslaw->getType());
  isNsLawAllocatedIn = true;
  *nslaw = *newI.getNonSmoothLawPtr(); // \todo: add a copy constructor in NonSmoothLaw class.

  // Relation
  string relationType = newI.getRelationPtr()->getType();
  // -> call copy constructor
  if (relationType == RELATION)
    relation = new Relation(*(newI.getRelationPtr()), this);

  if (relationType == LINEARTIRELATION)
    relation = new LinearTIR(*(newI.getRelationPtr()), this);

  else if (relationType == LAGRANGIANLINEARRELATION)
    relation = new LagrangianLinearR(*(newI.getRelationPtr()), this);

  else if (relationType == LAGRANGIANRELATION)
    relation = new LagrangianR(*(newI.getRelationPtr()), this);

  else RuntimeException::selfThrow("Interaction::copy constructor, unknown relation type " + relation->getType());
  isRelationAllocatedIn = true;

  // \remark FP: we do not link xml object in the copy
}

// --- XML constructor ---
Interaction::Interaction(InteractionXML* interxml, NonSmoothDynamicalSystem * nsds):
  id("undefined"), number(0), interactionSize(0), numberOfRelations(0), sizeOfDS(0),
  nslaw(NULL), relation(NULL), NSDS(nsds), interactionxml(interxml),
  isRelationAllocatedIn(false), isNsLawAllocatedIn(false)
{
  if (interactionxml == NULL)
    RuntimeException::selfThrow("Interaction::xml constructor, xmlfile = NULL");

  // id and number
  if (interactionxml->hasId()) id = interactionxml->getId();
  number = interactionxml->getNumber();

  // interaction size
  interactionSize = interactionxml->getSize();

  // Memory allocation for y and lambda
  initialize();

  // download xml values for y[0] and lambda[0] (optional)
  if (interactionxml->hasY()) *(y[0]) = interactionxml->getY();
  if (interactionxml->hasLambda()) *(lambda[0]) = interactionxml->getLambda();

  // Old values are initialized with current values (for y and lambda)
  swapInMemory();

  // --- Dynamical Systems ---
  unsigned int sizeDS ;
  if (nsds != NULL)
  {
    // Get a list of DS concerned from xml

    if (interactionxml->hasAllDS())
      involvedDS = nsds->getDynamicalSystems(); // this results in a call to the copy constructor. Do not understand why ??

    else
    {
      // get numbers of DS involved in the interaction from the xml input file.
      vector<int> listDS;
      interactionxml->getDSNumbers(listDS);

      // get corresponding DS and insert them into the involvedDS set.
      sizeDS = listDS.size();
      for (unsigned int i = 0; i < sizeDS; i++)
        involvedDS.insert(nsds->getDynamicalSystemPtrNumber(listDS[i]));
    }

    computeSizeOfDS();
  }
  else cout << "Interaction constructor, warning: no dynamical systems linked to the interaction!" << endl;

  // --- Non smooth law ---
  string NslawType = interactionxml->getNonSmoothLawXML()->getType();
  // ComplementarityConditionNSL
  if (NslawType == COMPLEMENTARITY_CONDITION_NSLAW_TAG)
    nslaw = new ComplementarityConditionNSL(interactionxml->getNonSmoothLawXML());
  // RelayNSL
  else if (NslawType == RELAY_NSLAW_TAG)
    nslaw = new RelayNSL(interactionxml->getNonSmoothLawXML());
  // NewtonImpactNSL
  else if (NslawType == NEWTON_IMPACT_NSLAW_TAG)
    nslaw = new NewtonImpactNSL(interactionxml->getNonSmoothLawXML());
  // Newton impact friction law
  else if (NslawType == NEWTON_IMPACT_FRICTION_NSLAW_TAG)
    nslaw = new NewtonImpactFrictionNSL(interactionxml->getNonSmoothLawXML());
  else RuntimeException::selfThrow("Interaction::xml constructor, unknown NSLAW type :" + nslaw->getType());
  isNsLawAllocatedIn = true;

  // --- Relation ---
  string relationType = interactionxml->getRelationXML()->getType();
  // general relation
  if (relationType == RELATION_TAG)
    relation = new Relation(interactionxml->getRelationXML(), this);

  // Linear relation
  else if (relationType == LINEAR_TIME_INVARIANT_RELATION_TAG)
    relation = new LinearTIR(interactionxml->getRelationXML(), this);

  // Lagrangian non-linear relation
  else if (relationType == LAGRANGIAN_RELATION_TAG)
    relation = new LagrangianR(interactionxml->getRelationXML(), this);

  // Lagrangian linear relation
  else if (relationType == LAGRANGIAN_LINEAR_RELATION_TAG)
    relation = new LagrangianLinearR(interactionxml->getRelationXML(), this);
  else RuntimeException::selfThrow("Interaction::xml constructor, unknown relation type " + relation->getType());
  isRelationAllocatedIn = true;

  // check coherence between interactionSize and nsLawSize
  if ((interactionSize % nslaw->getNsLawSize()) != 0)
    RuntimeException::selfThrow("Interaction::xml constructor, inconsistency between interaction size and non smooth law size.");
  // compute number of relations.
  numberOfRelations = interactionSize / nslaw->getNsLawSize();
}

// --- Constructor from a set of data ---

Interaction::Interaction(const string& newId, DSSet& dsConcerned, const int& newNumber, const int& nInter):
  id(newId), number(newNumber), interactionSize(nInter), numberOfRelations(1), sizeOfDS(0), involvedDS(),
  nslaw(NULL), relation(NULL), NSDS(NULL), interactionxml(NULL), isRelationAllocatedIn(false), isNsLawAllocatedIn(false)
{
  // Memory allocation and initialization for y and lambda
  initialize();
  involvedDS = dsConcerned; // !! this keeps pointers link between DS in the set !!
  computeSizeOfDS();

  cout << "Warning, Interaction(id,number, ...) constructor: neither relation nor nslaw has been set. This should be done later." << endl;
  // Remark(FP): neither nslaw nor relation are created in this constructor -> todo?
  // this means that numberOfRelations is not computed -> set to one by default.
}

// --- DESTRUCTOR ---
Interaction::~Interaction()
{
  for (unsigned int i = 0; i < y.size(); i++)
  {
    if (isYAllocatedIn[i]) delete y[i];
    y[i] = NULL;
    if (isYOldAllocatedIn[i]) delete yOld[i];
    yOld[i] = NULL;
    if (isLambdaAllocatedIn[i]) delete lambda[i];
    lambda[i] = NULL;
    if (isLambdaOldAllocatedIn[i]) delete lambdaOld[i];
    lambdaOld[i] = NULL;
  }
  y.clear();
  yOld.clear();
  lambda.clear();
  lambdaOld.clear();

  if (isRelationAllocatedIn) delete relation;
  relation = NULL;
  if (isNsLawAllocatedIn) delete nslaw;
  nslaw = NULL;
  NSDS = NULL;
}

void Interaction::initialize()
{
  // Memory allocation for y and lambda

  // \todo : compute numberOfDerivatives using relative degree
  // for the moment, numberOfDerivatives = 2: we save y and yDot
  unsigned int numberOfDerivatives = 2;
  y.resize(numberOfDerivatives) ;
  yOld.resize(numberOfDerivatives);
  lambda.resize(numberOfDerivatives);
  lambdaOld.resize(numberOfDerivatives);
  for (unsigned int i = 0; i < numberOfDerivatives ; i++)
  {
    y[i] = new SimpleVector(interactionSize);
    yOld[i] = new SimpleVector(interactionSize);
    lambda[i] = new SimpleVector(interactionSize);
    lambdaOld[i] = new SimpleVector(interactionSize);
  }

  isYAllocatedIn.resize(numberOfDerivatives, true);
  isYOldAllocatedIn.resize(numberOfDerivatives, true);
  isLambdaAllocatedIn.resize(numberOfDerivatives, true);
  isLambdaOldAllocatedIn.resize(numberOfDerivatives, true);

  //  initializeVectors(y);
  //  initializeVectors(yOld);
  //  initializeVectors(lambda);
  //  initializeVectors(lambdaOld);

}

// vector<SimpleVector*> initialization
// Initialization is required to avoid uninit memory reading (-> insure tests: READ_UNINIT_MEM)
void Interaction::initializeVectors(vector<SimpleVector*> inputVector)
{
  vector<SimpleVector*>::iterator iter;
  for (iter = inputVector.begin(); iter != inputVector.end(); ++iter)
    (*iter)->zero();
}

// --- GETTERS/SETTERS ---

void Interaction::setY(const vectorOfSiconosVector& newVector)
{
  // clear y
  for (unsigned int i = 0; i < y.size(); i++)
  {
    if (isYAllocatedIn[i]) delete y[i];
    y[i] = NULL;
  }
  y.clear();
  unsigned int size = newVector.size();
  y.resize(size, NULL);

  for (unsigned int i = 0; i < size; i++)
    y[i] = new SimpleVector(*(newVector[i])); // -> copy !
  isYAllocatedIn.resize(size, true);
}

void Interaction::setYPtr(const vectorOfSiconosVector& newVector)
{
  // clear y
  for (unsigned int i = 0; i < y.size(); i++)
  {
    if (isYAllocatedIn[i]) delete y[i];
    y[i] = NULL;
  }
  y.clear();

  // copy
  y = newVector; // warning: pointer equality between y[i] and newVector[i]
  isYAllocatedIn.resize(y.size(), false);
}

void Interaction::setY(const unsigned int & index, const SimpleVector& newY)
{
  if (y.size() <= index)
    RuntimeException::selfThrow("Interaction::setY, index out of range ");

  // set y[index]
  if (y[index] == NULL)
  {
    y[index] = new SimpleVector(newY);
    isYAllocatedIn[index] = true ;
  }
  else
  {
    if (y[index]->size() != newY.size())
      RuntimeException::selfThrow("Interaction::setY(index,newY), inconsistent sizes between y(index) and newY ");
    *(y[index]) = newY;
  }
}

void Interaction::setYPtr(const unsigned int & index, SimpleVector* newY)
{
  if (y.size() <= index)
    RuntimeException::selfThrow("Interaction::setYPtr, index out of range ");
  if (newY->size() != interactionSize)
    RuntimeException::selfThrow("Interaction::setYPtr, interactionSize differs from newY vector size ");

  // set y[index]
  if (isYAllocatedIn[index]) delete y[index];
  y[index] = newY;
  isYAllocatedIn[index] = false ;
}

void Interaction::setYOld(const vectorOfSiconosVector& newVector)
{
  // clear yOld
  for (unsigned int i = 0; i < yOld.size(); i++)
  {
    if (isYOldAllocatedIn[i]) delete yOld[i];
    yOld[i] = NULL;
  }
  yOld.clear();
  unsigned int size = newVector.size();
  yOld.resize(size, NULL);

  for (unsigned int i = 0; i < size; i++)
    yOld[i] = new SimpleVector(*(newVector[i])); // -> copy !
  isYOldAllocatedIn.resize(size, true);
}

void Interaction::setYOldPtr(const vectorOfSiconosVector& newVector)
{
  // clear yOld
  for (unsigned int i = 0; i < yOld.size(); i++)
  {
    if (isYOldAllocatedIn[i]) delete yOld[i];
    yOld[i] = NULL;
  }
  yOld.clear();

  // copy
  yOld = newVector; // warning: pointer equalityOld between yOld[i] and newVector[i]
  isYOldAllocatedIn.resize(yOld.size(), false);
}

void Interaction::setYOld(const unsigned int & index, const SimpleVector& newYOld)
{
  if (yOld.size() <= index)
    RuntimeException::selfThrow("Interaction::setYOld, index out of range ");

  // set yOld[index]
  if (yOld[index] == NULL)
  {
    yOld[index] = new SimpleVector(newYOld);
    isYOldAllocatedIn[index] = true ;
  }
  else
  {
    if (yOld[index]->size() != newYOld.size())
      RuntimeException::selfThrow("Interaction::setYOld(index,newYOld), inconsistent sizes between yOld(index) and newYOld ");
    *(yOld[index]) = newYOld;
  }
}

void Interaction::setYOldPtr(const unsigned int & index, SimpleVector* newYOld)
{
  if (yOld.size() <= index)
    RuntimeException::selfThrow("Interaction::setYOldPtr, index out of range ");
  if (newYOld->size() != interactionSize)
    RuntimeException::selfThrow("Interaction::setYOldPtr, interactionSize differs from newYOld vector size ");

  // set yOld[index]
  if (isYOldAllocatedIn[index]) delete yOld[index];
  yOld[index] = newYOld;
  isYOldAllocatedIn[index] = false ;
}

void Interaction::setLambda(const vectorOfSiconosVector& newVector)
{
  // clear lambda
  for (unsigned int i = 0; i < lambda.size(); i++)
  {
    if (isLambdaAllocatedIn[i]) delete lambda[i];
    lambda[i] = NULL;
  }
  lambda.clear();
  unsigned int size = newVector.size();
  lambda.resize(size, NULL);

  for (unsigned int i = 0; i < size; i++)
    lambda[i] = new SimpleVector(*(newVector[i])); // -> copy !
  isLambdaAllocatedIn.resize(size, true);
}

void Interaction::setLambdaPtr(const vectorOfSiconosVector& newVector)
{
  // clear lambda
  for (unsigned int i = 0; i < lambda.size(); i++)
  {
    if (isLambdaAllocatedIn[i]) delete lambda[i];
    lambda[i] = NULL;
  }
  lambda.clear();

  // copy
  lambda = newVector; // warning: pointer equality between lambda[i] and newVector[i]
  isLambdaAllocatedIn.resize(lambda.size(), false);
}

void Interaction::setLambda(const unsigned int & index, const SimpleVector& newLambda)
{
  if (lambda.size() <= index)
    RuntimeException::selfThrow("Interaction::setLambda, index out of range ");

  // set lambda[index]
  if (lambda[index] == NULL)
  {
    lambda[index] = new SimpleVector(newLambda);
    isLambdaAllocatedIn[index] = true ;
  }
  else
  {
    if (lambda[index]->size() != newLambda.size())
      RuntimeException::selfThrow("Interaction::setLambda(index,newLambda), inconsistent sizes between lambda(index) and newLambda ");
    *(lambda[index]) = newLambda;
  }
}

void Interaction::setLambdaPtr(const unsigned int & index, SimpleVector* newLambda)
{
  if (lambda.size() <= index)
    RuntimeException::selfThrow("Interaction::setLambdaPtr, index out of range ");
  if (newLambda->size() != interactionSize)
    RuntimeException::selfThrow("Interaction::setLambdaPtr, interactionSize differs from newLambda vector size ");

  // set lambda[index]
  if (isLambdaAllocatedIn[index]) delete lambda[index];
  lambda[index] = newLambda;
  isLambdaAllocatedIn[index] = false ;
}

void Interaction::setLambdaOld(const vectorOfSiconosVector& newVector)
{
  // clear lambdaOld
  for (unsigned int i = 0; i < lambdaOld.size(); i++)
  {
    if (isLambdaOldAllocatedIn[i]) delete lambdaOld[i];
    lambdaOld[i] = NULL;
  }
  lambdaOld.clear();
  unsigned int size = newVector.size();
  lambdaOld.resize(size, NULL);

  for (unsigned int i = 0; i < size; i++)
    lambdaOld[i] = new SimpleVector(*(newVector[i])); // -> copy !
  isLambdaOldAllocatedIn.resize(size, true);
}

void Interaction::setLambdaOldPtr(const vectorOfSiconosVector& newVector)
{
  // clear lambdaOld
  for (unsigned int i = 0; i < lambdaOld.size(); i++)
  {
    if (isLambdaOldAllocatedIn[i]) delete lambdaOld[i];
    lambdaOld[i] = NULL;
  }
  lambdaOld.clear();

  // copy
  lambdaOld = newVector; // warning: pointer equality between lambdaOld[i] and newVector[i]
  isLambdaOldAllocatedIn.resize(lambdaOld.size(), false);
}

void Interaction::setLambdaOld(const unsigned int & index, const SimpleVector& newLambdaOld)
{
  if (lambdaOld.size() <= index)
    RuntimeException::selfThrow("Interaction::setLambdaOld, index out of range ");

  // set lambdaOld[index]
  if (lambdaOld[index] == NULL)
  {
    lambdaOld[index] = new SimpleVector(newLambdaOld);
    isLambdaOldAllocatedIn[index] = true ;
  }
  else
  {
    if (lambdaOld[index]->size() != newLambdaOld.size())
      RuntimeException::selfThrow("Interaction::setLambdaOld(index,newLambdaOld), inconsistent sizes between lambdaOld(index) and newLambdaOld ");
    *(lambdaOld[index]) = newLambdaOld;
  }
}

void Interaction::setLambdaOldPtr(const unsigned int & index, SimpleVector* newLambdaOld)
{
  if (lambdaOld.size() <= index)
    RuntimeException::selfThrow("Interaction::setLambdaOldPtr, index out of range ");
  if (newLambdaOld->size() != interactionSize)
    RuntimeException::selfThrow("Interaction::setLambdaOldPtr, interactionSize differs from newLambdaOld vector size ");

  // set lambdaOld[index]
  if (isLambdaOldAllocatedIn[index]) delete lambdaOld[index];
  lambdaOld[index] = newLambdaOld;
  isLambdaOldAllocatedIn[index] = false ;
}


void Interaction::setDynamicalSystems(const DSSet& newSet)
{
  involvedDS = newSet; // !! Pointers links between ds !!
  computeSizeOfDS();
}

DynamicalSystem* Interaction::getDynamicalSystemPtr(const int& nb)
{
  if (! involvedDS.isDSIn(nb)) // if ds number nb is not in the set ...
    RuntimeException::selfThrow("Interaction::getDynamicalSystemPtr(nb), DS number nb is not in the set.");
  return involvedDS.getDynamicalSystem(number);
}

DynamicalSystem Interaction::getDynamicalSystem(const int& nb)
{
  if (! involvedDS.isDSIn(nb)) // if ds number nb is not in the set ...
    RuntimeException::selfThrow("Interaction::getDynamicalSystem(nb), DS number nb is not in the set.");
  return *(involvedDS.getDynamicalSystem(nb));
}

void Interaction::setRelationPtr(Relation* newRelation)
{
  if (isRelationAllocatedIn) delete relation;
  relation = newRelation;
  isRelationAllocatedIn = false;
  newRelation->setInteractionPtr(this);
}

void Interaction::setNonSmoothLawPtr(NonSmoothLaw* newNslaw)
{
  if (isNsLawAllocatedIn) delete nslaw;
  nslaw = newNslaw;
  isNsLawAllocatedIn = false;
}

// --- OTHER FUNCTIONS ---

void Interaction::computeSizeOfDS()
{
  sizeOfDS = 0;
  DSIterator it;
  for (it = involvedDS.begin(); it != involvedDS.end(); it++)
    sizeOfDS += (*it)->getDim();
}

void Interaction::swapInMemory()
{
  for (unsigned int i = 0; i < y.size() ; i++)
  {
    *(yOld[i]) = *(y[i]) ;
    *(lambdaOld[i]) = *(lambda[i]);
    lambda[i]->zero();
  }
}

void Interaction::display() const
{
  cout << "======= Interaction display =======" << endl;
  cout << "| id : " << id << endl;
  cout << "| number : " << number << endl;
  involvedDS.display();
  cout << "| y : " << endl;
  if (y[0] != NULL) y[0]->display();
  else cout << "->NULL" << endl;
  cout << "| yDot : " << endl;
  if (y[1] != NULL) y[1]->display();
  else cout << "->NULL" << endl;
  cout << "| yOld : " << endl;
  if (yOld[0] != NULL) yOld[0]->display();
  else cout << "->NULL" << endl;
  cout << "| yDotOld : " << endl;
  if (yOld[1] != NULL) yOld[1]->display();
  else cout << "->NULL" << endl;
  cout << "| lambda : " << endl;
  if (lambda[0] != NULL) lambda[0]->display();
  else cout << "->NULL" << endl;
  cout << "| lambdaDot : " << endl;
  if (lambda[1] != NULL) lambda[1]->display();
  else cout << "->NULL" << endl;
  cout << "===================================" << endl;
}

NonSmoothLaw* Interaction::createComplementarityConditionNSL()
{
  if (isNsLawAllocatedIn) delete nslaw;
  nslaw = new ComplementarityConditionNSL();
  isNsLawAllocatedIn = true;
  return nslaw;
}

NonSmoothLaw* Interaction::createRelayNSL(const double& c, const double& d)
{
  if (isNsLawAllocatedIn) delete nslaw;
  nslaw = new RelayNSL(c, d);
  isNsLawAllocatedIn = true;
  return nslaw;
}

NonSmoothLaw* Interaction::createNewtonImpactNSL(const double& e)
{
  if (isNsLawAllocatedIn) delete nslaw;
  nslaw = new NewtonImpactNSL(e);
  isNsLawAllocatedIn = true;
  return nslaw;
}

NonSmoothLaw* Interaction::createNewtonImpactFrictionNSL(const double& en, const double& et, const double& mu)
{
  if (isNsLawAllocatedIn) delete nslaw;
  nslaw = new NewtonImpactFrictionNSL(en, et, mu);
  isNsLawAllocatedIn = true;
  return nslaw;
}

// --- XML RELATED FUNCTIONS ---

void Interaction::saveInteractionToXML()
{
  /*
   * save the data of the Interaction
   */

  if (interactionxml != NULL)
  {
    //  interactionxml->setDSConcerned( involvedDS );
    interactionxml->setId(id);
    interactionxml->setNumber(number);
    interactionxml->setSize(interactionSize);
    interactionxml->setY(*(y[0]));
    interactionxml->setLambda(*(lambda[0]));
  }
  else RuntimeException::selfThrow("Interaction::saveInteractionToXML - object InteractionXML does not exist");

  /*
   * save the data of the Relation
   */
  if (relation->getType() == RELATION)
    relation->saveRelationToXML();
  else  if (relation->getType() == LINEARTIRELATION)
    (static_cast<LinearTIR*>(relation))->saveRelationToXML();
  else if (relation->getType() == LAGRANGIANLINEARRELATION)
    (static_cast<LagrangianLinearR*>(relation))->saveRelationToXML();
  else if (relation->getType() == LAGRANGIANRELATION)
    (static_cast<LagrangianR*>(relation))->saveRelationToXML();
  else RuntimeException::selfThrow("Interaction::saveInteractionToXML - bad kind of Relation :" + relation->getType());
  /*
   * save the data of the NonSmoothLaw
   */

  if (nslaw->getType() == COMPLEMENTARITYCONDITIONNSLAW)
    (static_cast<ComplementarityConditionNSL*>(nslaw))->saveNonSmoothLawToXML();
  else if (nslaw->getType() == RELAYNSLAW)
    (static_cast<RelayNSL*>(nslaw))->saveNonSmoothLawToXML();
  else if (nslaw->getType() == NEWTONIMPACTNSLAW)
    (static_cast<NewtonImpactNSL*>(nslaw))->saveNonSmoothLawToXML();
  else if (nslaw->getType() == NEWTONIMPACTFRICTIONNSLAW)
    (static_cast<NewtonImpactFrictionNSL*>(nslaw))->saveNonSmoothLawToXML();
  else RuntimeException::selfThrow("Interaction::saveInteractionToXML - bad kind of NonSmoothLaw : " + nslaw->getType());
}

// Default (private) constructor
Interaction::Interaction():
  id("none"), number(0), interactionSize(0), sizeOfDS(0), nslaw(NULL), relation(NULL), interactionxml(NULL),
  isRelationAllocatedIn(false), isNsLawAllocatedIn(false)
{}
