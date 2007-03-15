/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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
#include "InteractionXML.h"
#include "NonSmoothLawXML.h"
#include "RelationXML.h"
#include "FirstOrderLinearTIR.h"
#include "LagrangianLinearR.h"
#include "LagrangianScleronomousR.h"
#include "LagrangianRheonomousR.h"
#include "LagrangianCompliantR.h"
#include "ComplementarityConditionNSL.h"
#include "RelayNSL.h"
#include "NewtonImpactNSL.h"
#include "NewtonImpactFrictionNSL.h"
#include "NonSmoothDynamicalSystem.h"
#include "DynamicalSystem.h"

using namespace std;

// --- CONSTRUCTORS ---

// Default (private) constructor
Interaction::Interaction():
  id("none"), number(0), interactionSize(0), sizeOfDS(0), sizeZ(0), nslaw(NULL), relation(NULL), interactionxml(NULL)
{}

// --- XML constructor ---
Interaction::Interaction(InteractionXML* interxml, NonSmoothDynamicalSystem * nsds):
  id("undefined"), number(0), interactionSize(0), numberOfRelations(0), sizeOfDS(0), sizeZ(0),
  nslaw(NULL), relation(NULL), NSDS(nsds), interactionxml(interxml)
{
  if (interactionxml == NULL)
    RuntimeException::selfThrow("Interaction::xml constructor, xmlfile = NULL");

  // id and number
  if (interactionxml->hasId()) id = interactionxml->getId();
  number = interactionxml->getNumber();

  // interaction size
  interactionSize = interactionxml->getSize();

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
  isAllocatedIn["nsLaw"] = true;

  // --- Dynamical Systems ---
  unsigned int sizeDS ;
  if (nsds != NULL)
  {
    // Get a list of DS concerned from xml

    if (interactionxml->hasAllDS())
      involvedDS = nsds->getDynamicalSystems();

    else
    {
      // get numbers of DS involved in the interaction from the xml input file.
      std::vector<int> listDS;
      interactionxml->getDSNumbers(listDS);

      // get corresponding DS and insert them into the involvedDS set.
      sizeDS = listDS.size();
      for (unsigned int i = 0; i < sizeDS; i++)
        involvedDS.insert(nsds->getDynamicalSystemPtrNumber(listDS[i]));
    }
  }
  else cout << "Interaction constructor, warning: no dynamical systems linked to the interaction!" << endl;

  // --- Relation ---
  string relationType = interactionxml->getRelationXML()->getType();

  // First Order Non Linear Relation
  if (relationType == "FirstOrderRelation")
    relation = new FirstOrderR(interactionxml->getRelationXML());

  // Linear relation
  else if (relationType == "FirstOrderLinearTimeInvariantRelation")
    relation = new FirstOrderLinearTIR(interactionxml->getRelationXML());

  // Lagrangian non-linear relation
  else if (relationType == "LagrangianRelation")
  {
    string relationSubType = interactionxml->getRelationXML()->getSubType();
    // \todo create a factory to avoid "if" list for Relation construction according to subType.
    if (relationSubType == "Scleronomous")
      relation = new LagrangianScleronomousR(interactionxml->getRelationXML());
    else if (relationSubType == "Rheonomous")
      relation = new LagrangianRheonomousR(interactionxml->getRelationXML());
    else if (relationSubType == "Compliant")
      relation = new LagrangianCompliantR(interactionxml->getRelationXML());
  }
  // Lagrangian linear relation
  else if (relationType == "LagrangianLinearRelation")
    relation = new LagrangianLinearR(interactionxml->getRelationXML());
  else RuntimeException::selfThrow("Interaction::xml constructor, unknown relation type " + relation->getType());
  isAllocatedIn["relation"] = true;

  // check coherence between interactionSize and nsLawSize
  if ((interactionSize % nslaw->getNsLawSize()) != 0)
    RuntimeException::selfThrow("Interaction::xml constructor, inconsistency between interaction size and non smooth law size.");

  // download xml values for y[0] and lambda[0] forbidden since y and lambda are not allocated yet (this is
  // done in initialize function, called by simulation !!)
  if (interactionxml->hasY() ||  interactionxml->hasLambda())
    RuntimeException::selfThrow("Interaction::xml constructor, y or lambda download is forbidden.");
}

// --- Constructor from a set of data ---

Interaction::Interaction(const string& newId, DynamicalSystemsSet& dsConcerned, int newNumber, int nInter, NonSmoothLaw* newNSL, Relation* newRel):
  id(newId), number(newNumber), interactionSize(nInter), numberOfRelations(1), sizeOfDS(0), sizeZ(0), involvedDS(),
  nslaw(newNSL), relation(newRel), NSDS(NULL), interactionxml(NULL)
{
  // Memory allocation and initialization for y and lambda
  involvedDS = dsConcerned; // !! this keeps pointers link between DS in the set !!
  isAllocatedIn["relation"] = false;
  isAllocatedIn["nsLaw"] = false;
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

  if (isAllocatedIn["relation"]) delete relation;
  relation = NULL;
  if (isAllocatedIn["nsLaw"]) delete nslaw;
  nslaw = NULL;
  NSDS = NULL;
}

void Interaction::initialize(double t0, unsigned int level)
{
  if (relation == NULL)
    RuntimeException::selfThrow("Interaction::initialize failed, relation == NULL");
  if (nslaw == NULL)
    RuntimeException::selfThrow("Interaction::initialize failed, non smooth law == NULL");

  computeSizeOfDS();

  relation->setInteractionPtr(this);
  relation->initialize();

  // compute number of relations.
  numberOfRelations = interactionSize / nslaw->getNsLawSize();

  initializeMemory(level);

  // Compute y values for t0
  for (unsigned int i = 0; i < level; ++i)
  {
    computeOutput(t0, i);
    //      computeInput(t0,i);
  }
}

// Initialize and InitializeMemory are separated in two functions since we need to know the relative degree to know "numberOfDerivatives",
// while numberOfRelations and the size of the non smooth law are required inputs to compute the relative degree.
void Interaction::initializeMemory(unsigned int numberOfDerivatives)
{
  // Warning: this function is called from Simulation initialize, since we need to know the relative degree and
  // the type of simulation to size Y and Lambda.

  // Memory allocation for y and lambda

  // Note that numberOfDerivatives depends on the type of simulation and on the relative degree.

  y.resize(numberOfDerivatives) ;
  yOld.resize(numberOfDerivatives);
  lambda.resize(numberOfDerivatives);
  lambdaOld.resize(numberOfDerivatives);

  // get the dimension of the non smooth law, ie the size of a unitary blocks (one per relation)
  unsigned int nslawSize = nslaw->getNsLawSize();

  for (unsigned int i = 0; i < numberOfDerivatives ; i++)
  {

    y[i] = new BlockVector();
    yOld[i] = new BlockVector();
    lambda[i] = new BlockVector();
    lambdaOld[i] = new BlockVector();
    for (unsigned int j = 0; j < numberOfRelations; ++j)
    {
      y[i]->addPtr(new SimpleVector(nslawSize));
      yOld[i]->addPtr(new SimpleVector(nslawSize));
      lambda[i]->addPtr(new SimpleVector(nslawSize));
      lambdaOld[i]->addPtr(new SimpleVector(nslawSize));
    }
  }

  isYAllocatedIn.resize(numberOfDerivatives, true);
  isYOldAllocatedIn.resize(numberOfDerivatives, true);
  isLambdaAllocatedIn.resize(numberOfDerivatives, true);
  isLambdaOldAllocatedIn.resize(numberOfDerivatives, true);
}

// --- GETTERS/SETTERS ---

void Interaction::setY(const VectorOfVectors& newVector)
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
    y[i] = new BlockVector(*(newVector[i])); // -> copy !
  isYAllocatedIn.resize(size, true);
}

void Interaction::setYPtr(const VectorOfVectors& newVector)
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

void Interaction::setY(const unsigned int  index, const BlockVector& newY)
{
  if (y.size() <= index)
    RuntimeException::selfThrow("Interaction::setY, index out of range ");

  // set y[index]
  if (y[index] == NULL)
  {
    y[index] = new BlockVector(newY);
    isYAllocatedIn[index] = true ;
  }
  else
  {
    if (y[index]->size() != newY.size())
      RuntimeException::selfThrow("Interaction::setY(index,newY), inconsistent sizes between y(index) and newY ");
    *(y[index]) = newY;
  }
}

void Interaction::setYPtr(const unsigned int  index, SiconosVector* newY)
{
  if (y.size() <= index)
    RuntimeException::selfThrow("Interaction::setYPtr, index out of range ");
  if (newY->size() != interactionSize)
    RuntimeException::selfThrow("Interaction::setYPtr, interactionSize differs from newY vector size ");
  if (!newY->isBlock())
    RuntimeException::selfThrow("Interaction::setYPtr(newY), newY is not a block vector! ");

  // set y[index]
  if (isYAllocatedIn[index]) delete y[index];

  y[index] = static_cast<BlockVector*>(newY);
  isYAllocatedIn[index] = false ;
}

void Interaction::setYOld(const VectorOfVectors& newVector)
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
    yOld[i] = new BlockVector(*(newVector[i])); // -> copy !
  isYOldAllocatedIn.resize(size, true);
}

void Interaction::setYOldPtr(const VectorOfVectors& newVector)
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

void Interaction::setYOld(const unsigned int  index, const BlockVector& newYOld)
{
  if (yOld.size() <= index)
    RuntimeException::selfThrow("Interaction::setYOld, index out of range ");

  // set yOld[index]
  if (yOld[index] == NULL)
  {
    yOld[index] = new BlockVector(newYOld);
    isYOldAllocatedIn[index] = true ;
  }
  else
  {
    if (yOld[index]->size() != newYOld.size())
      RuntimeException::selfThrow("Interaction::setYOld(index,newYOld), inconsistent sizes between yOld(index) and newYOld ");
    *(yOld[index]) = newYOld;
  }
}

void Interaction::setYOldPtr(const unsigned int  index, SiconosVector* newYOld)
{
  if (yOld.size() <= index)
    RuntimeException::selfThrow("Interaction::setYOldPtr, index out of range ");
  if (newYOld->size() != interactionSize)
    RuntimeException::selfThrow("Interaction::setYOldPtr, interactionSize differs from newYOld vector size ");
  if (!newYOld->isBlock())
    RuntimeException::selfThrow("Interaction::setYOldPtr(newY), newY is not a block vector! ");

  // set yOld[index]
  if (isYOldAllocatedIn[index]) delete yOld[index];
  yOld[index] = static_cast<BlockVector*>(newYOld);
  isYOldAllocatedIn[index] = false ;
}

void Interaction::setLambda(const VectorOfVectors& newVector)
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
    lambda[i] = new BlockVector(*(newVector[i])); // -> copy !
  isLambdaAllocatedIn.resize(size, true);
}

void Interaction::setLambdaPtr(const VectorOfVectors& newVector)
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

void Interaction::setLambda(const unsigned int  index, const BlockVector& newLambda)
{
  if (lambda.size() <= index)
    RuntimeException::selfThrow("Interaction::setLambda, index out of range ");

  // set lambda[index]
  if (lambda[index] == NULL)
  {
    lambda[index] = new BlockVector(newLambda);
    isLambdaAllocatedIn[index] = true ;
  }
  else
  {
    if (lambda[index]->size() != newLambda.size())
      RuntimeException::selfThrow("Interaction::setLambda(index,newLambda), inconsistent sizes between lambda(index) and newLambda ");
    *(lambda[index]) = newLambda;
  }
}

void Interaction::setLambdaPtr(const unsigned int  index, SiconosVector* newLambda)
{
  if (lambda.size() <= index)
    RuntimeException::selfThrow("Interaction::setLambdaPtr, index out of range ");
  if (newLambda->size() != interactionSize)
    RuntimeException::selfThrow("Interaction::setLambdaPtr, interactionSize differs from newLambda vector size ");
  if (!newLambda->isBlock())
    RuntimeException::selfThrow("Interaction::setLambdaPtr(newLambda), newLambda is not a block vector! ");

  // set lambda[index]
  if (isLambdaAllocatedIn[index]) delete lambda[index];
  lambda[index] = static_cast<BlockVector*>(newLambda);
  isLambdaAllocatedIn[index] = false ;
}

void Interaction::setLambdaOld(const VectorOfVectors& newVector)
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
    lambdaOld[i] = new BlockVector(*(newVector[i])); // -> copy !
  isLambdaOldAllocatedIn.resize(size, true);
}

void Interaction::setLambdaOldPtr(const VectorOfVectors& newVector)
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

void Interaction::setLambdaOld(const unsigned int  index, const BlockVector& newLambdaOld)
{
  if (lambdaOld.size() <= index)
    RuntimeException::selfThrow("Interaction::setLambdaOld, index out of range ");

  // set lambdaOld[index]
  if (lambdaOld[index] == NULL)
  {
    lambdaOld[index] = new BlockVector(newLambdaOld);
    isLambdaOldAllocatedIn[index] = true ;
  }
  else
  {
    if (lambdaOld[index]->size() != newLambdaOld.size())
      RuntimeException::selfThrow("Interaction::setLambdaOld(index,newLambdaOld), inconsistent sizes between lambdaOld(index) and newLambdaOld ");
    *(lambdaOld[index]) = newLambdaOld;
  }
}

void Interaction::setLambdaOldPtr(const unsigned int  index, SiconosVector* newLambdaOld)
{
  if (lambdaOld.size() <= index)
    RuntimeException::selfThrow("Interaction::setLambdaOldPtr, index out of range ");
  if (newLambdaOld->size() != interactionSize)
    RuntimeException::selfThrow("Interaction::setLambdaOldPtr, interactionSize differs from newLambdaOld vector size ");
  if (!newLambdaOld->isBlock())
    RuntimeException::selfThrow("Interaction::setLambdaOldPtr(newLambda), newLambda is not a block vector! ");

  // set lambdaOld[index]
  if (isLambdaOldAllocatedIn[index]) delete lambdaOld[index];
  lambdaOld[index] = static_cast<BlockVector*>(newLambdaOld);
  isLambdaOldAllocatedIn[index] = false ;
}


void Interaction::setDynamicalSystems(const DynamicalSystemsSet& newSet)
{
  involvedDS = newSet; // !! Pointers links between ds !!
  computeSizeOfDS();
}

DynamicalSystem* Interaction::getDynamicalSystemPtr(int nb)
{
  if (! involvedDS.isDynamicalSystemIn(nb)) // if ds number nb is not in the set ...
    RuntimeException::selfThrow("Interaction::getDynamicalSystemPtr(nb), DS number nb is not in the set.");
  return involvedDS.getDynamicalSystemPtr(number);
}

void Interaction::getDynamicalSystem(int nb, DynamicalSystem& ds)
{
  // This function is useless in C++ but maybe required in Python? To be checked.
  // DS is a parameter, since it can be returned, DynamicalSystem being an abstract class.
  if (! involvedDS.isDynamicalSystemIn(nb)) // if ds number nb is not in the set ...
    RuntimeException::selfThrow("Interaction::getDynamicalSystem(nb), DS number nb is not in the set.");
  ds = *(involvedDS.getDynamicalSystemPtr(nb));
}

void Interaction::setRelationPtr(Relation* newRelation)
{
  if (isAllocatedIn["relation"]) delete relation;
  relation = newRelation;
  isAllocatedIn["relation"] = false;
}

void Interaction::setNonSmoothLawPtr(NonSmoothLaw* newNslaw)
{
  if (isAllocatedIn["nsLaw"]) delete nslaw;
  nslaw = newNslaw;
  isAllocatedIn["nsLaw"] = false;
}

// --- OTHER FUNCTIONS ---

void Interaction::computeSizeOfDS()
{
  sizeOfDS = 0;
  sizeZ = 0;
  DSIterator it;
  for (it = involvedDS.begin(); it != involvedDS.end(); it++)
  {
    sizeOfDS += (*it)->getDim();
    sizeZ += (*it)->getZPtr()->size();
  }
}

void Interaction::swapInMemory()
{
  // i corresponds to the derivative number and j the relation number.
  for (unsigned int i = 0; i < y.size() ; i++)
  {
    for (unsigned int j = 0; j < numberOfRelations; ++j)
    {
      *(yOld[i]->getVectorPtr(j)) = *(y[i]->getVectorPtr(j)) ;
      *(lambdaOld[i]->getVectorPtr(j)) = *(lambda[i]->getVectorPtr(j));
    }
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

void Interaction::computeOutput(double time, unsigned int level)
{
  relation->computeOutput(time, level);
}

void Interaction::computeFreeOutput(double time, unsigned int level)
{
  relation->computeFreeOutput(time, level);
}

void Interaction::computeInput(double time, unsigned int level)
{
  relation->computeInput(time, level);
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
  // Main type of the relation: FirstOrder or Lagrangian
  string type = relation->getType();
  // Subtype of the relation
  string subType = relation->getSubType();

  if (type == "FirstOrderRelation")
  {
    if (subType == "NonLinearR")
      relation->saveRelationToXML();
    else if (subType == "LinearTIR")
      (static_cast<FirstOrderLinearTIR*>(relation))->saveRelationToXML();
    else
      RuntimeException::selfThrow("Interaction::saveInteractionToXML - Unknown relation subtype: " + subType);
  }
  else if (type == "Lagrangian")
  {
    if (subType == "LinearR")
      (static_cast<LagrangianLinearR*>(relation))->saveRelationToXML();
    else
      RuntimeException::selfThrow("Interaction::saveInteractionToXML - Not yet implemented for relation subtype " + subType);
  }
  else
    RuntimeException::selfThrow("Interaction::saveInteractionToXML - Unknown relation type: " + type);

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

