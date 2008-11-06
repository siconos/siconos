/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
#include <assert.h>
#include "Interaction.h"
#include "InteractionXML.h"
#include "NonSmoothLawXML.h"
#include "RelationXML.h"
#include "RelationTypes.h"
#include "ComplementarityConditionNSL.h"
#include "RelayNSL.h"
#include "NewtonImpactNSL.h"
#include "NewtonImpactFrictionNSL.h"
#include "DynamicalSystem.h"

using namespace std;
using namespace RELATION;

// --- CONSTRUCTORS ---

// --- XML constructor ---
Interaction::Interaction(SP::InteractionXML interxml, SP::DynamicalSystemsSet nsdsSet):
  id("undefined"), number(0), interactionSize(0), numberOfRelations(0),
  sizeOfDS(0), sizeZ(0), interactionxml(interxml)
{
  assert(interactionxml && "NULL pointer");

  // id and number
  if (interactionxml->hasId()) id = interactionxml->getId();
  number = interactionxml->getNumber();

  // interaction size
  interactionSize = interactionxml->getSize();

  // --- Non smooth law ---
  string NslawType = interactionxml->getNonSmoothLawXML()->getType();
  // ComplementarityConditionNSL
  if (NslawType == COMPLEMENTARITY_CONDITION_NSLAW_TAG)
  {
    nslaw.reset(new ComplementarityConditionNSL(interactionxml->getNonSmoothLawXML()));
  }

  // RelayNSL
  else if (NslawType == RELAY_NSLAW_TAG)
  {
    nslaw.reset(new RelayNSL(interactionxml->getNonSmoothLawXML()));
  }

  // NewtonImpactNSL
  else if (NslawType == NEWTON_IMPACT_NSLAW_TAG)
  {
    nslaw.reset(new NewtonImpactNSL(interactionxml->getNonSmoothLawXML()));
  }
  // Newton impact friction law
  else if (NslawType == NEWTON_IMPACT_FRICTION_NSLAW_TAG)
  {
    nslaw.reset(new NewtonImpactFrictionNSL(interactionxml->getNonSmoothLawXML()));
  }
  else RuntimeException::selfThrow("Interaction::xml constructor, unknown NSLAW type :" + nslaw->getType());

  // --- Dynamical Systems ---
  unsigned int sizeDS ;
  involvedDS.reset(new DynamicalSystemsSet());
  if (!nsdsSet->isEmpty())
  {
    // Get a list of DS concerned from xml

    if (interactionxml->hasAllDS())
      involvedDS->insert(nsdsSet->begin(), nsdsSet->end());

    else
    {
      // get numbers of DS involved in the interaction from the xml input file.
      std::vector<int> dsNumbers;
      interactionxml->getDSNumbers(dsNumbers);
      // get corresponding DS and insert them into the involvedDS set.
      for (vector<int>::iterator it = dsNumbers.begin(); it != dsNumbers.end(); ++it)
        involvedDS->insert(nsdsSet->getPtr(*it));
    }
  }
  else cout << "Interaction constructor, warning: no dynamical systems linked to the interaction!" << endl;

  // --- Relation ---
  RELATION::TYPES relationType = interactionxml->getRelationXML()->getType();

  // First Order Non Linear Relation
  if (relationType == FirstOrder)
  {
    RELATION::SUBTYPES relationSubType = interactionxml->getRelationXML()->getSubType();
    if (relationSubType == Type1R)
      relation.reset(new FirstOrderType1R(interactionxml->getRelationXML()));
    // Linear relation
    else if (relationSubType == LinearR)
      relation.reset(new FirstOrderLinearR(interactionxml->getRelationXML()));
    // Linear time-invariant coef. relation
    else if (relationSubType == LinearTIR)
      relation.reset(new FirstOrderLinearTIR(interactionxml->getRelationXML()));
  }
  // Lagrangian non-linear relation
  else if (relationType == Lagrangian)
  {
    RELATION::SUBTYPES relationSubType = interactionxml->getRelationXML()->getSubType();
    // \todo create a factory to avoid "if" list for Relation construction according to subType.
    if (relationSubType == ScleronomousR)
      relation.reset(new LagrangianScleronomousR(interactionxml->getRelationXML()));
    else if (relationSubType == RheonomousR)
      relation.reset(new LagrangianRheonomousR(interactionxml->getRelationXML()));
    else if (relationSubType == CompliantR)
      relation.reset(new LagrangianCompliantR(interactionxml->getRelationXML()));
    // Lagrangian linear relation
    else if (relationSubType == LinearR || relationSubType == LinearTIR)
      relation.reset(new LagrangianLinearTIR(interactionxml->getRelationXML()));
  }
  else RuntimeException::selfThrow("Interaction::xml constructor, unknown relation type " + relation->getType());

  // check coherence between interactionSize and nsLawSize
  if ((interactionSize % nslaw->getNsLawSize()) != 0)
    RuntimeException::selfThrow("Interaction::xml constructor, inconsistency between interaction size and non smooth law size.");

  // download xml values for y[0] and lambda[0] forbidden since y and lambda are not allocated yet (this is
  // done in initialize function, called by simulation !!)
  if (interactionxml->hasY() ||  interactionxml->hasLambda())
    RuntimeException::selfThrow("Interaction::xml constructor, y or lambda download is forbidden.");
}

// --- Constructors from a set of data ---

Interaction::Interaction(SP::DynamicalSystem ds, int newNumber, int nInter,
                         SP::NonSmoothLaw newNSL, SP::Relation newRel):
  id("none"), number(newNumber), interactionSize(nInter), numberOfRelations(1),
  sizeOfDS(0), sizeZ(0), y(1), nslaw(newNSL), relation(newRel)
{
  involvedDS.reset(new DynamicalSystemsSet());
  involvedDS->insert(ds); // Warning: insert pointer to DS!!

}
Interaction::Interaction(const string& newId, SP::DynamicalSystem ds,
                         int newNumber, int nInter, SP::NonSmoothLaw newNSL, SP::Relation newRel):
  id(newId), number(newNumber), interactionSize(nInter), numberOfRelations(1),
  sizeOfDS(0), sizeZ(0), y(1), nslaw(newNSL), relation(newRel)
{
  involvedDS.reset(new DynamicalSystemsSet());
  involvedDS->insert(ds); // Warning: insert pointer to DS!!

}


Interaction::Interaction(DynamicalSystemsSet& dsConcerned, int newNumber, int nInter,
                         SP::NonSmoothLaw newNSL, SP::Relation newRel):
  id("none"), number(newNumber), interactionSize(nInter), numberOfRelations(1),
  sizeOfDS(0), sizeZ(0), y(1), nslaw(newNSL), relation(newRel)
{
  involvedDS.reset(new DynamicalSystemsSet());
  DSIterator itDS;
  for (itDS = dsConcerned.begin(); itDS != dsConcerned.end(); ++itDS)
    involvedDS->insert(*itDS); // Warning: insert pointers to DS!!
}

Interaction::Interaction(const string& newId, DynamicalSystemsSet& dsConcerned, int newNumber,
                         int nInter, SP::NonSmoothLaw newNSL, SP::Relation newRel):
  id(newId), number(newNumber), interactionSize(nInter), numberOfRelations(1), sizeOfDS(0), sizeZ(0),
  y(1), nslaw(newNSL), relation(newRel)
{
  involvedDS.reset(new DynamicalSystemsSet());
  DSIterator itDS;
  for (itDS = dsConcerned.begin(); itDS != dsConcerned.end(); ++itDS)
    involvedDS->insert(*itDS); // Warning: insert pointers to DS!!
}

// --- DESTRUCTOR ---
Interaction::~Interaction()
{
}

void Interaction::initialize(double t0, unsigned int level)
{
  assert(relation && "Interaction::initialize failed, relation == NULL");

  assert(nslaw && "Interaction::initialize failed, non smooth law == NULL");

  computeSizeOfDS();

  relation->setInteractionPtr(shared_from_this());
  relation->initialize(shared_from_this());

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

    y[i].reset(new BlockVector());
    yOld[i].reset(new BlockVector());
    lambda[i].reset(new BlockVector());
    lambdaOld[i].reset(new BlockVector());
    for (unsigned int j = 0; j < numberOfRelations; ++j)
    {
      y[i]->insertPtr(SP::SimpleVector(new SimpleVector(nslawSize)));
      yOld[i]->insertPtr(SP::SimpleVector(new SimpleVector(nslawSize)));
      lambda[i]->insertPtr(SP::SimpleVector(new SimpleVector(nslawSize)));
      lambdaOld[i]->insertPtr(SP::SimpleVector(new SimpleVector(nslawSize)));
    }
  }

}

// --- GETTERS/SETTERS ---

void Interaction::setY(const VectorOfVectors& newVector)
{

  unsigned int size = newVector.size();

  y.clear();
  y.resize(size);

  for (unsigned int i = 0; i < size; i++)
    y[i].reset(new BlockVector(*(newVector[i]))); // -> copy !
}

void Interaction::setYPtr(const VectorOfVectors& newVector)
{
  y.clear();

  // copy
  y = newVector; // smart ptr
}

void Interaction::setY(const unsigned int  index, const BlockVector& newY)
{
  assert(y.size() > index &&
         "Interaction::setY, index out of range ");

  // set y[index]
  if (! y[index])
  {
    y[index].reset(new BlockVector(newY));
  }
  else
  {
    assert(y[index]->size() == newY.size() &&
           "Interaction::setY(index,newY), inconsistent sizes between y(index) and newY ");
    *(y[index]) = newY;
  }
}

void Interaction::setYPtr(const unsigned int  index, SP::SiconosVector newY)
{
  assert(y.size() > index &&
         "Interaction::setYPtr, index out of range");

  assert(newY->size() == interactionSize &&
         "Interaction::setYPtr, interactionSize differs from newY vector size");

  assert(newY->isBlock() &&
         "Interaction::setYPtr(newY), newY is not a block vector!");

  y[index] = boost::static_pointer_cast<BlockVector>(newY);
}

void Interaction::setYOld(const VectorOfVectors& newVector)
{
  unsigned int size = newVector.size();
  yOld.clear();
  yOld.resize(size);

  for (unsigned int i = 0; i < size; i++)
    yOld[i].reset(new BlockVector(*(newVector[i]))); // -> copy !
}

void Interaction::setYOldPtr(const VectorOfVectors& newVector)
{
  // clear yOld

  yOld.clear();

  // copy
  yOld = newVector; // smart ptr
}

void Interaction::setYOld(const unsigned int  index, const BlockVector& newYOld)
{
  if (yOld.size() <= index)
    RuntimeException::selfThrow("Interaction::setYOld, index out of range ");

  // set yOld[index]
  if (! yOld[index])
  {
    yOld[index].reset(new BlockVector(newYOld));
  }
  else
  {
    assert(yOld[index]->size() == newYOld.size() &&
           "Interaction::setYOld(index,newYOld), inconsistent sizes between yOld(index) and newYOld");
    *(yOld[index]) = newYOld;
  }
}

void Interaction::setYOldPtr(const unsigned int  index, SP::SiconosVector newYOld)
{
  assert(yOld.size() > index &&
         "Interaction::setYOldPtr, index out of range");

  assert(newYOld->size() == interactionSize &&
         "Interaction::setYOldPtr, interactionSize differs from newYOld vector size");

  assert((! newYOld->isBlock()) &&
         "Interaction::setYOldPtr(newY), newY is not a block vector!");

  // set yOld[index]
  yOld[index] = boost::static_pointer_cast<BlockVector>(newYOld);
}

void Interaction::setLambda(const VectorOfVectors& newVector)
{
  unsigned int size = newVector.size();
  lambda.clear();
  lambda.resize(size);

  for (unsigned int i = 0; i < size; i++)
    lambda[i].reset(new BlockVector(*(newVector[i]))); // -> copy !
}

void Interaction::setLambdaPtr(const VectorOfVectors& newVector)
{
  lambda.clear();

  lambda = newVector; // smart ptr
}

void Interaction::setLambda(const unsigned int  index, const BlockVector& newLambda)
{
  assert(lambda.size() <= index &&
         "Interaction::setLambda, index out of range");

  // set lambda[index]
  if (! lambda[index])
  {
    lambda[index].reset(new BlockVector(newLambda));
  }
  else
  {
    assert(lambda[index]->size() == newLambda.size() &&
           "Interaction::setLambda(index,newLambda), inconsistent sizes between lambda(index) and newLambda");
    *(lambda[index]) = newLambda;
  }
}

void Interaction::setLambdaPtr(const unsigned int  index, SP::SiconosVector newLambda)
{
  assert(lambda.size() > index &&
         "Interaction::setLambdaPtr, index out of range ");

  assert(newLambda->size() == interactionSize &&
         "Interaction::setLambdaPtr, interactionSize differs from newLambda vector size ");

  assert(newLambda->isBlock() &&
         "Interaction::setLambdaPtr(newLambda), newLambda is not a block vector! ");

  // set lambda[index]
  lambda[index] = boost::static_pointer_cast<BlockVector>(newLambda);
}

void Interaction::setLambdaOld(const VectorOfVectors& newVector)
{
  unsigned int size = newVector.size();

  // clear lambdaOld
  lambdaOld.clear();
  lambdaOld.resize(size);

  for (unsigned int i = 0; i < size; i++)
    lambdaOld[i].reset(new BlockVector(*(newVector[i]))); // -> copy !
}

void Interaction::setLambdaOldPtr(const VectorOfVectors& newVector)
{
  // clear lambdaOld
  lambdaOld.clear();

  // copy
  lambdaOld = newVector; // smart ptrs
}

void Interaction::setLambdaOld(const unsigned int  index, const BlockVector& newLambdaOld)
{
  assert(lambdaOld.size() > index &&
         "Interaction::setLambdaOld, index out of range ");

  // set lambdaOld[index]
  if (! lambdaOld[index])
  {
    lambdaOld[index].reset(new BlockVector(newLambdaOld));
  }
  else
  {
    if (lambdaOld[index]->size() != newLambdaOld.size())
      RuntimeException::selfThrow("Interaction::setLambdaOld(index,newLambdaOld), inconsistent sizes between lambdaOld(index) and newLambdaOld ");
    *(lambdaOld[index]) = newLambdaOld;
  }
}

void Interaction::setLambdaOldPtr(const unsigned int  index, SP::SiconosVector newLambdaOld)
{
  if (lambdaOld.size() > index)
    RuntimeException::selfThrow("Interaction::setLambdaOldPtr, index out of range ");
  if (newLambdaOld->size() != interactionSize)
    RuntimeException::selfThrow("Interaction::setLambdaOldPtr, interactionSize differs from newLambdaOld vector size ");
  if (!newLambdaOld->isBlock())
    RuntimeException::selfThrow("Interaction::setLambdaOldPtr(newLambda), newLambda is not a block vector! ");

  // set lambdaOld[index]
  lambdaOld[index] = boost::static_pointer_cast<BlockVector>(newLambdaOld);
}


void Interaction::setDynamicalSystems(const DynamicalSystemsSet& newSet)
{
  DSIterator itDS;
  for (itDS = newSet.begin(); itDS != newSet.end(); ++itDS)
    involvedDS->insert(*itDS); // smart ptrs

  computeSizeOfDS();
}

SP::DynamicalSystem Interaction::getDynamicalSystemPtr(int nb)
{
  assert(involvedDS->isIn(nb) &&  // if ds number nb is not in the set ...
         "Interaction::getDynamicalSystemPtr(nb), DS number nb is not in the set.");
  return involvedDS->getPtr(number);
}

void Interaction::setRelationPtr(SP::Relation newRelation)
{
  relation = newRelation;
}

void Interaction::setNonSmoothLawPtr(SP::NonSmoothLaw newNslaw)
{
  nslaw = newNslaw;
}

// --- OTHER FUNCTIONS ---

void Interaction::computeSizeOfDS()
{
  sizeOfDS = 0;
  sizeZ = 0;
  DSIterator it;
  SP::SiconosVector ZP;
  for (it = involvedDS->begin(); it != involvedDS->end(); it++)
  {
    sizeOfDS += (*it)->getDim();
    ZP = (*it)->getZPtr();
    if (ZP) sizeZ += ZP->size();
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
  }
}

void Interaction::display() const
{
  cout << "======= Interaction display =======" << endl;
  cout << "| id : " << id << endl;
  cout << "| number : " << number << endl;
  involvedDS->display();
  cout << "| y : " << endl;
  if (y[0]) y[0]->display();
  else cout << "->NULL" << endl;
  cout << "| yDot : " << endl;
  if (y[1]) y[1]->display();
  else cout << "->NULL" << endl;
  cout << "| yOld : " << endl;
  if (yOld[0]) yOld[0]->display();
  else cout << "->NULL" << endl;
  cout << "| yDotOld : " << endl;
  if (yOld[1]) yOld[1]->display();
  else cout << "->NULL" << endl;
  cout << "| lambda : " << endl;
  if (lambda[0]) lambda[0]->display();
  else cout << "->NULL" << endl;
  cout << "| lambdaDot : " << endl;
  if (lambda[1]) lambda[1]->display();
  else cout << "->NULL" << endl;
  cout << "===================================" << endl;
}

void Interaction::computeOutput(double time, unsigned int level)
{
  relation->computeOutput(time, level);
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

  if (interactionxml)
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
  RELATION::TYPES type = relation->getType();
  // Subtype of the relation
  RELATION::SUBTYPES subType = relation->getSubType();

  if (type == FirstOrder)
  {
    if (subType == NonLinearR)
      relation->saveRelationToXML();
    else if (subType == LinearR)
      (boost::static_pointer_cast<FirstOrderLinearR>(relation))->saveRelationToXML();
    else if (subType == LinearTIR)
      (boost::static_pointer_cast<FirstOrderLinearTIR>(relation))->saveRelationToXML();
    else
      RuntimeException::selfThrow("Interaction::saveInteractionToXML - Unknown relation subtype: " + subType);
  }
  else if (type == Lagrangian)
  {
    if (subType == LinearTIR)
      (boost::static_pointer_cast<LagrangianLinearTIR>(relation))->saveRelationToXML();
    else
      RuntimeException::selfThrow("Interaction::saveInteractionToXML - Not yet implemented for relation subtype " + subType);
  }
  else
    RuntimeException::selfThrow("Interaction::saveInteractionToXML - Unknown relation type: " + type);

  /*
   * save the data of the NonSmoothLaw
   */

  if (nslaw->getType() == COMPLEMENTARITYCONDITIONNSLAW)
  {
    (boost::static_pointer_cast<ComplementarityConditionNSL>(nslaw))->saveNonSmoothLawToXML();
  }
  else if (nslaw->getType() == RELAYNSLAW)
  {
    (boost::static_pointer_cast<RelayNSL>(nslaw))->saveNonSmoothLawToXML();
  }
  else if (nslaw->getType() == NEWTONIMPACTNSLAW)
  {
    (boost::static_pointer_cast<NewtonImpactNSL>(nslaw))->saveNonSmoothLawToXML();
  }

  else if (nslaw->getType() == NEWTONIMPACTFRICTIONNSLAW)
  {
    (boost::static_pointer_cast<NewtonImpactFrictionNSL>(nslaw))->saveNonSmoothLawToXML();
  }

  else RuntimeException::selfThrow("Interaction::saveInteractionToXML - bad kind of NonSmoothLaw : " + nslaw->getType());
}

