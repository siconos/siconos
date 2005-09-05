#include "Interaction.h"

// includes to be deleted thanks to factories
#include "LinearTIR.h"
#include "LagrangianLinearR.h"
#include "LagrangianNonLinearR.h"
#include "ComplementarityConditionNSL.h"
#include "RelayNSL.h"
#include "NewtonImpactLawNSL.h"
#include "NewtonImpactFrictionNSL.h"

using namespace std;

// --- CONSTRUCTORS ---

// Copy constructor
Interaction::Interaction(const Interaction& newI):
  id(newI.getId()), number(newI.getNumber()), nInteraction(newI.getNInteraction()),
  nslaw(NULL), relation(NULL), interactionxml(NULL),
  isRelationAllocatedIn(true), isNsLawAllocatedIn(true)
{
  // Memory allocation and copy for simple vectors
  unsigned int size = newI.getY().size();
  y.resize(size, NULL);
  yOld.resize(size, NULL);
  lambda.resize(size, NULL);
  lambdaOld.resize(size, NULL);

  for (unsigned int i = 0; i < size; i++)
  {
    y[i] = new SimpleVector(*(newI.getYPtr(i)));
    yOld[i] = new SimpleVector(*(newI.getYOldPtr(i)));
    lambda[i] = new SimpleVector(*(newI.getLambdaPtr(i)));
    lambdaOld[i] = new SimpleVector(*(newI.getLambdaOldPtr(i)));
  }

  isYAllocatedIn.resize(size, true);
  isYOldAllocatedIn.resize(size, true);
  isLambdaAllocatedIn.resize(size, true);
  isLambdaOldAllocatedIn.resize(size, true);

  vectorDS.clear();
  vectorDS = newI.getDynamicalSystems();

  // Nslaw (warning! nslaw is an abstract class)
  string NslawType = newI.getNonSmoothLawPtr()->getType();
  if (NslawType ==  COMPLEMENTARITYCONDITIONNSLAW)
    nslaw = new ComplementarityConditionNSL();
  else if (NslawType == NEWTONIMPACTLAWNSLAW)
    nslaw = new NewtonImpactLawNSL();
  else if (NslawType == NEWTONIMPACTFRICTIONNSLAW)
    nslaw = new NewtonImpactFrictionNSL();
  else if (NslawType == RELAYNSLAW)
    nslaw = new  RelayNSL();
  else RuntimeException::selfThrow("Interaction::copy constructor, unknown NSLAW type :" + nslaw->getType());

  *nslaw = *newI.getNonSmoothLawPtr();

  // Relation
  string relationType = newI.getRelationPtr()->getType();
  // -> call copy constructor
  if (relationType == LINEARTIRELATION)
    relation = new LinearTIR(*(newI.getRelationPtr()), this);

  else if (relationType == LAGRANGIANLINEARRELATION)
    relation = new LagrangianLinearR(*(newI.getRelationPtr()), this);

  else if (relationType == LAGRANGIANNONLINEARRELATION)
    relation = new LagrangianNonLinearR(*(newI.getRelationPtr()));

  else RuntimeException::selfThrow("Interaction::copy constructor, unknown relation type " + relation->getType());

  // \remark FP: we do not link xml object in the copy
}

// --- XML constructor ---
Interaction::Interaction(InteractionXML* interxml, NonSmoothDynamicalSystem * nsds):
  id("none"), number(0), nInteraction(0), nslaw(NULL), relation(NULL), interactionxml(interxml),
  isRelationAllocatedIn(true), isNsLawAllocatedIn(true)
{
  if (interactionxml != NULL)
  {
    if (interactionxml->hasId()) id = interactionxml->getId();
    number = interactionxml->getNumber();
    nInteraction = interactionxml->getNInter();
    // Memory allocation for simple vectors

    // \todo : compute size using relative degree
    // for the moment, size = 2: we save y and yDot
    unsigned int size = 2;
    y.resize(size) ;
    yOld.resize(size);
    lambda.resize(size) ;
    lambdaOld.resize(size);
    for (unsigned int i = 0; i < size ; i++)
    {
      y[i] = new SimpleVector(nInteraction);
      yOld[i] = new SimpleVector(nInteraction);
      lambda[i] = new SimpleVector(nInteraction);
      lambdaOld[i] = new SimpleVector(nInteraction);
    }

    isYAllocatedIn.resize(size, true);
    isYOldAllocatedIn.resize(size, true);
    isLambdaAllocatedIn.resize(size, true);
    isLambdaOldAllocatedIn.resize(size, true);

    if (interactionxml->hasY()) *(y[0]) = interactionxml->getY();
    if (interactionxml->hasLambda()) *(lambda[0]) = interactionxml->getLambda();

    // Old values are initialized with current values
    swapInMemory();

    // --- Dynamical Systems ---

    if (nsds != NULL)
    {
      // Get a list of DS concerned from xml
      if (interactionxml->hasAll())
        vectorDS = nsds->getDynamicalSystems();
      else
      {
        vector<int> listDS = interactionxml->getDSConcerned();
        unsigned int size = listDS.size();
        vectorDS.resize(size);
        for (unsigned int i = 0; i < size; i++)
          vectorDS[i] = nsds->getDynamicalSystemPtrNumber(listDS[i]);
      }
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
    // NewtonImpactLawNSL
    else if (NslawType == NEWTON_IMPACT_LAW_NSLAW_TAG)
      nslaw = new NewtonImpactLawNSL(interactionxml->getNonSmoothLawXML());
    // Newton impact friction law
    else if (NslawType == NEWTON_IMPACT_FRICTION_NSLAW_TAG)
      nslaw = new NewtonImpactFrictionNSL(interactionxml->getNonSmoothLawXML());
    else RuntimeException::selfThrow("Interaction::xml constructor, unknown NSLAW type :" + nslaw->getType());

    // --- Relation ---
    string relationType = interactionxml->getRelationXML()->getType();
    // Linear relation
    if (relationType == LINEAR_TIME_INVARIANT_RELATION_TAG)
      relation = new LinearTIR(interactionxml->getRelationXML(), this);

    // Lagrangian linear relation
    else if (relationType == LAGRANGIAN_LINEAR_RELATION_TAG)
      relation = new LagrangianLinearR(interactionxml->getRelationXML(), this);

    // Lagrangian non-linear relation
    else if (relationType == LAGRANGIAN_NON_LINEAR_RELATION_TAG)
      relation = new LagrangianNonLinearR(interactionxml->getRelationXML(), this);

    else RuntimeException::selfThrow("Interaction::xml constructor, unknown relation type " + relation->getType());
  }
  else RuntimeException::selfThrow("Interaction::xml constructor, xmlfile = NULL");
}

// --- Constructor from a set of data ---

Interaction::Interaction(const string& newId, const int& newNumber, const int& nInter,
                         vector<DynamicalSystem*> *dsConcerned):
  id(newId), number(newNumber), nInteraction(nInter), nslaw(NULL),
  relation(NULL), interactionxml(NULL), isRelationAllocatedIn(false), isNsLawAllocatedIn(false)
{
  // Memory allocation for simple vectors

  // \todo : compute size using relative degree
  // for the moment, size = 2: we save y and yDot
  unsigned int size = 2;
  y.resize(size) ;
  yOld.resize(size);
  lambda.resize(size);
  lambdaOld.resize(size);
  for (unsigned int i = 0; i < size ; i++)
  {
    y[i] = new SimpleVector(nInteraction);
    yOld[i] = new SimpleVector(nInteraction);
    lambda[i] = new SimpleVector(nInteraction);
    lambdaOld[i] = new SimpleVector(nInteraction);
  }

  isYAllocatedIn.resize(size, true);
  isYOldAllocatedIn.resize(size, true);
  isLambdaAllocatedIn.resize(size, true);
  isLambdaOldAllocatedIn.resize(size, true);

  vectorDS.clear();
  if (dsConcerned != NULL) vectorDS = *dsConcerned;
  else RuntimeException::selfThrow("Interaction::createInteraction - The dsConcerned are not given");

  // Remark(FP): neither nslaw nor relation are created in this constructor -> todo?
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
}

// --- GETTERS/SETTERS ---

void Interaction::setY(const std::vector<SimpleVector*>& newVector)
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

void Interaction::setYPtr(const std::vector<SimpleVector*>& newVector)
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
  if (newY->size() != nInteraction)
    RuntimeException::selfThrow("Interaction::setYPtr, nInteraction differs from newY vector size ");

  // set y[index]
  if (isYAllocatedIn[index]) delete y[index];
  y[index] = newY;
  isYAllocatedIn[index] = false ;
}

void Interaction::setYOld(const std::vector<SimpleVector*>& newVector)
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

void Interaction::setYOldPtr(const std::vector<SimpleVector*>& newVector)
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
  if (newYOld->size() != nInteraction)
    RuntimeException::selfThrow("Interaction::setYOldPtr, nInteraction differs from newYOld vector size ");

  // set yOld[index]
  if (isYOldAllocatedIn[index]) delete yOld[index];
  yOld[index] = newYOld;
  isYOldAllocatedIn[index] = false ;
}

void Interaction::setLambda(const std::vector<SimpleVector*>& newVector)
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

void Interaction::setLambdaPtr(const std::vector<SimpleVector*>& newVector)
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
  if (newLambda->size() != nInteraction)
    RuntimeException::selfThrow("Interaction::setLambdaPtr, nInteraction differs from newLambda vector size ");

  // set lambda[index]
  if (isLambdaAllocatedIn[index]) delete lambda[index];
  lambda[index] = newLambda;
  isLambdaAllocatedIn[index] = false ;
}

void Interaction::setLambdaOld(const std::vector<SimpleVector*>& newVector)
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

void Interaction::setLambdaOldPtr(const std::vector<SimpleVector*>& newVector)
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
  if (newLambdaOld->size() != nInteraction)
    RuntimeException::selfThrow("Interaction::setLambdaOldPtr, nInteraction differs from newLambdaOld vector size ");

  // set lambdaOld[index]
  if (isLambdaOldAllocatedIn[index]) delete lambdaOld[index];
  lambdaOld[index] = newLambdaOld;
  isLambdaOldAllocatedIn[index] = false ;
}


void Interaction::setDynamicalSystems(const std::vector<DynamicalSystem*>& newVector)
{
  vectorDS.clear();
  vectorDS = newVector;
}

DynamicalSystem* Interaction::getDynamicalSystemPtr(const int& number)
{
  DynamicalSystem * tmpDS = NULL;
  vector<DynamicalSystem*>::iterator it;
  for (it = vectorDS.begin(); it != vectorDS.end(); ++it)
    if ((*it)->getNumber() == number) tmpDS = (*it);

  if (tmpDS == NULL)
    RuntimeException::selfThrow("Interaction::getDynamicalSystemPtr(number), there is no DS which number is " + number);

  return tmpDS;
}

DynamicalSystem Interaction::getDynamicalSystem(const int& number)
{
  return *getDynamicalSystemPtr(number);
}

void Interaction::setRelationPtr(Relation* newRelation)
{
  if (isRelationAllocatedIn) delete relation;
  relation = newRelation;
  isRelationAllocatedIn = false;
}

void Interaction::setNonSmoothLawPtr(NonSmoothLaw* newNslaw)
{
  if (isNsLawAllocatedIn) delete nslaw;
  nslaw = newNslaw;
  isNsLawAllocatedIn = false;
}

// --- OTHER FUNCTIONS ---

void Interaction::swapInMemory()
{
  IN("Interaction::swapInMemory(void)\n");

  for (unsigned int i = 0; i < y.size() ; i++)
  {
    *(yOld[i]) = *(y[i]) ;
    *(lambdaOld[i]) = *(lambda[i]);
    lambda[i]->zero();
  }

  OUT("Interaction::swapInMemory(void)\n");
}

void Interaction::display() const
{
  cout << "======= Interaction display =======" << endl;
  cout << "| id : " << id << endl;
  cout << "| number : " << number << endl;
  cout << "| Dynamical Systems linked to this Interaction : " << endl;
  for (unsigned int i = 0; i < vectorDS.size() ; i++) cout << vectorDS[i] << endl;
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
  nslaw = new ComplementarityConditionNSL();
  return nslaw;
}

NonSmoothLaw* Interaction::createRelayNSL(const double& c, const double& d)
{
  nslaw = new RelayNSL(c, d);
  return nslaw;
}

NonSmoothLaw* Interaction::createNewtonImpactLawNSL(const double& e)
{
  nslaw = new NewtonImpactLawNSL(e);
  return nslaw;
}

NonSmoothLaw* Interaction::createNewtonImpactFrictionNSL(const double& en, const double& et, const double& mu)
{
  nslaw = new NewtonImpactFrictionNSL(en, et, mu);
  return nslaw;
}

// --- XML RELATED FUNCTIONS ---

void Interaction::saveInteractionToXML()
{
  IN("Interaction::saveInteractionToXML\n");
  /*
   * save the data of the Interaction
   */

  if (interactionxml != NULL)
  {
    //  interactionxml->setDSConcerned( vectorDS );
    interactionxml->setId(id);
    interactionxml->setNumber(number);
    interactionxml->setNInter(nInteraction);
    interactionxml->setY(*(y[0]));
    interactionxml->setLambda(*(lambda[0]));
  }
  else RuntimeException::selfThrow("Interaction::saveInteractionToXML - object InteractionXML does not exist");

  /*
   * save the data of the Relation
   */
  if (relation->getType() == LINEARTIRELATION)
    (static_cast<LinearTIR*>(relation))->saveRelationToXML();
  else if (relation->getType() == LAGRANGIANLINEARRELATION)
    (static_cast<LagrangianLinearR*>(relation))->saveRelationToXML();
  else if (relation->getType() == LAGRANGIANNONLINEARRELATION)
    (static_cast<LagrangianNonLinearR*>(relation))->saveRelationToXML();
  else RuntimeException::selfThrow("Interaction::saveInteractionToXML - bad kind of Relation :" + relation->getType());
  /*
   * save the data of the NonSmoothLaw
   */

  if (nslaw->getType() == COMPLEMENTARITYCONDITIONNSLAW)
    (static_cast<ComplementarityConditionNSL*>(nslaw))->saveNonSmoothLawToXML();
  else if (nslaw->getType() == RELAYNSLAW)
    (static_cast<RelayNSL*>(nslaw))->saveNonSmoothLawToXML();
  else if (nslaw->getType() == NEWTONIMPACTLAWNSLAW)
    (static_cast<NewtonImpactLawNSL*>(nslaw))->saveNonSmoothLawToXML();
  else if (nslaw->getType() == NEWTONIMPACTFRICTIONNSLAW)
    (static_cast<NewtonImpactFrictionNSL*>(nslaw))->saveNonSmoothLawToXML();
  else RuntimeException::selfThrow("Interaction::saveInteractionToXML - bad kind of NonSmoothLaw : " + nslaw->getType());

  OUT("Interaction::saveInteractionToXML\n");
}

// Default (private) constructor
Interaction::Interaction():
  id("none"), number(0), nInteraction(0), nslaw(NULL), relation(NULL), interactionxml(NULL),
  isRelationAllocatedIn(false), isNsLawAllocatedIn(false)
{}
