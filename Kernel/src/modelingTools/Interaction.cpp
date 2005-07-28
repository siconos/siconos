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
  lambda(NULL), lambdaOld(NULL), nslaw(NULL), relation(NULL), interactionxml(NULL),
  isLambdaAllocatedIn(true), isLambdaOldAllocatedIn(true),
  isRelationAllocatedIn(true), isNsLawAllocatedIn(true)
{
  // Memory allocation and copy for simple vectors
  unsigned int size = newI.getY().size();
  y.resize(size, NULL);
  yOld.resize(size, NULL);

  for (unsigned int i = 0; i < size; i++)
  {
    y[i] = new SimpleVector(*(newI.getYPtr(i)));
    yOld[i] = new SimpleVector(*(newI.getYOldPtr(i)));
  }
  isYAllocatedIn.resize(size, true);
  isYOldAllocatedIn.resize(size, true);

  lambda = new SimpleVector(*(newI.getLambdaPtr()));
  lambdaOld = new SimpleVector(*(newI.getLambdaOldPtr()));

  status = newI.getStatus();

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
  id("none"), number(0), nInteraction(0), lambda(NULL), lambdaOld(NULL),
  nslaw(NULL), relation(NULL), interactionxml(interxml),
  isLambdaAllocatedIn(true), isLambdaOldAllocatedIn(true),
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
    for (unsigned int i = 0; i < size ; i++)
    {
      y[i] = new SimpleVector(nInteraction);
      yOld[i] = new SimpleVector(nInteraction);
    }
    isYAllocatedIn.resize(size, true);
    isYOldAllocatedIn.resize(size, true);

    lambda = new SimpleVector(nInteraction);
    lambdaOld = new SimpleVector(nInteraction);

    if (interactionxml->hasY()) *(y[0]) = interactionxml->getY();
    if (interactionxml->hasLambda()) *lambda = interactionxml->getLambda();

    // Old values are initialized with current values
    swapInMemory();

    // --- Status ---
    status.reserve(nInteraction);
    status = interactionxml->getStatus();

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
                         vector<int>* newStatus, vector<DynamicalSystem*> *dsConcerned):
  id(newId), number(newNumber), nInteraction(nInter), lambda(NULL), lambdaOld(NULL), nslaw(NULL),
  relation(NULL), interactionxml(NULL), isLambdaAllocatedIn(true), isLambdaOldAllocatedIn(true),
  isRelationAllocatedIn(false), isNsLawAllocatedIn(false)
{
  // Memory allocation for simple vectors

  // \todo : compute size using relative degree
  // for the moment, size = 2: we save y and yDot
  unsigned int size = 2;
  y.resize(size) ;
  yOld.resize(size);
  for (unsigned int i = 0; i < size ; i++)
  {
    y[i] = new SimpleVector(nInteraction);
    yOld[i] = new SimpleVector(nInteraction);
  }
  isYAllocatedIn.resize(size, true);
  isYOldAllocatedIn.resize(size, true);

  lambda = new SimpleVector(nInteraction);
  lambdaOld = new SimpleVector(nInteraction);

  status = *newStatus;

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
  }
  y.clear();
  yOld.clear();
  if (isLambdaAllocatedIn) delete lambda;
  lambda = NULL;
  if (isLambdaOldAllocatedIn) delete lambdaOld;
  lambdaOld = NULL;
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
    *(yOld[i]) = *(y[i]) ;

  *lambdaOld = *lambda;
  OUT("Interaction::swapInMemory(void)\n");
}

void Interaction::check(const double& time, const double& pasH)
{
  IN("Interaction::check(void)\n");
  unsigned int i;
  // relation->computeOutput(time);

  // Compute yp, predicted value for constrained variable, for contact detection
  // if contact (yp<0), status=1, else equal 0.

  if (nslaw->getType() == NEWTONIMPACTLAWNSLAW || nslaw->getType() == NEWTONIMPACTFRICTIONNSLAW)
  {
    SimpleVector *yDetection = new SimpleVector(*(yOld[0]));
    *yDetection += 0.5 * pasH * *(yOld[1]);
    for (i = 0; i < nInteraction; i++)
    {
      if ((*yDetection)(i) < 0.0) status[i] = 1;
    }
    delete yDetection;
  }
  else if (nslaw->getType() == COMPLEMENTARITYCONDITIONNSLAW)
  {
    for (i = 0; i < nInteraction; i++) status[i] = 1;
  }
  else
    RuntimeException::selfThrow("Interaction::check - not yet implemented for this NSLAW type :" + nslaw->getType());
  OUT("Interaction::checkInteraction(void)\n");
}


void Interaction::update(const double& time, const double& pasH)
{
  IN("Interaction::update(void)\n");
  unsigned int i;

  // Status update
  if (nslaw->getType() == NEWTONIMPACTLAWNSLAW || nslaw->getType() == NEWTONIMPACTFRICTIONNSLAW)
  {
    SimpleVector *yDetection = new SimpleVector(*(yOld[0]));
    *yDetection += 0.5 * pasH * *(yOld[1]);
    for (i = 0; i < nInteraction; i++)
    {
      if ((status[i] == 1) && ((*yDetection)(i) > 0.0)) status[i] = 0;
    }
    delete yDetection;
  }
  else if (nslaw->getType() == COMPLEMENTARITYCONDITIONNSLAW)
    // to do : check the relative degree
  {}
  else
    RuntimeException::selfThrow("Interaction::update - not yet implemented for this NSLAW type :" + nslaw->getType());
  OUT("Interaction::update(void)\n");
}

void Interaction::display() const
{
  cout << "======= Interaction display =======" << endl;
  cout << "| id : " << id << endl;
  cout << "| number : " << number << endl;
  cout << "| status : ";
  for (unsigned int i = 0; i < nInteraction; i++) cout << status[i] << " ";
  cout << endl;
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
  if (lambda != NULL) lambda->display();
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
    interactionxml->setStatus(status);
    interactionxml->setNInter(nInteraction);
    interactionxml->setY(*(y[0]));
    interactionxml->setLambda(*lambda);
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
  id("none"), number(0), nInteraction(0), lambda(NULL), lambdaOld(NULL), nslaw(NULL), relation(NULL), interactionxml(NULL),
  isLambdaAllocatedIn(false), isLambdaOldAllocatedIn(false),
  isRelationAllocatedIn(false), isNsLawAllocatedIn(false)
{}
