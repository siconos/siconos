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
  y(NULL), yDot(NULL), lambda(NULL), yOld(NULL), yDotOld(NULL), lambdaOld(NULL),
  nslaw(NULL), relation(NULL), interactionxml(NULL), isYAllocatedIn(true), isYDotAllocatedIn(true),
  isLambdaAllocatedIn(true), isYOldAllocatedIn(true), isYDotOldAllocatedIn(true), isLambdaOldAllocatedIn(true),
  isRelationAllocatedIn(true), isNsLawAllocatedIn(true)
{
  // Memory allocation and copy for simple vectors
  y = new SimpleVector(nInteraction);
  yDot = new SimpleVector(nInteraction);
  lambda = new SimpleVector(nInteraction);
  yOld = new SimpleVector(nInteraction);
  yDotOld = new SimpleVector(nInteraction);
  lambdaOld = new SimpleVector(nInteraction);
  *y = *newI.getYPtr();
  *yDot = *newI.getYDotPtr();
  *lambda = *newI.getLambdaPtr();
  *yOld = *newI.getYOldPtr();
  *yDotOld = *newI.getYDotOldPtr();
  *lambdaOld = *newI.getLambdaOldPtr();

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
  if (relationType == LINEARTIRELATION)
    relation = new LinearTIR();
  else if (relationType == LAGRANGIANLINEARRELATION)
    relation = new LagrangianLinearR();
  else if (relationType == LAGRANGIANNONLINEARRELATION)
    relation = new LagrangianNonLinearR();
  else RuntimeException::selfThrow("Interaction::copy constructor, unknown relation type " + relation->getType());
  *relation = *newI.getRelationPtr();

  // Remark (FP): no link with newI xml object -> useless and avoid new/delete bugs
}

// --- XML constructor ---
Interaction::Interaction(InteractionXML* interxml):
  id("none"), number(0), nInteraction(0), y(NULL), yDot(NULL), lambda(NULL), yOld(NULL),
  yDotOld(NULL), lambdaOld(NULL), nslaw(NULL), relation(NULL), interactionxml(interxml),
  isYAllocatedIn(true), isYDotAllocatedIn(true), isLambdaAllocatedIn(true), isYOldAllocatedIn(true),
  isYDotOldAllocatedIn(true), isLambdaOldAllocatedIn(true), isRelationAllocatedIn(true), isNsLawAllocatedIn(true)
{
  if (interactionxml != NULL)
  {
    vectorDS.clear();
    if (interactionxml->hasId()) id = interactionxml->getId();
    number = interactionxml->getNumber();
    nInteraction = interactionxml->getNInter();
    // Memory allocation for simple vectors
    y = new SimpleVector(nInteraction);
    yDot = new SimpleVector(nInteraction);
    lambda = new SimpleVector(nInteraction);
    yOld = new SimpleVector(nInteraction);
    yDotOld = new SimpleVector(nInteraction);
    lambdaOld = new SimpleVector(nInteraction);
    if (interactionxml->hasY()) *y = interactionxml->getY();
    if (interactionxml->hasLambda()) *lambda = interactionxml->getLambda();
    // Remark: no yDot in xml

    // Old values are initialized with current values
    swapInMemory();

    // --- Status ---
    status.resize(nInteraction);
    status = interactionxml->getStatus();

    // --- Dynamical Systems ---
    // do nothing. The way DS are loaded from xml is to be reviewed

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
    {
      relation = new LinearTIR(interactionxml->getRelationXML());
      relation->setInteraction(this);
    }
    // Lagrangian linear relation
    else if (relationType == LAGRANGIAN_LINEAR_RELATION_TAG)
    {
      relation = new LagrangianLinearR(interactionxml->getRelationXML());
      relation->setInteraction(this);
    }
    // Lagrangian non-linear relation
    else if (relationType == LAGRANGIAN_NON_LINEAR_RELATION_TAG)
    {
      relation = new LagrangianNonLinearR(interactionxml->getRelationXML());
      relation->setInteraction(this);
    }
    else RuntimeException::selfThrow("Interaction::xml constructor, unknown relation type " + relation->getType());
  }
  else RuntimeException::selfThrow("Interaction::xml constructor, xmlfile = NULL");
}

// --- Constructor from a set of data ---

Interaction::Interaction(const string& newId, const int& newNumber, const int& nInter,
                         vector<int>* newStatus,
                         vector<DynamicalSystem*> *dsConcerned):
  id(newId), number(newNumber), nInteraction(nInter), y(NULL), yDot(NULL),
  lambda(NULL), yOld(NULL), yDotOld(NULL), lambdaOld(NULL), nslaw(NULL),
  relation(NULL), interactionxml(NULL), isYAllocatedIn(true), isYDotAllocatedIn(true),
  isLambdaAllocatedIn(true), isYOldAllocatedIn(true), isYDotOldAllocatedIn(true),
  isLambdaOldAllocatedIn(true), isRelationAllocatedIn(false), isNsLawAllocatedIn(false)
{
  // Memory allocation and copy for simple vectors
  y = new SimpleVector(nInteraction);
  yDot = new SimpleVector(nInteraction);
  lambda = new SimpleVector(nInteraction);
  yOld = new SimpleVector(nInteraction);
  yDotOld = new SimpleVector(nInteraction);
  lambdaOld = new SimpleVector(nInteraction);

  status.clear();
  status = *newStatus;

  vectorDS.clear();
  if (dsConcerned != NULL) vectorDS = *dsConcerned;
  else RuntimeException::selfThrow("Interaction::createInteraction - The dsConcerned are not given");

  // Remark(FP): neither nslaw nor relation are created in this constructor -> todo?
}

// --- DESTRUCTOR ---
Interaction::~Interaction()
{
  if (isYAllocatedIn)
  {
    delete y;
    y = NULL;
  }
  if (isYDotAllocatedIn)
  {
    delete yDot;
    yDot = NULL;
  }
  if (isLambdaAllocatedIn)
  {
    delete lambda;
    lambda = NULL;
  }
  if (isYOldAllocatedIn)
  {
    delete yOld;
    yOld = NULL;
  }
  if (isYDotOldAllocatedIn)
  {
    delete yDotOld;
    yDotOld = NULL;
  }
  if (isLambdaOldAllocatedIn)
  {
    delete lambdaOld;
    lambdaOld = NULL;
  }
  if (isRelationAllocatedIn)
  {
    delete relation;
    relation = NULL;
  }
  if (isNsLawAllocatedIn)
  {
    delete nslaw;
    nslaw = NULL;
  }
}

// --- GETTERS/SETTERS ---
void Interaction::setDynamicalSystems(DynamicalSystem* ds1, DynamicalSystem* ds2)
{
  vectorDS.clear();
  vectorDS.push_back(ds1);
  vectorDS.push_back(ds2);
}

DynamicalSystem* Interaction::getDynamicalSystemPtr(const int& number)
{
  if (vectorDS[0]->getNumber() == number) return vectorDS[0];
  else if (vectorDS[1]->getNumber() == number) return vectorDS[1];
  else return NULL;
}

DynamicalSystem Interaction::getDynamicalSystem(const int& number)
{
  if (number == 0 || number == 1) return *vectorDS[number];
  else RuntimeException::selfThrow("Interaction: getDynamicalSystem, unknown DS number");
}

// --- OTHER FUNCTIONS ---

void Interaction::swapInMemory()
{
  IN("Interaction::swapInMemory(void)\n");
  *yOld = *y;
  *lambdaOld = *lambda;
  *yDotOld = *yDot;
  OUT("Interaction::swapInMemory(void)\n");
}

void Interaction::check(const double& time, const double& pasH)
{
  IN("Interaction::check(void)\n");
  int i;
  relation->computeOutput(time);

  // Compute yp, predicted value for constrained variable, for contact detection

  if (nslaw->getType() == COMPLEMENTARITYCONDITIONNSLAW || nslaw->getType() == NEWTONIMPACTLAWNSLAW || nslaw->getType() == NEWTONIMPACTFRICTIONNSLAW)
  {
    SimpleVector *yDetection = new SimpleVector(yOld->size());
    *yDetection = *yOld + pasH * *yDotOld;
    for (i = 0; i < nInteraction; i++)
    {
      if ((*yDetection)(i) < 0.0) status[i] = 1;
    }
    delete yDetection;
  }
  else
    RuntimeException::selfThrow("Interaction::check - not yet implemented for this NSLAW type :" + nslaw->getType());
  OUT("Interaction::checkInteraction(void)\n");
}


void Interaction::update(const double& time, const double& pasH)
{
  IN("Interaction::update(void)\n");
  int i;
  // Status update
  if (nslaw->getType() == COMPLEMENTARITYCONDITIONNSLAW || nslaw->getType() == NEWTONIMPACTLAWNSLAW || nslaw->getType() == NEWTONIMPACTFRICTIONNSLAW)
  {
    SimpleVector *yDetection = new SimpleVector(y->size());
    *yDetection = *yOld + pasH * *yDotOld;
    for (i = 0; i < nInteraction; i++)
    {
      if ((status[i] == 1) && ((*yDetection)(i) > 0.0)) status[i] = 0;
    }
    delete yDetection;
  }
  else
    RuntimeException::selfThrow("Interaction::update - not yet implemented for this NSLAW type :" + nslaw->getType());
  OUT("Interaction::update(void)\n");
}

void Interaction::display() const
{
  cout << "-----------------------------------------------------" << endl;
  cout << "____ data of the Interaction " << endl;
  cout << "| id : " << id << endl;
  cout << "| number : " << number << endl;
  cout << "| status : ";
  for (int i = 0; i < nInteraction; i++) cout << status[i] << " ";
  cout << endl;
  cout << "| DS linked to this Interaction : ";
  for (int i = 0; i < nInteraction; i++) cout << vectorDS[i]->getNumber() << " ";
  cout << endl;
  cout << "| y : " << endl;
  if (y != NULL) y->display();
  else cout << "->NULL" << endl;
  cout << "| yDot : " << endl;
  if (yDot != NULL) yDot->display();
  else cout << "->NULL" << endl;
  cout << "| yOld : " << endl;
  if (yOld != NULL) yOld->display();
  else cout << "->NULL" << endl;
  cout << "| yDotOld : " << endl;
  if (yDotOld != NULL) yDotOld->display();
  else cout << "->NULL" << endl;
  cout << "| lambda : " << endl;
  if (lambda != NULL) lambda->display();
  else cout << "->NULL" << endl;
  cout << "-----------------------------------------------------" << endl << endl;
}

Relation* Interaction::createLagrangianLinearR(SiconosMatrix* H, SimpleVector* b)
{
  relation = new LagrangianLinearR(H, b);
  relation->setInteraction(this);
  return relation;
}

Relation* Interaction::createLagrangianNonLinearR(const string& computeInput, const string& computeOutput)
{
  relation = new LagrangianNonLinearR(computeInput, computeOutput);
  relation->setInteraction(this);
  return relation;
}

Relation* Interaction::createLinearTIR(SiconosMatrix* C, SiconosMatrix* D,
                                       SiconosMatrix* E, SiconosVector* a)
{
  relation = new LinearTIR(C, D, E, a);
  relation->setInteraction(this);
  return relation;
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
    interactionxml->setDSConcerned(vectorDS);

    interactionxml->setId(id);
    interactionxml->setNumber(number);
    interactionxml->setStatus(status);
    interactionxml->setNInter(nInteraction);

    interactionxml->setY(*y);
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
  id("none"), number(0), nInteraction(0), y(NULL), yDot(NULL), lambda(NULL), yOld(NULL),
  yDotOld(NULL), lambdaOld(NULL), nslaw(NULL), relation(NULL), interactionxml(NULL),  isYAllocatedIn(false), isYDotAllocatedIn(false),
  isLambdaAllocatedIn(false), isYOldAllocatedIn(false), isYDotOldAllocatedIn(false), isLambdaOldAllocatedIn(false),
  isRelationAllocatedIn(false), isNsLawAllocatedIn(false)
{}
