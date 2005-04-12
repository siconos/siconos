
#include "Interaction.h"

// include here the Relations for the static_casts
#include "LinearTIR.h"
#include "LagrangianLinearR.h"
#include "LagrangianNonLinearR.h"

// include here the NonSmoothLaws for the static_casts
#include "ComplementarityConditionNSL.h"
#include "RelayNSL.h"
#include "NewtonImpactLawNSL.h"
#include "NewtonImpactFrictionNSL.h"

#include "check.h"

Interaction::Interaction()
{
  this->id = "none";
  this->nInteraction = 0;

  this->vectorDS.clear();
  this->nslaw = NULL;
  this->relation = NULL;

  this->interactionxml = NULL;
}

Interaction::Interaction(InteractionXML* interxml)
{
  this->id = "none";
  this->nInteraction = 0;

  this->vectorDS.clear();
  this->nslaw = NULL;
  this->relation = NULL;

  this->interactionxml = interxml;
}

Interaction::~Interaction()
{
  if (relation != NULL)
  {
    delete relation;
  }
  if (nslaw != NULL)
  {
    delete nslaw;
  }
}


SimpleVector* Interaction::getYPtr(void)
{
  return &this->y;
}


SimpleVector* Interaction::getLambdaPtr(void)
{
  return &this->lambda;
}

SimpleVector* Interaction::getYOldPtr(void)
{
  return &this->yOld;
}

SimpleVector* Interaction::getLambdaOldPtr(void)
{
  return &this->lambdaOld;
}


DynamicalSystem* Interaction::getDynamicalSystem(int number)
{
  if (this->vectorDS[0]->getNumber() == number) return this->vectorDS[0];
  else if (this->vectorDS[1]->getNumber() == number) return this->vectorDS[1];
  else return NULL;

}

void Interaction::setDynamicalSystems(DynamicalSystem* ds1, DynamicalSystem* ds2)
{
  this->vectorDS.clear();
  this->vectorDS.push_back(ds1);
  this->vectorDS.push_back(ds2);
}


////////////////////////

void Interaction::initialize()
{
  IN("Interaction::initialize()\n");

  this->y = SimpleVector(this->nInteraction);
  this->lambda = SimpleVector(this->nInteraction);
  this->yOld = SimpleVector(this->nInteraction);
  this->yDot = SimpleVector(this->nInteraction);
  this->yDotOld = SimpleVector(this->nInteraction);
  this->lambdaOld = SimpleVector(this->nInteraction);

  this->display();

  OUT("Interaction::initialize()\n");
}


void Interaction::swapInMemory(void)
{
  IN("Interaction::swapInMemory(void)\n");
  yOld = y;
  lambdaOld = lambda;
  yDotOld = yDot;
  OUT("Interaction::swapInMemory(void)\n");
}




void Interaction::check(double time)
{
  IN("Interaction::check(void)\n");
  int i;
  if (this->nslaw->getType() == COMPLEMENTARITYCONDITIONNSLAW)
  {
    this->relation->computeOutput(time);
    for (i = 0; i < this->nInteraction; i++)
    {
      if (((this->yOld)(i) < 0.0) || (this->status[i] == 1))
      {
        this->status[i] = 1;
      }
    }
  }
  else if (this->nslaw->getType() == NEWTONIMPACTLAWNSLAW)
  {
    this->relation->computeOutput(time);
    for (i = 0; i < this->nInteraction; i++)
    {
      if (((this->yOld)(i) < 0.0) || (this->status[i] == 1))
      {
        this->status[i] = 1;
      }
    }
  }
  else if (this->nslaw->getType() == NEWTONIMPACTFRICTIONNSLAW)
  {
    this->relation->computeOutput(time);
    for (i = 0; i < this->nInteraction; i++)
    {
      if (((this->yOld)(i) < 0.0) || (this->status[i] == 1))
      {
        this->status[i] = 1;
      }
    }
  }
  else
    RuntimeException::selfThrow("Interaction::check - not yet implemented for this NSLAW type :" + nslaw->getType());


  OUT("Interaction::checkInteraction(void)\n");
}


void Interaction::update(double time)
{
  IN("Interaction::update(void)\n");
  int i;
  // Status update
  if (this->nslaw->getType() == COMPLEMENTARITYCONDITIONNSLAW)
  {
    for (i = 0; i < this->nInteraction; i++)
    {
      if ((status[i] == 1) && ((this->yDot)(i) > 0.0))
      {
        this->status[i] = 0;
      }
    }
  }
  else if (this->nslaw->getType() == NEWTONIMPACTLAWNSLAW)
  {
    for (i = 0; i < this->nInteraction; i++)
    {
      if ((status[i] == 1) && ((this->yDot)(i) > 0.0))
      {
        this->status[i] = 0;
      }
    }
  }
  else if (this->nslaw->getType() == NEWTONIMPACTFRICTIONNSLAW)
  {
    for (i = 0; i < this->nInteraction; i++)
    {
      if ((status[i] == 1) && ((this->yDot)(i) > 0.0))
      {
        this->status[i] = 0;
      }
    }
  }

  /*
   * only done when the interaction is active
   */
  // Input update
  //this->relation->computeInput( time );
  OUT("Interaction::update(void)\n");
}


void Interaction::fillInteractionWithInteractionXML()
{
  IN("Interaction::fillInteractionWithInteractionXML\n");
  if (this->interactionxml != NULL)
  {
    if (this->interactionxml->hasId()) this->id = this->interactionxml->getId();
    this->number = this->interactionxml->getNumber();

    this->nInteraction = this->interactionxml->getNInter();

    this->status.resize(this->nInteraction);
    this->status = this->interactionxml->getStatus();

    if (this->interactionxml->hasY()) this->y = this->interactionxml->getY();
    else this->y = /*SiconosVector*/SimpleVector(this->nInteraction);

    if (this->interactionxml->hasLambda()) this->lambda = this->interactionxml->getLambda();
    else this->lambda = /*SiconosVector*/SimpleVector(this->nInteraction);

    //this->display();
  }
  else RuntimeException::selfThrow("Interaction::fillInteractionWithInteractionXML - InteractionXML object not exists");
  OUT("Interaction::fillInteractionWithInteractionXML\n");
}

void Interaction::display() const
{
  cout << "-----------------------------------------------------" << endl;
  cout << "____ data of the Interaction " << endl;
  cout << "| id : " << this->id << endl;
  cout << "| number : " << this->number << endl;
  cout << "| status : ";
  for (int i = 0; i < nInteraction; i++) cout << status[i] << " ";
  cout << endl;
  cout << "| DS linked to this Interaction : ";
  for (int i = 0; i < nInteraction; i++) cout << vectorDS[i]->getNumber() << " ";
  cout << endl;
  cout << "| y : " << endl;
  this->y.display();
  cout << "| yDot : " << endl;
  this->yDot.display();
  cout << "| yOld : " << endl;
  this->yOld.display();
  cout << "| yDotOld : " << endl;
  this->yDotOld.display();
  cout << "| lambda : " << endl;
  this->lambda.display();
  cout << "-----------------------------------------------------" << endl << endl;
}

void Interaction::linkInteractionWithInteractionXML()
{
  IN("Interaction::linkInteractionXML\n");
  // get the RelationXML object then create the Relation for this RelationXML and add this Relation to the Interaction

  // with the data of the XML object, we know the type of Relation, so we can instanciate the right type of Relation
  // Relation - LinearTIR
  if (this->interactionxml->getRelationXML()->getType() == LINEAR_TIME_INVARIANT_RELATION_TAG)
  {
    cout << "Interaction - Relation LTI ..." << endl;
    // creation of the LinearTIR with this constructor and call of a method to fill
    this->relation = new LinearTIR();
    static_cast<LinearTIR*>(this->relation)->createRelation(static_cast<LinearTIRXML*>(this->interactionxml->getRelationXML()));
    this->relation->setInteraction(this);
  }
  // Relation - LagrangianLinearR
  else if (this->interactionxml->getRelationXML()->getType() == LAGRANGIAN_LINEAR_RELATION_TAG)
  {
    cout << "Interaction - Relation LagrangianLinear ..." << endl;
    // creation of the LagrangianLinearR with this constructor and call of a method to fill
    this->relation = new LagrangianLinearR();
    static_cast<LagrangianLinearR*>(this->relation)->createRelation(static_cast<LagrangianLinearRXML*>(this->interactionxml->getRelationXML()));
    this->relation->setInteraction(this);
  }

  // Relation - LagrangianNonLinearR
  else if (this->interactionxml->getRelationXML()->getType() == LAGRANGIAN_NON_LINEAR_RELATION_TAG)
  {
    cout << "Interaction - Relation LagrangianNonLinear ..." << endl;
    // creation of the LagrangianNonLinearR with this constructor and call of a method to fill
    this->relation = new LagrangianNonLinearR();
    static_cast<LagrangianNonLinearR*>(this->relation)->createRelation(static_cast<LagrangianNonLinearRXML*>(this->interactionxml->getRelationXML()));
    this->relation->setInteraction(this);
  }

  // for the other Relation, we must have the xxxRelationXML.h .cpp objects needed, and the objects in the Relations inherited of the platform
  /*
   *  other "if" to create LagrangianNonLinearR and LagrangianLinearR
   *
   */

  /* ================================================ */

  // with the data of the XML object, we know the type of NonSmoothLaw, so we can instanciate the right type of NonSmoothLaw

  // NonSmoothLaw - ComplementarityConditionNSL
  if (this->interactionxml->getNonSmoothLawXML()->getType() == COMPLEMENTARITY_CONDITION_NSLAW_TAG)
  {
    // creation of the ComplementarityConditionNSL with this constructor and call of a method to fill
    this->nslaw = new ComplementarityConditionNSL();
    static_cast<ComplementarityConditionNSL*>(this->nslaw)->createNonSmoothLaw(static_cast<ComplementarityConditionNSLXML*>(this->interactionxml->getNonSmoothLawXML()));
  }
  // NonSmoothLaw - RelayNSL
  else if (this->interactionxml->getNonSmoothLawXML()->getType() == RELAY_NSLAW_TAG)
  {
    // creation of the RelayNSL with this constructor and call of a method to fill
    this->nslaw = new RelayNSL();
    static_cast<RelayNSL*>(this->nslaw)->createNonSmoothLaw(static_cast<RelayNSLXML*>(this->interactionxml->getNonSmoothLawXML()));
  }
  // NonSmoothLaw - NewtonImpactLawNSL
  else if (this->interactionxml->getNonSmoothLawXML()->getType() == NEWTON_IMPACT_LAW_NSLAW_TAG)
  {
    // creation of the RelayNSL with this constructor and call of a method to fill
    this->nslaw = new NewtonImpactLawNSL();
    static_cast<NewtonImpactLawNSL*>(this->nslaw)->createNonSmoothLaw(static_cast<NewtonImpactLawNSLXML*>(this->interactionxml->getNonSmoothLawXML()));
  }
  else if (this->interactionxml->getNonSmoothLawXML()->getType() == NEWTON_IMPACT_FRICTION_NSLAW_TAG)
  {
    // creation of the RelayNSL with this constructor and call of a method to fill
    this->nslaw = new NewtonImpactFrictionNSL();
    static_cast<NewtonImpactFrictionNSL*>(this->nslaw)->createNonSmoothLaw(static_cast<NewtonImpactFrictionNSLXML*>(this->interactionxml->getNonSmoothLawXML()));
  }
  OUT("Interaction::linkInteractionXML\n");
}

void Interaction::saveInteractionToXML()
{
  IN("Interaction::saveInteractionToXML\n");
  cout << "##### Interaction::saveInteractionToXML" << endl;;
  /*
   * save the data of the Interaction
   */
  if (this->interactionxml != NULL)
  {
    this->interactionxml->setDSConcerned(this->vectorDS);

    this->interactionxml->setId(this->id);
    this->interactionxml->setNumber(this->number);
    //this->interactionxml->setStatus( this->status[0] );
    this->interactionxml->setStatus(this->status);
    this->interactionxml->setNInter(this->nInteraction);

    this->interactionxml->setY(&(this->y));
    this->interactionxml->setLambda(&(this->lambda));
  }
  else RuntimeException::selfThrow("Interaction::saveInteractionToXML - object InteractionXML does not exist");

  /*
   * save the data of the Relation
   */
  if (this->relation->getType() == LINEARTIRELATION)
    (static_cast<LinearTIR*>(this->relation))->saveRelationToXML();
  else if (this->relation->getType() == LAGRANGIANLINEARRELATION)
    (static_cast<LagrangianLinearR*>(this->relation))->saveRelationToXML();
  else if (this->relation->getType() == LAGRANGIANNONLINEARRELATION)
    (static_cast<LagrangianNonLinearR*>(this->relation))->saveRelationToXML();
  else RuntimeException::selfThrow("Interaction::saveInteractionToXML - bad kind of Relation :" + this->relation->getType());
  /*
   * save the data of the NonSmoothLaw
   */

  cout << "    save of the NSLaw : " << this->nslaw << endl;
  cout << nslaw->getType() << endl;
  if (this->nslaw->getType() == COMPLEMENTARITYCONDITIONNSLAW)
    (static_cast<ComplementarityConditionNSL*>(this->nslaw))->saveNonSmoothLawToXML();
  else if (this->nslaw->getType() == RELAYNSLAW)
    (static_cast<RelayNSL*>(this->nslaw))->saveNonSmoothLawToXML();
  else if (this->nslaw->getType() == NEWTONIMPACTLAWNSLAW)
    (static_cast<NewtonImpactLawNSL*>(this->nslaw))->saveNonSmoothLawToXML();
  else if (this->nslaw->getType() == NEWTONIMPACTFRICTIONNSLAW)
    (static_cast<NewtonImpactFrictionNSL*>(this->nslaw))->saveNonSmoothLawToXML();
  else RuntimeException::selfThrow("Interaction::saveInteractionToXML - bad kind of NonSmoothLaw : " + this->nslaw->getType());

  OUT("Interaction::saveInteractionToXML\n");
}

void Interaction::createInteraction(InteractionXML * interactionXML, int number, int nInter,
                                    vector<int>* status, vector<DynamicalSystem*>* dsConcerned)//, NonSmoothDynamicalSystem * nsds)
{
  IN("Interaction::createInteraction\n");
  if (interactionXML != NULL)
  {
    this->id = "none";
    this->nInteraction = 0;

    this->vectorDS.clear();
    this->nslaw = NULL;
    this->relation = NULL;

    this->interactionxml = interactionXML;

    this->fillInteractionWithInteractionXML();
    this->linkInteractionWithInteractionXML();
  }
  else
  {
    this->number = number;
    this->nInteraction = nInter;
    this->status = *status;
    if (dsConcerned != NULL) this->vectorDS = *dsConcerned;
    else RuntimeException::selfThrow("Interaction::createInteraction - The dsConcerned are not given");
  }
  OUT("Interaction::createInteraction\n");
}


Relation* Interaction::createLagrangianLinearR(SiconosMatrix* H, SiconosVector* b)
{
  this->relation = new LagrangianLinearR();
  static_cast<LagrangianLinearR*>(this->relation)->createRelation(NULL, H, b);
  this->relation->setInteraction(this);
  return this->relation;
}

Relation* Interaction::createLagrangianNonLinearR(string computeInput, string computeOutput)
{
  this->relation = new LagrangianNonLinearR();
  static_cast<LagrangianNonLinearR*>(this->relation)->createRelation(NULL, computeInput, computeOutput);
  this->relation->setInteraction(this);
  return this->relation;
}

Relation* Interaction::createLinearTIR(SiconosMatrix* C, SiconosMatrix* D,
                                       SiconosMatrix* E, SiconosVector* a)
{
  this->relation = new LinearTIR();
  static_cast<LinearTIR*>(this->relation)->createRelation(NULL, C, D, E, a);
  this->relation->setInteraction(this);
  return this->relation;
}


NonSmoothLaw* Interaction::createComplementarityConditionNSL()
{
  this->nslaw = new ComplementarityConditionNSL();
  static_cast<ComplementarityConditionNSL*>(this->nslaw)->createNonSmoothLaw(NULL);
  return this->nslaw;
}

NonSmoothLaw* Interaction::createRelayNSL(double c, double d)
{
  this->nslaw = new RelayNSL();
  static_cast<RelayNSL*>(this->nslaw)->createNonSmoothLaw(NULL, c, d);
  return this->nslaw;
}

NonSmoothLaw* Interaction::createNewtonImpactLawNSL(double e)
{
  this->nslaw = new NewtonImpactLawNSL();
  static_cast<NewtonImpactLawNSL*>(this->nslaw)->createNonSmoothLaw(NULL, e);
  return this->nslaw;
}

NonSmoothLaw* Interaction::createNewtonImpactFrictionNSL(double en, double et, double mu)
{
  this->nslaw = new NewtonImpactFrictionNSL();
  static_cast<NewtonImpactFrictionNSL*>(this->nslaw)->createNonSmoothLaw(NULL, en, et, mu);
  return this->nslaw;
}

