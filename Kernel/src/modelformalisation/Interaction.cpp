
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

//$Log: Interaction.cpp,v $
//Revision 1.57  2005/03/22 15:55:04  jbarbier
//- class NewtonImpactFriction non smooth law added to the kernel
//
//- xml schema modified for this new class
//- xml schema modified to accept a "joker" for further use of a LMGC90 mechanical plugin
//
//- new test added for the loading/saving of a NewtonImpactFrictionNSL
//
//Revision 1.56  2005/03/08 12:41:36  jbarbier
//- constant variables files modified :
//Some constants added in SiconosConst
//
//all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//
//Revision 1.55  2005/03/07 13:17:19  jbarbier
//- new test : Ball2D, with a ball moving in a 2D system
//
//- another constant variables moved/refactored in XMLTagsName
//- making uniform the name of the constant variables
//
//Revision 1.54  2005/03/01 10:38:19  jbarbier
//- RollingBalls sample is OK
//
//Revision 1.53  2005/02/15 15:15:32  charlety
//
//_ modified some very slow functions to increase performance
//
//Revision 1.52  2005/02/04 14:52:44  jbarbier
//- Rolling balls in progress (contact is detected)
//
//- time data is given in parameter for computInput and Output in the Relation. Parameter is sent by methods of the OneStepNSProblem
//
//Revision 1.51  2005/02/01 11:08:41  charlety
//
//_ some displays of values during computations suppressed.
//
//Revision 1.50  2005/01/20 14:44:48  jbarbier
//- NSDS class renamed NonSmoothDynamicalSystem
//
//- code reduce, some comments remove
//
//Revision 1.49  2005/01/10 17:06:37  jbarbier
//- attribute "size" is now unused in the code
//
//- xml schema v1.2 is in progress
//
//Revision 1.48  2004/12/08 12:49:37  jbarbier
//- changes in the XML Schema, respect of the recommandations of the W3C
//version 1.1
//
//- changes in all balises DS, Relation, NSLaw, OneStepIntegrator, OneStepNSProblem
//in the XML files into specific names like LagrangianNLDS, LinearSystemDS, ...
//for the DS
//
//Revision 1.47  2004/09/28 08:21:27  jbarbier
//
//- manual creation of the BouncingBall example successful
//
//Revision 1.46  2004/09/23 14:09:23  jbarbier
//- modification of the integrators, the attribute r is always optional.
//
//- modification of the LagrangianNonLinearR. computeInput and computeOutput are
//required.
//
//Revision 1.45  2004/09/22 11:16:28  charlety
//
//_ revision of Doxygen comments in modelformalisation
//
//Revision 1.44  2004/09/21 11:49:09  jbarbier
//- correction in the XML save for a manual construction of the platform :
//    DS_Concerned of the Interaction
//    DS_Concerned of the Integrator
//
//- test updated for these changes
//
//Revision 1.43  2004/09/10 11:26:10  charlety
//
//_ Integration of the new version of the SiconosVector in the platform. the class simpleVector is used mostly to replace old SiconosVector. When a vector can be composite or simple, like the state of a dynamical system, a pointer on SiconosVector is used, and the vector is initialized simple or composite when the system is initialized.
//
//_ All the tests which worked with the previous version of the vector are OK with the new version.
//
//_ Example SICONOS and bouncingBall are OK
//
//_ some comments have still to be adapted to NewSiconosVector .
//
//_ This version of NewSiconosVector could be called 0.9. some details have to be fixed, it will be done before the end of September.
//
//Revision 1.42  2004/09/03 14:41:41  jbarbier
//- new functions to create the boundary condition of the dynamical systems
//- new functions to add an interaction to a NonSmoothDynamicalSystem
//- new functions to create the relation and the non-smooth law of an interaciton
//
//Revision 1.41  2004/08/18 14:37:18  jbarbier
//- creation of Model, NonSmoothDynamicalSystem, Strategy(TimeStepping and EventDriven) and
//DynamicalSystem available when the creation is in a command program
//
//Revision 1.40  2004/08/17 15:12:37  jbarbier
//- methods createDynamicalSystem, createBoundaryCondition, createInteraction,
//createRelation and createNSLaw completed with the required attributes
//
//Revision 1.39  2004/08/12 11:55:14  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
//Revision 1.38  2004/08/03 12:07:11  jbarbier
//- all test on th eModel are successfull
//
//- new tests on the Model with the opening of XML file
//
//- link TimeDiscretisation -> Strategy
//
//- attribute T of the Model is now optional
//
//Revision 1.37  2004/07/29 14:25:36  jbarbier
