//$Id: Strategy.cpp,v 1.52 2005/03/11 15:06:20 jbarbier Exp $

#include "Strategy.h"
#include "Moreau.h"
#include "Lsodar.h"
#include "Adams.h"
#include "LagrangianTIDS.h"

#include "LCP.h"
#include "QP.h"
#include "Relay.h"

#include "check.h"

Strategy::Strategy()
{
  IN("Strategy::Strategy()\n");
  this->timeDiscretisation = NULL;
  this->integratorVector.clear();
  this->nsProblem = NULL;

  this->strategyxml = NULL;
  this->model = NULL;
  OUT("Strategy::Strategy()\n");
}

Strategy::Strategy(Strategy* str)
{
  IN("Strategy::Strategy( Strategy* str )\n");
  this->strategyType = str->getType();
  this->timeDiscretisation = str->getTimeDiscretisation();
  this->integratorVector = str->getOneStepIntegrators();
  this->nsProblem = str->getOneStepNSProblem();
  this->strategyxml = str->getStrategyXML();
  this->model = str->getModel();
  OUT("Strategy::Strategy( Strategy* str )\n");
}

Strategy::Strategy(StrategyXML* strxml, Model *model)
{
  this->timeDiscretisation = NULL;
  this->integratorVector.clear();
  this->nsProblem = NULL;

  this->strategyxml = strxml;
  this->model = model;
}

Strategy::~Strategy()
{
  if (this->nsProblem != NULL) delete this->nsProblem;
  if (this->integratorVector.size() > 0)
  {
    for (int i = 0; i < this->integratorVector.size(); i++)
    {
      delete this->integratorVector[i];
    }
    this->integratorVector.clear();
  }
  if (this->timeDiscretisation != NULL)
    delete this->timeDiscretisation;
}

OneStepIntegrator* Strategy::getOneStepIntegrator(int nb) const
{
  if (nb < this->integratorVector.size())
  {
    return this->integratorVector[nb];
  }
  RuntimeException::selfThrow("Strategy - getIntegrator : \'nb\' is out of range");
}



void Strategy::computeFreeState(void)
{
  IN("Strategy::computeFreeState\n");
  //cout<<"this->integratorVector.size()" <<this->integratorVector.size()<<endl;
  for (int i = 0; i < this->integratorVector.size(); i++)
  {
    this->integratorVector[i]->computeFreeState();
  }

  OUT("Strategy::computeFreeState\n");
}

void Strategy::nextStep(void)
{
  this->timeDiscretisation->increment();

  for (int i = 0; i < this->integratorVector.size(); i++)
  {
    this->integratorVector[i]->nextStep();
  }
  if (this->nsProblem != NULL)
    this->nsProblem->nextStep();

}


void Strategy::formaliseOneStepNSProblem()
{
  // formalise the OneStepNSProblem
  if (this->nsProblem != NULL)this->nsProblem->formalize(this->model->getCurrentT());
}

void Strategy::computeOneStepNSProblem(void)
{
  // compute the OneStepNSProblem
  if (this->nsProblem != NULL)this->nsProblem->compute();
}

void Strategy::updateState()
{
  // compute the OneStepNSProblem
  if (this->nsProblem != NULL)this->nsProblem->updateState();


  for (int i = 0; i < this->integratorVector.size(); i++)
  {
    //    cout<<"### Strategy::updateState (Ufree) :"<<endl;
    //    static_cast<LagrangianTIDS*>(this->integratorVector[i]->getDynamicalSystem())->getQ().display();
    //    static_cast<LagrangianTIDS*>(this->integratorVector[i]->getDynamicalSystem())->getVelocity().display();

    this->integratorVector[i]->updateState();

    //    cout<<"### Strategy::updateState (U) :"<<endl;
    //    static_cast<LagrangianTIDS*>(this->integratorVector[i]->getDynamicalSystem())->getQ().display();
    //    static_cast<LagrangianTIDS*>(this->integratorVector[i]->getDynamicalSystem())->getVelocity().display();
  }
  //  cout<<"                    <<Press Enter>>"<<endl;
  //  getchar();
}

void Strategy::initialize()
{
  // initialization of the TimeDiscretisation
  this->timeDiscretisation->init(this->model->getT0(), this->model->getFinalT());
  this->timeDiscretisation->display();

  // initialization of the OneStepIntegrators
  for (int i = 0; i < this->integratorVector.size(); i++)
  {
    this->integratorVector[i]->initialize();
  }


  // initialization of the OneStepNSProblem
  if (this->nsProblem != NULL)
    this->nsProblem->initialize();

  cout << "timediscretisation after init" << endl;
  this->timeDiscretisation->display();

}



OneStepIntegrator* Strategy::getIntegratorOfDS(int numberDS)
{
  for (int i = 0; i < this->integratorVector.size(); i++)
  {
    if (this->integratorVector[i]->getDynamicalSystem()->getNumber() == numberDS)
    {
      return this->integratorVector[i];
      break;
    }
  }
  return NULL;
}



void Strategy::linkStrategyXML()
{
  int i = 0;

  IN("Strategy::linkStrategyXML\n");

  this->timeDiscretisation = new TimeDiscretisation();
  this->timeDiscretisation->setStrategy(this);
  this->timeDiscretisation->createTimeDiscretisation(this->strategyxml->getTimeDiscretisationXML());

  /*
   * for the moment, all the OneStepIntegrators have the same TimeDiscretisation
   * so, each Integrator is built with this piece of information
   */
  int dsNb;
  DynamicalSystem *dsPtr;
  OneStepIntegrator *integrator;

  // get all the OneStepIntegratorXML objects then create the OneStepIntegrator for this OneStepIntegratorXML and add this OneStepIntergator to the vector of OneStepIntergator of the Strategy
  vector<OneStepIntegratorXML*> osiXMLVector = this->strategyxml->getOneStepIntegratorXML();
  for (i = 0; i < osiXMLVector.size(); i++)
  {
    // with the data of the XML object, we know the type of OneStepIntegrator, so we can instanciate the right type of OneStepIntegrator

    /*
     *  make the link with the DynamicalSystem which must be integrated by this Integrator
     */
    // we get the number of the DynamicalSystem to link
    dsNb = (osiXMLVector[i]->getDSConcerned())[0];

    // we get the address of this DynamicalSystem
    vector<DynamicalSystem*> vds = this->model->getNonSmoothDynamicalSystem()->getDynamicalSystems();

    dsPtr = this->model->getNonSmoothDynamicalSystem()->getDynamicalSystemOnNumber(dsNb);
    if (dsPtr == NULL)
      RuntimeException::selfThrow("Strategy::linkStrategyXML - dsPtr NULL");
    // OneStepIntegrator - Moreau
    if (osiXMLVector[i]->getType() == MOREAU_TAG)
    {
      // creation of the Moreau OneStepIntegrator with this constructor and call of a method to fill
      integrator = new Moreau();
      static_cast<Moreau*>(integrator)->createOneStepIntegrator(osiXMLVector[i], this->timeDiscretisation, dsPtr);
      this->integratorVector.push_back(integrator);
    }
    // OneStepIntegrator - Lsodar
    else if (osiXMLVector[i]->getType() == LSODAR_TAG)
    {
      integrator = new Lsodar();
      static_cast<Lsodar*>(integrator)->createOneStepIntegrator(osiXMLVector[i], this->timeDiscretisation, dsPtr);
      this->integratorVector.push_back(integrator);
    }
    // OneStepIntegrator - Adams
    else if (osiXMLVector[i]->getType() == ADAMS_TAG)
    {
      integrator = new Adams();
      static_cast<Adams*>(integrator)->createOneStepIntegrator(osiXMLVector[i], this->timeDiscretisation, dsPtr);
      this->integratorVector.push_back(integrator);
    }
    else RuntimeException::selfThrow("Strategy::linkStrategyXML - bad kind of Integrator");

    // for the other OneStepIntegrator, we must have the xxxXML.h .cpp objects needed, and the objects in the OneStepIntegrators inherited of the platform
    /*
     *  other "if" to create other Integrators
     *
     */
  }

  // get the OneStepNSProblemXML object then create the OneStepNSProblem for this OneStepNSProblemXML
  // with the data of the XML object, we know the type of OneStepNSProblem, so we can instanciate the right type of OneStepNSProblem

  /*
   *  make the link with the DynamicalSystem which must be integrated by this Integrator
   */
  if (this->strategyxml->hasOneStepNSProblemXML())
  {
    // we get all the numbers of the Interactions to link
    vector<int> interactionNumbers = this->strategyxml->getOneStepNSProblemXML()->getInteractionConcerned();

    // OneStepNSProblem - LCP
    if (this->strategyxml->getOneStepNSProblemXML()->getType() == LCP_TAG)
    {
      // creation of the LCP OneStepNSProblem with this constructor and call of a method to fill
      this->nsProblem = new LCP();
      static_cast<LCP*>(this->nsProblem)->createOneStepNSProblem(this->strategyxml->getOneStepNSProblemXML());
      for (i = 0; i < interactionNumbers.size(); i++)
        this->nsProblem->addInteraction(this->model->getNonSmoothDynamicalSystem()->getInteractionOnNumber(interactionNumbers[i]));
      this->nsProblem->setStrategy(this);
    }
    // OneStepNSProblem - QP
    else if (this->strategyxml->getOneStepNSProblemXML()->getType() == QP_TAG)
    {
      // creation of the QP OneStepNSProblem with this constructor and call of a method to fill
      this->nsProblem = new QP();
      static_cast<QP*>(this->nsProblem)->createOneStepNSProblem(this->strategyxml->getOneStepNSProblemXML());
      for (i = 0; i < interactionNumbers.size(); i++)
        this->nsProblem->addInteraction(this->model->getNonSmoothDynamicalSystem()->getInteractionOnNumber(interactionNumbers[i]));
      this->nsProblem->setStrategy(this);
    }

    // OneStepNSProblem - Relay
    else if (this->strategyxml->getOneStepNSProblemXML()->getType() == RELAY_TAG)
    {
      // creation of the Relay OneStepNSProblem with this constructor and call of a method to fill
      this->nsProblem = new Relay();
      static_cast<Relay*>(this->nsProblem)->createOneStepNSProblem(this->strategyxml->getOneStepNSProblemXML());
      for (i = 0; i < interactionNumbers.size(); i++)
        this->nsProblem->addInteraction(this->model->getNonSmoothDynamicalSystem()->getInteractionOnNumber(interactionNumbers[i]));
      this->nsProblem->setStrategy(this);

    }
    else RuntimeException::selfThrow("Strategy::LinkStrategyXML - bad kind of NSProblem");
  }
  else cout << "Warning : There's no OneStepNSProblem defined, this is an optional attribute." << endl;

  OUT("Strategy::linkStrategyXML\n");
}

void Strategy::fillStrategyWithStrategyXML()
{
  IN("Strategy::fillStrategyWithStrategyXML\n");
  if (this->strategyxml != NULL)
  {}
  else RuntimeException::selfThrow("Strategy::fillStrategyWithStrategyXML - StrategyXML object not exists");
  OUT("Strategy::fillStrategyWithStrategyXML\n");
}

void Strategy::saveStrategyToXML()
{
  IN("Strategy::saveStrategyToXML\n");
  if (this->strategyxml != NULL)
  {
    int size, i;

    size = this->integratorVector.size();
    for (i = 0; i < size; i++)
    {
      if (this->integratorVector[i]->getType() == MOREAU_INTEGRATOR)
        (static_cast<Moreau*>(this->integratorVector[i]))->saveIntegratorToXML();
      else if (this->integratorVector[i]->getType() == ADAMS_INTEGRATOR)
        (static_cast<Adams*>(this->integratorVector[i]))->saveIntegratorToXML();
      else if (this->integratorVector[i]->getType() == LSODAR_INTEGRATOR)
        (static_cast<Lsodar*>(this->integratorVector[i]))->saveIntegratorToXML();
      else RuntimeException::selfThrow("Model::saveToXML - bad kind of OneStepIntegrator");
    }

    if (this->getStrategyXML()->hasOneStepNSProblemXML())
    {
      if (this->nsProblem->getType() == LCP_OSNSP)
        (static_cast<LCP*>(this->nsProblem))->saveNSProblemToXML();
      else if (this->nsProblem->getType() == QP_OSNSP)
        (static_cast<QP*>(this->nsProblem))->saveNSProblemToXML();
      else if (this->nsProblem->getType() == RELAY_OSNSP)
        (static_cast<Relay*>(this->nsProblem))->saveNSProblemToXML();
      else RuntimeException::selfThrow("Model::saveToXML - bad kind of OneStepNSProblem");
    }
  }
  else RuntimeException::selfThrow("Strategy::saveStrategyToXML - StrategyXML object not exists");
  OUT("Strategy::saveStrategyToXML\n");
}

TimeDiscretisation* Strategy::createTimeDiscretisation(double h, int N, SimpleVector * tk,
    double hMin, double hMax, bool constant)
{
  this->timeDiscretisation = new TimeDiscretisation();
  this->timeDiscretisation->createTimeDiscretisation(NULL, h, N, tk, hMin, hMax, constant, this);
  return this->timeDiscretisation;
}

//=========================================================
OneStepNSProblem* Strategy::createLCP()
{
  this->nsProblem = new LCP();
  static_cast<LCP*>(this->nsProblem)->createOneStepNSProblem(NULL, this);
  return this->nsProblem;
}

OneStepNSProblem* Strategy::createQP()
{
  this->nsProblem = new QP();
  static_cast<QP*>(this->nsProblem)->createOneStepNSProblem(NULL, this);
  return this->nsProblem;
}

OneStepNSProblem* Strategy::createRelay()
{
  this->nsProblem = new Relay();
  static_cast<Relay*>(this->nsProblem)->createOneStepNSProblem(NULL, this);
  return this->nsProblem;
}

OneStepIntegrator* Strategy::addAdams(TimeDiscretisation* td, DynamicalSystem* ds)
{
  if (!this->hasDynamicalSystemIntegrator(ds))
  {
    OneStepIntegrator* osi;
    osi = new Adams();
    static_cast<Adams*>(osi)->createOneStepIntegrator(NULL, td, ds);//, this);
    this->integratorVector.push_back(osi);
    return osi;
  }
  else RuntimeException::selfThrow("Strategy::addAdams : Error - The DynamicalSystem of this OneStepIntegrator has already an integrator.");
}

OneStepIntegrator* Strategy::addMoreau(TimeDiscretisation* td, DynamicalSystem* ds,
                                       /*int r,*/ double theta)
{
  if (!this->hasDynamicalSystemIntegrator(ds))
  {
    OneStepIntegrator* osi;
    osi = new Moreau();
    static_cast<Moreau*>(osi)->createOneStepIntegrator(NULL, td, ds/*, r*/, theta);//, this);
    this->integratorVector.push_back(osi);
    return osi;
  }
  else RuntimeException::selfThrow("Strategy::addAdams : Error - The DynamicalSystem of this OneStepIntegrator has already an integrator.");
}

OneStepIntegrator* Strategy::addLsodar(TimeDiscretisation* td, DynamicalSystem* ds)
{
  if (!this->hasDynamicalSystemIntegrator(ds))
  {
    OneStepIntegrator* osi;
    osi = new Lsodar();
    static_cast<Lsodar*>(osi)->createOneStepIntegrator(NULL, td, ds);//, this);
    this->integratorVector.push_back(osi);
    return osi;
  }
  else RuntimeException::selfThrow("Strategy::addAdams : Error - The DynamicalSystem of this OneStepIntegrator has already an integrator.");
}

bool Strategy::hasDynamicalSystemIntegrator(DynamicalSystem* ds)
{
  for (int i = 0; i < integratorVector.size(); i++)
  {
    if (ds == integratorVector[i]->getDynamicalSystem()) return true;
  }
  return false;
}

//$Log: Strategy.cpp,v $
//Revision 1.52  2005/03/11 15:06:20  jbarbier
//- save to XML methods of EqualityConstraint and DSInputOutput added
//
//- XML loading process modified : Model loads NSDS, then NSDS loads the DynamicalSystems, EqualityConstraints, Interactions; Modle loads Strategy, then Strategy loads TimeDiscretisation, then the Integrators, then the OneStepNSProblem
//
//Revision 1.51  2005/03/08 14:23:44  jbarbier
//- modification of constant variables :
//in the XML module, main tags of the XML objects of the strategy are in XMLTagsName.h
//
//in simualtion tools, some constants have been moved to SiconosConst.h
//
//Revision 1.50  2005/03/02 16:06:34  jbarbier
//- DoubleContact sample runnig successfully!
//
//- computeM and computeQ of LCP fixed
//
//Revision 1.49  2005/03/01 15:53:09  jbarbier
//- new sample in progress : 3 balls with 2 balls which are going to touch th third ball at the same time
//
//Revision 1.48  2005/03/01 10:38:20  jbarbier
//- RollingBalls sample is OK
//
//Revision 1.47  2005/02/14 09:52:21  charlety
//_ getters / setters put inline
//
//Revision 1.46  2005/02/04 14:52:44  jbarbier
//- Rolling balls in progress (contact is detected)
//
//- time data is given in parameter for computInput and Output in the Relation. Parameter is sent by methods of the OneStepNSProblem
//
//Revision 1.45  2005/02/04 07:46:21  jbarbier
//- last modification for RollingBalls
//
//Revision 1.44  2005/02/01 11:08:42  charlety
//
//_ some displays of values during computations suppressed.
//
//Revision 1.43  2005/01/20 14:44:49  jbarbier
//- NSDS class renamed NonSmoothDynamicalSystem
//
//- code reduce, some comments remove
//
//Revision 1.42  2005/01/14 09:02:32  jbarbier
//- renaming function
//
//Revision 1.41  2004/12/06 10:10:34  jbarbier
//- integration of Numerics and use of Numerics on the bouncing ball sample
//
//- Numerics is now needed to run the bouncing ball sample!
//
//Revision 1.40  2004/09/27 13:27:14  jbarbier
//
//- Siconos schema renamed : SiconosModelSchema-V1.0.xsd
//
//- new required tags of the model : title, author, description, date, xmlSchema.
//They replace previous attributes author, description and date of the Model.
//
//Revision 1.39  2004/09/23 14:09:24  jbarbier
//- modification of the integrators, the attribute r is always optional.
//
//- modification of the LagrangianNonLinearR. computeInput and computeOutput are
//required.
//
//Revision 1.38  2004/09/16 11:35:25  jbarbier
//- save of the TimeDiscretisation in a XML file in manual creation of the
//platform which was forgotten is now available.
//
//- the save of the platform's data can be done when the platform is created with
//an XML input file and completed with dynmical systems, interactions, one-step
//non smooth problem and one-step integrator.
//
//Revision 1.37  2004/09/15 13:23:13  jbarbier
//- corrections in the OneStepNSProblem, for the XML save. The list of interaction
//linked to the onestepnsproblem is now saved correctly. It is updated before
//during the creation process.
//
//Revision 1.36  2004/09/14 13:49:55  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.35  2004/09/10 11:26:17  charlety
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
//Revision 1.34  2004/09/09 08:57:44  jbarbier
//- functions createLCP, createQP, createRelay, addMoreau, addAdams, addLsodar,
//createTimeDiscretisation of the Strategy done.
//
//=> all functions to create manually the objects of the platform are done
//
//Revision 1.33  2004/08/18 14:37:19  jbarbier
//- creation of Model, NSDS, Strategy(TimeStepping and EventDriven) and
//DynamicalSystem available when the creation is in a command program
//
//Revision 1.32  2004/08/12 11:55:19  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
//Revision 1.31  2004/08/05 14:33:21  charlety
//
//_ class LSODAR is now named Lsodar
//
//Revision 1.30  2004/08/05 12:44:43  jbarbier
//- loading XML file with no OneStepNSProblem succesfull
//
//- NonLinearSystemDS is now available
//
//Revision 1.29  2004/08/03 12:07:11  jbarbier
//- all test on th eModel are successfull
//
//- new tests on the Model with the opening of XML file
//
//- link TimeDiscretisation -> Strategy
//
//- attribute T of the Model is now optional
//
//Revision 1.28  2004/07/29 14:25:40  jbarbier
//- $Log: Strategy.cpp,v $
//- Revision 1.52  2005/03/11 15:06:20  jbarbier
//- - save to XML methods of EqualityConstraint and DSInputOutput added
//-
//- - XML loading process modified : Model loads NSDS, then NSDS loads the DynamicalSystems, EqualityConstraints, Interactions; Modle loads Strategy, then Strategy loads TimeDiscretisation, then the Integrators, then the OneStepNSProblem
//-
//- Revision 1.51  2005/03/08 14:23:44  jbarbier
//- - modification of constant variables :
//- in the XML module, main tags of the XML objects of the strategy are in XMLTagsName.h
//-
//- in simualtion tools, some constants have been moved to SiconosConst.h
//-
//- Revision 1.50  2005/03/02 16:06:34  jbarbier
//- - DoubleContact sample runnig successfully!
//-
//- - computeM and computeQ of LCP fixed
//-
//- Revision 1.49  2005/03/01 15:53:09  jbarbier
//- - new sample in progress : 3 balls with 2 balls which are going to touch th third ball at the same time
//-
//- Revision 1.48  2005/03/01 10:38:20  jbarbier
//- - RollingBalls sample is OK
//-
//- Revision 1.47  2005/02/14 09:52:21  charlety
//- _ getters / setters put inline
//-
//- Revision 1.46  2005/02/04 14:52:44  jbarbier
//- - Rolling balls in progress (contact is detected)
//-
//- - time data is given in parameter for computInput and Output in the Relation. Parameter is sent by methods of the OneStepNSProblem
//-
//- Revision 1.45  2005/02/04 07:46:21  jbarbier
//- - last modification for RollingBalls
//-
//- Revision 1.44  2005/02/01 11:08:42  charlety
//-
//- _ some displays of values during computations suppressed.
//-
//- Revision 1.43  2005/01/20 14:44:49  jbarbier
//- - NSDS class renamed NonSmoothDynamicalSystem
//-
//- - code reduce, some comments remove
//-
//- Revision 1.42  2005/01/14 09:02:32  jbarbier
//- - renaming function
//-
//- Revision 1.41  2004/12/06 10:10:34  jbarbier
//- - integration of Numerics and use of Numerics on the bouncing ball sample
//-
//- - Numerics is now needed to run the bouncing ball sample!
//-
//- Revision 1.40  2004/09/27 13:27:14  jbarbier
//-
//- - Siconos schema renamed : SiconosModelSchema-V1.0.xsd
//-
//- - new required tags of the model : title, author, description, date, xmlSchema.
//- They replace previous attributes author, description and date of the Model.
//-
//- Revision 1.39  2004/09/23 14:09:24  jbarbier
//- - modification of the integrators, the attribute r is always optional.
//-
//- - modification of the LagrangianNonLinearR. computeInput and computeOutput are
//- required.
//-
//- Revision 1.38  2004/09/16 11:35:25  jbarbier
//- - save of the TimeDiscretisation in a XML file in manual creation of the
//- platform which was forgotten is now available.
//-
//- - the save of the platform's data can be done when the platform is created with
//- an XML input file and completed with dynmical systems, interactions, one-step
//- non smooth problem and one-step integrator.
//-
//- Revision 1.37  2004/09/15 13:23:13  jbarbier
//- - corrections in the OneStepNSProblem, for the XML save. The list of interaction
//- linked to the onestepnsproblem is now saved correctly. It is updated before
//- during the creation process.
//-
//- Revision 1.36  2004/09/14 13:49:55  jbarbier
//- - files added in sample/ to run run the main_siconos test program
//-
//- - all the platform can now be saved in an XML file when it is created manually
//-
//- Revision 1.35  2004/09/10 11:26:17  charlety
//-
//- _ Integration of the new version of the SiconosVector in the platform. the class simpleVector is used mostly to replace old SiconosVector. When a vector can be composite or simple, like the state of a dynamical system, a pointer on SiconosVector is used, and the vector is initialized simple or composite when the system is initialized.
//-
//- _ All the tests which worked with the previous version of the vector are OK with the new version.
//-
//- _ Example SICONOS and bouncingBall are OK
//-
//- _ some comments have still to be adapted to NewSiconosVector .
//-
//- _ This version of NewSiconosVector could be called 0.9. some details have to be fixed, it will be done before the end of September.
//-
//- Revision 1.34  2004/09/09 08:57:44  jbarbier
//- - functions createLCP, createQP, createRelay, addMoreau, addAdams, addLsodar,
//- createTimeDiscretisation of the Strategy done.
//-
//- => all functions to create manually the objects of the platform are done
//-
//- Revision 1.33  2004/08/18 14:37:19  jbarbier
//- - creation of Model, NSDS, Strategy(TimeStepping and EventDriven) and
//- DynamicalSystem available when the creation is in a command program
//-
//- Revision 1.32  2004/08/12 11:55:19  jbarbier
//- - new methods createModel, createNSDS, createStrategy, ...
//- they now allow to make the link with upper objects of the platform
//- it will be used for the creation of the platform without XML input file
//-
//- - the createModel method is finished but the attributes of the other objects
//- of the platform are missing for the conctruction
//-
//- Revision 1.31  2004/08/05 14:33:21  charlety
//-
//- _ class LSODAR is now named Lsodar
//-
//- Revision 1.30  2004/08/05 12:44:43  jbarbier
//- - loading XML file with no OneStepNSProblem succesfull
//-
//- - NonLinearSystemDS is now available
//-
//- Revision 1.29  2004/08/03 12:07:11  jbarbier
//- - all test on th eModel are successfull
//-
//- - new tests on the Model with the opening of XML file
//-
//- - link TimeDiscretisation -> Strategy
//-
//- - attribute T of the Model is now optional
//- and $Id: Strategy.cpp,v 1.52 2005/03/11 15:06:20 jbarbier Exp $ added
//
