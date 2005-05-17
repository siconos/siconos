
#include "Strategy.h"
#include "Moreau.h"
#include "Lsodar.h"
#include "Adams.h"
#include "LagrangianLinearTIDS.h"

#include "LCP.h"
#include "CFD.h"
#include "QP.h"
#include "Relay.h"

#include "check.h"

Strategy::Strategy()
{
  IN("Strategy::Strategy()\n");
  this->timeDiscretisation = 0;
  this->integratorVector.clear();
  this->nsProblem = 0;

  this->strategyxml = 0;
  this->model = 0;
  OUT("Strategy::Strategy()\n");
}

Strategy::Strategy(Strategy* str)
{
  IN("Strategy::Strategy( Strategy* str )\n");
  this->strategyType = str->getType();
  this->timeDiscretisation = str->getTimeDiscretisationPtr();
  this->integratorVector = str->getOneStepIntegrators();
  this->nsProblem = str->getOneStepNSProblem();
  this->strategyxml = str->getStrategyXML();
  this->model = str->getModel();
  OUT("Strategy::Strategy( Strategy* str )\n");
}

Strategy::Strategy(StrategyXML* strxml, Model *model)
{
  this->timeDiscretisation = 0;
  this->integratorVector.clear();
  this->nsProblem = 0;

  this->strategyxml = strxml;
  this->model = model;
}

Strategy::~Strategy()
{
  if (this->nsProblem != 0) delete this->nsProblem;
  if (this->integratorVector.size() > 0)
  {
    for (int i = 0; i < this->integratorVector.size(); i++)
    {
      delete this->integratorVector[i];
    }
    this->integratorVector.clear();
  }
  if (this->timeDiscretisation != 0)
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
  if (this->nsProblem != 0)
    this->nsProblem->nextStep();

}


void Strategy::formaliseOneStepNSProblem()
{
  // formalise the OneStepNSProblem
  if (this->nsProblem != 0)this->nsProblem->formalize(this->model->getCurrentT());
}

void Strategy::computeOneStepNSProblem(void)
{
  // compute the OneStepNSProblem
  if (this->nsProblem != 0)this->nsProblem->compute();
}

void Strategy::updateState()
{
  // compute the OneStepNSProblem
  if (this->nsProblem != 0)this->nsProblem->updateState();


  for (int i = 0; i < this->integratorVector.size(); i++)
  {
    this->integratorVector[i]->updateState();
  }

  this->model->setCurrentT(this->model->getCurrentT() + this->timeDiscretisation->getH());
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
  if (this->nsProblem != 0)
    this->nsProblem->initialize();
}



OneStepIntegrator* Strategy::getIntegratorOfDS(int numberDS)
{
  for (int i = 0; i < this->integratorVector.size(); i++)
  {
    if (this->integratorVector[i]->getDynamicalSystemPtr()->getNumber() == numberDS)
    {
      return this->integratorVector[i];
      break;
    }
  }
  return 0;
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
    if (dsPtr == 0)
      RuntimeException::selfThrow("Strategy::linkStrategyXML - dsPtr 0");
    // OneStepIntegrator - Moreau
    if (osiXMLVector[i]->getType() == MOREAU_TAG)
    {
      // creation of the Moreau OneStepIntegrator with this constructor and call of a method to fill
      integrator = new Moreau(osiXMLVector[i], this->timeDiscretisation, dsPtr);
      //static_cast<Moreau*>(integrator)->createOneStepIntegrator( osiXMLVector[i], this->timeDiscretisation, dsPtr );
      this->integratorVector.push_back(integrator);
    }
    // OneStepIntegrator - Lsodar
    else if (osiXMLVector[i]->getType() == LSODAR_TAG)
    {
      integrator = new Lsodar(osiXMLVector[i]);
      //static_cast<Lsodar*>(integrator)->createOneStepIntegrator( osiXMLVector[i], this->timeDiscretisation, dsPtr );
      this->integratorVector.push_back(integrator);
    }
    // OneStepIntegrator - Adams
    else if (osiXMLVector[i]->getType() == ADAMS_TAG)
    {
      integrator = new Adams(osiXMLVector[i]);
      //static_cast<Adams*>(integrator)->createOneStepIntegrator( osiXMLVector[i], this->timeDiscretisation, dsPtr );
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
    // OneStepNSProblem - CFD
    else if (this->strategyxml->getOneStepNSProblemXML()->getType() == CFD_TAG)
    {
      // creation of the CFD OneStepNSProblem with this constructor and call of a method to fill
      this->nsProblem = new CFD();
      static_cast<CFD*>(this->nsProblem)->createOneStepNSProblem(this->strategyxml->getOneStepNSProblemXML());
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
  if (this->strategyxml != 0)
  {}
  else RuntimeException::selfThrow("Strategy::fillStrategyWithStrategyXML - StrategyXML object not exists");
  OUT("Strategy::fillStrategyWithStrategyXML\n");
}

void Strategy::saveStrategyToXML()
{
  IN("Strategy::saveStrategyToXML\n");
  if (this->strategyxml != 0)
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
      else RuntimeException::selfThrow("Strategy::saveStrategyToXML - bad kind of OneStepIntegrator");
    }

    if (this->getStrategyXML()->hasOneStepNSProblemXML())
    {
      if (this->nsProblem->getType() == LCP_OSNSP)
        (static_cast<LCP*>(this->nsProblem))->saveNSProblemToXML();
      else if (this->nsProblem->getType() == CFD_OSNSP)
        (static_cast<CFD*>(this->nsProblem))->saveNSProblemToXML();
      else if (this->nsProblem->getType() == QP_OSNSP)
        (static_cast<QP*>(this->nsProblem))->saveNSProblemToXML();
      else if (this->nsProblem->getType() == RELAY_OSNSP)
        (static_cast<Relay*>(this->nsProblem))->saveNSProblemToXML();
      else
        RuntimeException::selfThrow("Strategy::saveStrategyToXML - bad kind of OneStepNSProblem");
    }
  }
  else RuntimeException::selfThrow("Strategy::saveStrategyToXML - StrategyXML object not exists");
  OUT("Strategy::saveStrategyToXML\n");
}

TimeDiscretisation* Strategy::createTimeDiscretisation(double h, int N, SimpleVector * tk,
    double hMin, double hMax, bool constant)
{
  this->timeDiscretisation = new TimeDiscretisation();
  this->timeDiscretisation->createTimeDiscretisation(0, h, N, tk, hMin, hMax, constant, this);
  return this->timeDiscretisation;
}

//=========================================================
OneStepNSProblem* Strategy::createLCP()
{
  this->nsProblem = new LCP();
  static_cast<LCP*>(this->nsProblem)->createOneStepNSProblem(0, this);
  return this->nsProblem;
}

OneStepNSProblem* Strategy::createQP()
{
  this->nsProblem = new QP();
  static_cast<QP*>(this->nsProblem)->createOneStepNSProblem(0, this);
  return this->nsProblem;
}

OneStepNSProblem* Strategy::createRelay()
{
  this->nsProblem = new Relay();
  static_cast<Relay*>(this->nsProblem)->createOneStepNSProblem(0, this);
  return this->nsProblem;
}

OneStepIntegrator* Strategy::addAdams(TimeDiscretisation* td, DynamicalSystem* ds)
{
  if (!this->hasDynamicalSystemIntegrator(ds))
  {
    OneStepIntegrator* osi;
    osi = new Adams(td, ds);
    this->integratorVector.push_back(osi);
    return osi;
  }
  else RuntimeException::selfThrow("Strategy::addAdams : Error - The DynamicalSystem of this OneStepIntegrator has already an integrator.");
}

OneStepIntegrator* Strategy::addMoreau(TimeDiscretisation* td, DynamicalSystem* ds, double theta)
{
  if (!this->hasDynamicalSystemIntegrator(ds))
  {
    OneStepIntegrator* osi;
    osi = new Moreau(td, ds, theta);
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
    osi = new Lsodar(td, ds);
    this->integratorVector.push_back(osi);
    return osi;
  }
  else RuntimeException::selfThrow("Strategy::addAdams : Error - The DynamicalSystem of this OneStepIntegrator has already an integrator.");
}

bool Strategy::hasDynamicalSystemIntegrator(DynamicalSystem* ds)
{
  for (int i = 0; i < integratorVector.size(); i++)
  {
    if (ds == integratorVector[i]->getDynamicalSystemPtr()) return true;
  }
  return false;
}

void Strategy::newtonSolve(double criterion, int maxStep)
{

  bool isNewtonConverge = false;
  long nbNewtonStep = 0; // number of Newton iterations
  while ((!isNewtonConverge) && (nbNewtonStep <= maxStep))
  {
    //s->newtonNextStep();
    nbNewtonStep++;
    cout << "nbStep:" << nbNewtonStep << endl ;
    computeFreeState();
    formaliseOneStepNSProblem();
    computeOneStepNSProblem();
    newtonUpdateState();
    isNewtonConverge = newtonCheckConvergence(criterion);
  }
  if (isNewtonConverge)
    cout << "Newton process: convergence reached after " << nbNewtonStep << " Newtons steps" << endl;
  else
    cout << "Newton process stopped: reach max step number" << endl ;


  // time step increment
  this->model->setCurrentT(this->model->getCurrentT() + this->timeDiscretisation->getH());
}

void Strategy::newtonNextStep()
{
}


void Strategy::newtonUpdateState()
{
  IN("Strategy::newtonUpdateState\n");
  // update NonSmooth problem
  if (this->nsProblem != 0)this->nsProblem->updateState();

  // update OneStep Integrators
  for (int i = 0; i < this->integratorVector.size(); i++)
  {
    this->integratorVector[i]->updateState();
  }

  OUT("Strategy::newtonUpdateState\n");

}

bool Strategy::newtonCheckConvergence(double criterion)
{
  bool checkConvergence = false;

  // get the non smooth dynamical system
  NonSmoothDynamicalSystem* nsds = this->model-> getNonSmoothDynamicalSystem();

  // get the nsds indicator of convergence
  double nsdsConverge = nsds -> nsdsConvergenceIndicator();
  if (nsdsConverge < criterion) checkConvergence = true ;

  return(checkConvergence);
}


