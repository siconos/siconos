#ifndef STRATEGY_H
#define STRATEGY_H

#include "OneStepIntegrator.h"

#include "OneStepNSProblem.h"
#include "TimeDiscretisation.h"
#include <iostream>
#include <vector>
#include "SiconosMatrix.h"
#include "NewSiconosVector.h"
#include "StrategyXML.h"
#include "SiconosConst.h"
#include "Model.h"


//using namespace std;

class Model;
class OneStepIntegrator;

class OneStepNSProblem;

class TimeDiscretisation;
class StrategyXML;

/** \class Strategy
 *  \brief It regroups all the elements to lead the resolution of the simulation
*  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Crestion) Apr 26, 2004
 *
 *
 */
class Strategy
{
public:
  /** \fn Strategy()
   *  \brief default constructor
   */
  Strategy();

  /** \fn Strategy()
   *  \brief constructor by copy
   */
  Strategy(Strategy*);

  /** \fn Strategy(StrategyXML*, Model*)
   *  \brief constructor with XML object of the Strategy
   *  \param StrategyXML* : the XML object corresponding
   *  \param Model* : the Model which contains the Strategy
   */
  Strategy(StrategyXML*, Model*);

  virtual ~Strategy();

  // getter/setter

  /** \fn Model* getModel()
   *  \brief allows to get the Model which contains the Strategy
   *  \return Model* : the Model which the Strategy
   */
  inline Model* getModel() const
  {
    return model;
  }

  /** \fn TimeDiscretisation* getTimeDiscretisationPtr()
   *  \brief allows to get the TimeDiscretisation of the Strategy
   *  \return the TimeDiscretisation
   */
  inline TimeDiscretisation* getTimeDiscretisationPtr() const
  {
    return timeDiscretisation;
  };

  /** \fn inline int getOneStepIntegratorVectorSize()
   *  \brief allows to get the size of the vector of Integrators
   *  \return int : the size of integratorVector
   */
  inline int getOneStepIntegratorVectorSize() const
  {
    return this->integratorVector.size();
  }

  /** \fn vector<OneStepIntegrator*> getOneStepIntegrators(void)
   *  \brief allows to get all the Integrators of the Strategy
   *  \return a vector of OneStepIntegrator*
   *  \exception RuntimeException
   */
  inline vector<OneStepIntegrator*> getOneStepIntegrators(void) const
  {
    return this->integratorVector;
  };

  /** \fn OneStepIntegrator* getOneStepIntegrator(int)
   *  \brief allows to get one Integrator of the Strategy
   *  \return one OneStepIntegrator* if it exists, else false
   */
  OneStepIntegrator* getOneStepIntegrator(int) const;

  /** \fn OneStepNSProblem* getOneStepNSProblem(void)
   *  \brief allows to get the OneStepNSProblem of the Strategy
   *  \return the OneStepNSProblem of the Strategy if it exists, else false
   */
  inline OneStepNSProblem* getOneStepNSProblem(void) const
  {
    return this->nsProblem;
  };


  /** \fn void setModel(Model* m)
   *  \brief allows to set the Model which contains the Strategy
   *  \param Model* : the Model to set
   */
  inline void setModel(Model* m)
  {
    this->model = m;
  }

  /** \fn void setTimeDiscretisation(TimeDiscretisation*)
   *  \brief allows to set timeDiscretisation
   *  \param the TimeDiscretisation to set
   */
  void setTimeDiscretisation(TimeDiscretisation* td)
  {
    this->timeDiscretisation = td;
  };

  /** \fn void setOneStepIntegrators(vector<OneStepIntegrator*>)
   *  \brief allows to set the vector of Integrator
   *  \param a vector of OneStepIntegrator to set
   */
  void setOneStepIntegrators(const vector<OneStepIntegrator*> vOSI)
  {
    this->integratorVector = vOSI;
  };

  /** \fn void addOneStepIntegrator(OneStepIntegrator*)
   *  \brief allows to add an Integrator to the vector of Integrator
   *  \param the OneStepIntegrator to add
   */
  void addOneStepIntegrator(OneStepIntegrator *osi)
  {
    this->integratorVector.push_back(osi);
  };

  /** \fn void setOneStepNSProblem(OneStepNSProblem*)
   *  \brief allows to set the OneStepNSProblem of the Strategy
   *  \param the OneStepNSProblem to set
   */
  void setOneStepNSProblem(OneStepNSProblem* nspb)
  {
    this->nsProblem = nspb;
  };

  /** \fn inline StrategyXML* getStrategyXML()
   *  \brief allows to get the StrategyXML* of the Strategy
   *  \return a pointer on the StrategyXML of the Strategy
   */
  inline StrategyXML* getStrategyXML() const
  {
    return this->strategyxml;
  }

  /** \fn inline setStrategyXML(StrategyXML* strxml)
   *  \brief allows to set the StrategyXML of the Strategy
   *  \param StrategyXML* : the pointer to set the StrategyXML
   */
  inline void setStrategyXML(StrategyXML* strxml)
  {
    this->strategyxml = strxml;
  }



  /** \fn inline string getType()
   *  \brief allows to get the type of the Strategy
   *  \return string : the type of the Strategy
   */
  inline string getType() const
  {
    return this->strategyType;
  }


  //////////////////////////

  /** \fn void computeFreeState(void)
   *  \brief integrates all the DynamicalSystem without taking care of the relation and NS Laws
   */
  virtual void computeFreeState(void);

  /** \fn void nextStep(void)
   *  \brief increments all the Integrators to go to the next step of the simulation
   */
  virtual void nextStep(void);


  /** \fn void formaliseOneStepNSProblem()
   *  \brief transforms the discretised problem to a numerical problem
   */
  virtual void formaliseOneStepNSProblem();

  /** \fn void predictAndFormalise(void)
   *  \brief begins the computations of the one step NS problem
   */
  virtual void computeOneStepNSProblem(void);


  /** \fn voir updateState()
   *  \brief updates the state of each DynamicalSystem
   */
  virtual void updateState();

  /** \fn void initialize()
   *  \brief executes the complete initialisation of Strategy (OneStepIntegrators, OneStepNSProblem, TImediscretisation) with the XML Object
   */
  virtual void initialize();

  /* \fn OneStepIntegrator* getIntegratorOfDS(int numberDS);
   * \brief searchs the integrator of the DS number "numberDS"
   *
   */
  OneStepIntegrator* getIntegratorOfDS(int numberDS);

  /** \fn void linkStrategyXML()
   *  \brief makes the links between the OneStepIntegratorXMLs, OneStepNSProblemXML, TimeDiscretisationXML of the StrategyXML of the Strategy and the OneStepIntegrators, OneStepNSProblem, TimeDiscretisation
   *  \return a pointer on the integrator of the DS, or NULL if not found
   */
  virtual void linkStrategyXML();

  /** \fn void fillStrategyWithStrategyXML()
   *  \brief uses the StrategyXML of the Strategy to fill the fields of this Strategy
   *  \exception RuntimeException
   */
  virtual void fillStrategyWithStrategyXML();

  /** \fn void saveStrategyToXML()
   *  \brief copys the data of the Strategy to the XML tree
   *  \exception RuntimeException
   */
  virtual void saveStrategyToXML();

  /** \fn void createTimeDiscretisation(double h, int N, SiconosVector * tk,
                double hMin, double hMax, bool constant)
   *  \brief allows to create the TimeDiscretisation of the Strategy
   *  \param double : the h value
   *  \param int : the N value
   *  \param SiconosVector* : the tk vector
   *  \param double : the hMin value
   *  \param double : the hMax value
   *  \param bool : the boolean which determines if the TimeDiscretisation is constant
   *  \return TimeDiscretisation* : the TimeDiscretisation created
   *  \exception RuntimeException
   */
  TimeDiscretisation* createTimeDiscretisation(double h, int N, SimpleVector * tk,
      double hMin, double hMax, bool constant);


  //===========================================
  /** \fn OneStepNSProblem* createLCP()
   *  \brief allows to create a LCP
   *  \return OneStepNSProblem* : the OneStepNSProblem created
   */
  OneStepNSProblem* createLCP();

  /** \fn OneStepNSProblem* createQP()
   *  \brief allows to create a QP
   *  \return OneStepNSProblem* : the OneStepNSProblem created
   */
  OneStepNSProblem* createQP();

  /** \fn OneStepNSProblem* createRelay()
   *  \brief allows to create a Relay NSProblem
   *  \return OneStepNSProblem* : the OneStepNSProblem created
   */
  OneStepNSProblem* createRelay();


  /** \fn OneStepIntegrator* addAdams(TimeDiscretisation* td, DynamicalSystem* ds)
   *  \brief allows to add an Adams integrator to the Strategy
   *  \param TimeDiscretisation* : the TimeDiscretisation of the OneStepIntegrator
   *  \param DynamicalSystem* : the DynamicalSystem that OneStepIntegrator must integrate
   *  \return OneStepIntegrator* : the OneStepIntegrator created
   */
  OneStepIntegrator* addAdams(TimeDiscretisation* td, DynamicalSystem* ds);

  /** \fn OneStepIntegrator* addMoreau(TimeDiscretisation* td, DynamicalSystem* ds,
   *                    int r, double theta)
   *  \brief allows to add an Moreau integrator to the Strategy
   *  \param TimeDiscretisation* : the TimeDiscretisation of the OneStepIntegrator
   *  \param DynamicalSystem* : the DynamicalSystem that OneStepIntegrator must integrate
   *  \param double : the theta value
   *  \return OneStepIntegrator* : the OneStepIntegrator created
   */
  OneStepIntegrator* addMoreau(TimeDiscretisation* td, DynamicalSystem* ds, double theta);

  /** \fn OneStepIntegrator* addLsodar(TimeDiscretisation* td, DynamicalSystem* ds)
   *  \brief allows to add an Lsodar integrator to the Strategy
   *  \param TimeDiscretisation* : the TimeDiscretisation of the OneStepIntegrator
   *  \param DynamicalSystem* : the DynamicalSystem that OneStepIntegrator must integrate
   *  \return OneStepIntegrator* : the OneStepIntegrator created
   */
  OneStepIntegrator* addLsodar(TimeDiscretisation* td, DynamicalSystem* ds);

  /** \fn bool hasDynamicalSystemIntegrator( DynamicalSystem* ds)
   *  \brief checks if a DynamicalSystem owns already an OneStepIntegrator
   *  \return bool : false if the DynamicalSystem has no OneStepIntegrator, else true
   */
  bool hasDynamicalSystemIntegrator(DynamicalSystem* ds);


  /** \fn void newtonSolve(double criterion , int maxStep)
   *  \brief newton algorithm
   *  \param double criterion: convergence criterion, int maxStep: maximum number of Newton steps
   */
  void newtonSolve(double criterion , int maxStep);

  /** */
  void Strategy::newtonNextStep();


  /** \fn newtonUpdateState()
   *  \brief update the state of the dynamical system at the end of Newton step
   */

  void Strategy::newtonUpdateState();


  /** \fn newtonCheckConvergence(double criterion);
   *  \brief check the convergence of Newton algorithm
   */
  bool Strategy::newtonCheckConvergence(double criterion);


protected:
  /** the name of the type of the Strategy */
  string strategyType;

  /** contains the data of the time discretisation */
  TimeDiscretisation *timeDiscretisation;
  //  vector<TimeDiscretisation*> timeDiscretisations;

  /** contains the integrators to integre the DS */
  vector<OneStepIntegrator*> integratorVector;

  /** contains the type of resolution (on en a 0..1 ou 1 ??) */
  OneStepNSProblem *nsProblem;

  /** the XML object linked to the Strategy to read XML data */
  StrategyXML *strategyxml;

  /** the Model which contains the Strategy */
  Model *model;
  //  NSDS *nsds;

};

#endif // STRATEGY_H
