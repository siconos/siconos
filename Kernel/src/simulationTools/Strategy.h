/* Siconos-Kernel version 1.1.4, Copyright INRIA 2005-2006.
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
#ifndef STRATEGY_H
#define STRATEGY_H

#include "OneStepIntegrator.h"
#include "OneStepNSProblem.h"
#include "TimeDiscretisation.h"
#include "SiconosMatrix.h"
#include "SiconosVector.h"
#include "StrategyXML.h"
#include "SiconosConst.h"
#include "Model.h"

#include "check.h"
#include <iostream>
#include <vector>
#include <deque>
#include <set>

class Model;
class OneStepIntegrator;
class OneStepNSProblem;
class TimeDiscretisation;
class StrategyXML;

/** \class Strategy
 *  \brief It regroups all the elements to lead the resolution of the simulation
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.1.4.
 *  \date (Crestion) Apr 26, 2004
 *
 * !!! This is an abstract class !!!
 *
 * The available strategies are time stepping and event driven. See derived classes for more details.
 *
 * Rules:
 *     - A Model* must be given to the constructor, else => exception.
 *     - Strategy must always be built before TimeDiscretisation, which
 *       means that a TimeDiscretisation constructor can not be called
 *       without a strategy as argument.
 *       Moreover, there is no strategy->setTimeDiscretisation() function
 */

/** vector of OneStepIntegrator */
typedef std::set<OneStepIntegrator*> OSISet;

/** iterator through vector of OSI*/
typedef OSISet::iterator OSIIterator;

/** const iterator through vector of OSI*/
typedef OSISet::const_iterator ConstOSIIterator;

/** return type value for insert function - bool = false if insertion failed. */
typedef std::pair<OSISet::iterator, bool> CheckInsertOSI;

class Strategy
{
protected:

  /** name or id of the Strategy */
  std::string name;

  /** the type of the Strategy */
  std::string strategyType;

  /** the time discretisation scheme */
  TimeDiscretisation *timeDiscretisation;

  /** the dynamical systems integrators */
  OSISet allOSI;

  /** inside-class allocation flags*/
  std::map<OneStepIntegrator*, bool> isOSIAllocatedIn;

  /** the non smooth problem */
  OneStepNSProblem *nsProblem;

  /** the XML object linked to the Strategy to read XML data */
  StrategyXML *strategyxml;

  /** A link to the Model which contains the Strategy */
  Model *model;

  /** Flags to check wheter pointers were allocated in the present class constructors or not */
  bool isTimeDiscretisationAllocatedIn;
  bool isNSProblemAllocatedIn;

public:

  /** \fn Strategy(Model* = NULL, const string& id = "undefined")
   *  \brief default constructor
   *  \param a pointer to the model that owns this strategy. NULL Model leads to exception
   *  \param string: strategy type, default = undefined
   */
  Strategy(Model* = NULL, const std::string& = "undefined");

  /** \fn Strategy(Model&,const std::string& = "undefined")
   *  \brief constructor from Model => avoid this function, prefer the one with Model*
   *  \param a Model.
   *  \param string: strategy type, default = undefined
   */
  Strategy(Model&, const std::string& = "undefined");

  /** \fn Strategy(const OSISet&, OneStepNSProblem*, Model*, const std::string& = "undefined")
   *  \brief constructor from a given set of data (1)
   *  \param a set of osi
   *  \param pointer to a OneStepNSProblem
   *  \param the model that owns this strategy
   *  \param string: strategy type, default = undefined
   */
  Strategy(const OSISet&, OneStepNSProblem *, Model*, const std::string& = "undefined");

  /** \fn Strategy(const OSISet&, Model*, const std::string& = "undefined")
   *  \brief constructor from a given set of data (2)
   *  \param a set of osi
   *  \param the model that owns this strategy
   *   with previous constructor)
   *  \param string: strategy type, default = undefined
   */
  Strategy(const OSISet&, Model*, const std::string& = "undefined");

  /** \fn Strategy(OneStepNSProblem *, Model* , const std::string& = "undefined")
   *  \brief constructor from a given set of data (3)
   *  \param pointer to a OneStepNSProblem
   *  \param the model that owns this strategy
   *  \param string: strategy type, default = undefined
   */
  Strategy(OneStepNSProblem *, Model*, const std::string& = "undefined");

  /** \fn Strategy(StrategyXML*, Model*, const std::string& = "undefined")
   *  \brief constructor with XML object of the Strategy
   *  \param StrategyXML* : the XML object corresponding
   *  \param the model which owns this strategy  (optional parameter)
   *  \param string: strategy type, default = undefined
   */
  Strategy(StrategyXML*, Model*, const std::string& = "undefined");

  /** \fn ~Strategy()
   *  \brief destructor
   */
  virtual ~Strategy();

  /** \fn isStrategyComplete()
   *  \brief test if all strategy members are complete or not.
   */
  bool isStrategyComplete() const;

  // GETTERS/SETTERS

  /** \fn inline string getName() const
   *  \brief get the name of the Strategy
   *  \return string : the name of the Strategy
   */
  inline const std::string getName() const
  {
    return name;
  }

  /** \fn inline string setName(const string&)
   *  \brief set the name of the Strategy
   */
  inline void setName(const std::string& newName)
  {
    name = newName;
  }

  /** \fn inline string getType() const
   *  \brief get the type of the Strategy
   *  \return string : the type of the Strategy
   */
  inline const std::string getType() const
  {
    return strategyType;
  }

  /** \fn inline string setType(const string&)
   *  \brief set the type of the Strategy
   */
  inline void setType(const std::string& newType)
  {
    strategyType = newType;
  }

  /** \fn TimeDiscretisation* getTimeDiscretisationPtr()
   *  \brief get the TimeDiscretisation of the Strategy
   *  \return the TimeDiscretisation
   */
  inline TimeDiscretisation* getTimeDiscretisationPtr() const
  {
    return timeDiscretisation;
  };

  /** \fn void setTimeDiscretisationPtr(TimeDiscretisation*)
   *  \brief set the TimeDiscretisation of the Strategy
   *  \param the TimeDiscretisation
   */
  void setTimeDiscretisationPtr(TimeDiscretisation*);

  /** \fn const OSISet getOneStepIntegrators()
   *  \brief get all the Integrators of the Strategy
   *  \return an OSISset
   */
  inline const OSISet getOneStepIntegrators() const
  {
    return allOSI;
  };

  /** \fn void setOneStepIntegrators(const OSISet&)
   *  \brief set the Integrators of the Strategy
   *  \param an OSISset
   */
  void setOneStepIntegrators(const OSISet&);

  /** \fn inline const unsigned int getNumberOfOSI()
   *  \brief get the number of OSIs in the Strategy (ie the size of allOSI)
   *  \return an unsigned int
   */
  inline const unsigned int getNumberOfOSI() const
  {
    return allOSI.size();
  }

  /** \fn OneStepNSProblem* getOneStepNSProblemPtr(void)
   *  \brief get the OneStepNSProblem of the Strategy
   *  \return the OneStepNSProblem of the Strategy if it exists
   */
  inline OneStepNSProblem* getOneStepNSProblemPtr() const
  {
    return nsProblem;
  };

  /** \fn void setOneStepNSProblemPtr(OneStepNSProblem*)
   *  \brief set the OneStepNSProblem of the Strategy
   *  \param the OneStepNSProblem to set
   */
  void setOneStepNSProblemPtr(OneStepNSProblem* nspb);

  /** \fn inline StrategyXML* getStrategyXMLPtr() const
   *  \brief get the StrategyXML* of the Strategy
   *  \return a pointer on the StrategyXML of the Strategy
   */
  inline StrategyXML* getStrategyXMLPtr() const
  {
    return strategyxml;
  }

  /** \fn inline setStrategyXMLPtr(StrategyXML* strxml)
   *  \brief set the StrategyXML of the Strategy
   *  \param StrategyXML* : the pointer to set the StrategyXML
   */
  inline void setStrategyXMLPtr(StrategyXML* strxml)
  {
    strategyxml = strxml;
  }

  /** \fn Model* getModelPtr()
   *  \brief get the Model which contains the Strategy
   *  \return Model* : the Model which the Strategy
   */
  inline Model* getModelPtr() const
  {
    return model;
  }

  /** \fn void setModelPtr(Model* m)
   *  \brief set the Model which contains the Strategy
   *  \param Model* : the Model to set
   */
  inline void setModelPtr(Model* m)
  {
    model = m;
  }

  /* \fn OneStepIntegrator* getIntegratorOfDSptr(const int& numberDS);
   * \brief searchs the integrator of the DS number "numberDS"
   *
   */
  OneStepIntegrator* getIntegratorOfDSPtr(const int&  numberDS) const ;

  /* \fn OneStepIntegrator* getIntegratorOfDSptr(DynamicalSystem * ds) ;
   * \brief get the integrator of ds
   * \param a pointer to DynamicalSystem
   */
  OneStepIntegrator* getIntegratorOfDSPtr(DynamicalSystem * ds) const ;

  // --- OTHER FUNCTIONS ---

  /** \fn void addOneStepIntegratorPtr(OneStepIntegrator*)
   *  \brief add an Integrator into allOSI (pointer link, no copy!)
   *  \param a pointer to a OneStepIntegrator
   */
  void addOneStepIntegratorPtr(OneStepIntegrator *);

  /** \fn void computeFreeState()
   *  \brief integrates all the DynamicalSystem taking not into account nslaw, reactions ...
   */
  virtual void computeFreeState();

  /** \fn void nextStep()
   *  \brief increments all the Integrators to next step of the simulation
   */
  virtual void nextStep();

  /** \fn void computeOneStepNSProblem()
   *  \brief computes the one step NS problem
   */
  virtual void computeOneStepNSProblem();

  /** \fn void update()
   *  \brief update input, state of each dynamical system and output
   */
  virtual void update();

  /** \fn void initialize()
   *  \brief executes the complete initialisation of Strategy (OneStepIntegrators, OneStepNSProblem, TImediscretisation) with the XML Object
   */
  virtual void initialize() = 0;

  /** \fn void run()
   *  \brief run the simulation, from t0 to T
   */
  virtual void run() = 0;

  /** \fn void computeOneStep()
   *  \brief run one step of the simulation
   */
  virtual void computeOneStep() = 0;

  /** \fn void newtonSolve(const double& criterion , const unsigned int& maxStep)
   *  \brief newton algorithm
   *  \param double criterion: convergence criterion, int maxStep: maximum number of Newton steps
   */
  void newtonSolve(const double&, const unsigned int&);

  /** \fn newtonCheckConvergence(const double& criterion);
   *  \brief check the convergence of Newton algorithm
   */
  bool Strategy::newtonCheckConvergence(const double& criterion);

  // --- XML RELATED FUNCTIONS ---

  /** \fn void saveStrategyToXML()
   *  \brief copys the data of the Strategy to the XML tree
   *  \exception RuntimeException
   */
  virtual void saveStrategyToXML();

  /** \fn bool hasDynamicalSystemIntegrator( DynamicalSystem* ds) const
   *  \brief check if ds is already present in an OSI of the strategy.
   *  \return a bool
   */
  bool hasDynamicalSystemIntegrator(DynamicalSystem*) const ;

  /** \fn bool hasDynamicalSystemIntegrator( OneStepIntegrator* osi) const
   *  \brief check if a DynamicalSystem of osi is already present in another OSI of the strategy.
   *  \return a bool
   */
  bool hasDynamicalSystemIntegrator(OneStepIntegrator*) const ;
};

#endif // STRATEGY_H
