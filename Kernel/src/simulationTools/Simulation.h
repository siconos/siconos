/* Siconos-Kernel version 1.2.0, Copyright INRIA 2005-2006.
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
#ifndef SIMULATION_H
#define SIMULATION_H

#include "OneStepIntegrator.h"
#include "OneStepNSProblem.h"
#include "TimeDiscretisation.h"
#include "SiconosMatrix.h"
#include "SiconosVector.h"
#include "SimulationXML.h"
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
class SimulationXML;

/** \class Simulation
 *  \brief It regroups all the elements to lead the resolution of the simulation
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.2.0.
 *  \date (Crestion) Apr 26, 2004
 *
 * !!! This is an abstract class !!!
 *
 * The available simulations are time stepping and event driven. See derived classes for more details.
 *
 * "Unitary block matrices" are saved in the simulation (using stl maps). They will be used to save block matrices corresponding
 * to diagonal or extra diagonal terms in OneStepNSProblem "complete" matrices, such as M in LCP formulation.
 * Each block corresponds to one UnitaryRelation (diagonal terms) or two UnitaryRelations (extra diagonal terms);
 * Their dimensions are NSLawSize X NSLawSize, NSLawSize being the dimension of the NonSmoothLaw corresponding to the concerned UnitaryRelation.
 *
 * Rules:
 *     - A Model* must be given to the constructor, else => exception.
 *     - Simulation must always be built before TimeDiscretisation, which
 *       means that a TimeDiscretisation constructor can not be called
 *       without a simulation as argument.
 *       Moreover, there is no simulation->setTimeDiscretisation() function
 */

/** vector of OneStepIntegrator */
typedef std::set<OneStepIntegrator*> OSISet;

/** iterator through vector of OSI*/
typedef OSISet::iterator OSIIterator;

/** const iterator through vector of OSI*/
typedef OSISet::const_iterator ConstOSIIterator;

/** return type value for insert function - bool = false if insertion failed. */
typedef std::pair<OSISet::iterator, bool> CheckInsertOSI;

/** vector that contains a sequel of sets of UnitaryRelations*/
typedef std::vector< UnitaryRelationsSet > VectorOfSetOfUnitaryRelations;

/** map of OSNS */
typedef std::map<std::string, OneStepNSProblem * > OneStepNSProblems;

/** and its corresponding iterator */
typedef OneStepNSProblems::iterator OSNSIterator;
typedef OneStepNSProblems::const_iterator ConstOSNSIterator;

/** tolerance value used in indexSets updating */
const double TOLERANCE = 1e-5;

/** default name for One Step NS Problem of the simulation */
const std::string DEFAULT_OSNS_NAME = "unamed";

class Simulation
{
protected:

  /** name or id of the Simulation */
  std::string name;

  /** the type of the Simulation (TimeStepping, EventDriven ...) */
  std::string simulationType;

  /** the time discretisation scheme */
  TimeDiscretisation *timeDiscretisation;

  /** Inside-class allocation flag for time discretisation */
  bool isTimeDiscretisationAllocatedIn;

  /** the dynamical systems integrators */
  OSISet allOSI;

  /** inside-class allocation flags*/
  std::map<OneStepIntegrator*, bool> isOSIAllocatedIn;

  /** index sets vector (indexSets[0] is the set where y[0]=0, indexSets[1] where y[0] = 0 and y[1]=0 and so on */
  VectorOfSetOfUnitaryRelations indexSets;

  /** the non smooth problems (each problem is identified thanks to its id) */
  OneStepNSProblems allNSProblems;

  /** Inside-class allocation flag for Non Smooth Problem(s) */
  std::map<OneStepNSProblem*, bool >  isNSProblemAllocatedIn;

  /** the XML object linked to the Simulation to read XML data */
  SimulationXML *simulationxml;

  /** A link to the Model which contains the Simulation */
  Model *model;

  /** int used to set the minimal derivative order used in the OSNS variables */
  unsigned int levelMin;

  /** int used to set the maximal derivative order used in the OSNS variables */
  unsigned int levelMax;

  /** \fn Simulation()
   *  \brief default constructor
   */
  Simulation(const std::string = "undefined");

public:

  /** \fn Simulation(Model*, const string = "undefined")
   *  \brief default constructor
   *  \param a pointer to the model that owns this simulation. NULL Model leads to exception
   *  \param string: simulation type, default = undefined
   */
  Simulation(Model*, const std::string = "undefined");

  /** \fn Simulation(SimulationXML*, Model*, const std::string = "undefined")
   *  \brief constructor with XML object of the Simulation
   *  \param SimulationXML* : the XML object corresponding
   *  \param the model which owns this simulation  (optional parameter)
   *  \param string: simulation type, default = undefined
   */
  Simulation(SimulationXML*, Model*, const std::string = "undefined");

  /** \fn ~Simulation()
   *  \brief destructor
   */
  virtual ~Simulation();

  // GETTERS/SETTERS

  /** \fn inline string getName() const
   *  \brief get the name of the Simulation
   *  \return string : the name of the Simulation
   */
  inline const std::string getName() const
  {
    return name;
  }

  /** \fn inline string setName(const string)
   *  \brief set the name of the Simulation
   */
  inline void setName(const std::string newName)
  {
    name = newName;
  }

  /** \fn inline string getType() const
   *  \brief get the type of the Simulation
   *  \return string : the type of the Simulation
   */
  inline const std::string getType() const
  {
    return simulationType;
  }

  /** \fn inline string setType(const string)
   *  \brief set the type of the Simulation
   */
  inline void setType(const std::string newType)
  {
    simulationType = newType;
  }

  /** \fn TimeDiscretisation* getTimeDiscretisationPtr()
   *  \brief get the TimeDiscretisation of the Simulation
   *  \return the TimeDiscretisation
   */
  inline TimeDiscretisation* getTimeDiscretisationPtr() const
  {
    return timeDiscretisation;
  };

  /** \fn void setTimeDiscretisationPtr(TimeDiscretisation*)
   *  \brief set the TimeDiscretisation of the Simulation
   *  \param the TimeDiscretisation
   */
  void setTimeDiscretisationPtr(TimeDiscretisation*);

  /** \fn const OSISet getOneStepIntegrators()
   *  \brief get all the Integrators of the Simulation
   *  \return an OSISset
   */
  inline const OSISet getOneStepIntegrators() const
  {
    return allOSI;
  };

  /** \fn void setOneStepIntegrators(const OSISet&)
   *  \brief set the Integrators of the Simulation
   *  \param an OSISset
   */
  void setOneStepIntegrators(const OSISet&);

  /* \fn OneStepIntegrator* getIntegratorOfDSptr(const int numberDS);
   * \brief searchs the integrator of the DS number "numberDS"
   *
   */
  OneStepIntegrator* getIntegratorOfDSPtr(const int numberDS) const ;

  /* \fn OneStepIntegrator* getIntegratorOfDSptr(DynamicalSystem * ds) ;
   * \brief get the integrator of ds
   * \param a pointer to DynamicalSystem
   */
  OneStepIntegrator* getIntegratorOfDSPtr(DynamicalSystem * ds) const ;

  /** \fn inline const unsigned int getNumberOfOSI()
   *  \brief get the number of OSIs in the Simulation (ie the size of allOSI)
   *  \return an unsigned int
   */
  inline const unsigned int getNumberOfOSI() const
  {
    return allOSI.size();
  }

  /** \fn void addOneStepIntegratorPtr(OneStepIntegrator*)
   *  \brief add an Integrator into allOSI (pointer link, no copy!)
   *  \param a pointer to a OneStepIntegrator
   */
  void addOneStepIntegratorPtr(OneStepIntegrator *);

  /** \fn const bool hasCommonDSInIntegrators( OneStepIntegrator* osi) const
   *  \brief check if a DynamicalSystem of osi is already present in another OSI of the simulation.
   *  \return a bool
   */
  const bool hasCommonDSInIntegrators(OneStepIntegrator*) const ;

  /** \fn const VectorOfSetOfUnitaryRelations getIndexSets() const
   *  \brief get indexSets (the whole vector)
   *  \return a VectorOfSetOfUnitaryRelations
   */
  inline const VectorOfSetOfUnitaryRelations getIndexSets() const
  {
    return indexSets;
  };

  /** \fn const UnitaryRelationsSet getIndexSet(const unsigned int i) const
   *  \brief get indexSets[i]
   *  \return a UnitaryRelationsSet
   */
  const UnitaryRelationsSet getIndexSet(const unsigned int) const ;

  /** \fn const OneStepNSProblems getOneStepNSProblems() const
   *  \brief get allNSProblems
   *  \return a OneStepNSProblems object (container of OneStepNSProblem*)
   */
  inline const OneStepNSProblems getOneStepNSProblems() const
  {
    return allNSProblems;
  };

  /** \fn OneStepNSProblem* getOneStepNSProblemPtr(const string name =  DEFAULT_OSNS_NAME) const
   *  \brief get allNSProblems[name], a specific OneStepNSProblem
   *  \param a string, the name of the osns, optional, default value = DEFAULT_OSNS_NAME
   *  \return a pointer to OneStepNSProblem
   */
  OneStepNSProblem* getOneStepNSProblemPtr(const std::string = DEFAULT_OSNS_NAME);

  /** \fn void setOneStepNSProblems(const OneStepNSProblems& allOSNS)
   *  \brief set allNSProblems map - Warning: no copy between OneStepNSProblem of each map, pointers links!
   *  \param a OneStepNSProblems object (map of OneStepNSProblem*)
   */
  void setOneStepNSProblems(const OneStepNSProblems&);

  /** \fn void clearOneStepNSProblems()
  *  \brief remove all OneStepNSProblem of the Simulation
  */
  void clearOneStepNSProblems();

  /** \fn const bool hasOneStepNSProblem(OneStepNSProblem * osns) const
   *  \brief check if a OneStepNSProblem osns is already in the map
   *  \param a pointer to OneStepNSProblem
   *  \return a bool
   */
  const bool hasOneStepNSProblem(OneStepNSProblem*) const ;

  /** \fn const bool hasOneStepNSProblem(const string id) const
   *  \brief check if a OneStepNSProblem named id is already in the map
   *  \param a string
   *  \return a bool
   */
  const bool hasOneStepNSProblem(const std::string) const ;

  /** \fn void addOneStepNSProblemPtr(OneStepNSProblem*)
   *  \brief add a OneStepNSProblem of the Simulation (if its not the first, it needs to have an id clearly defined)
   *  \param a pointer to OneStepNSProblem
   */
  virtual void addOneStepNSProblemPtr(OneStepNSProblem*);

  /** \fn inline SimulationXML* getSimulationXMLPtr() const
   *  \brief get the SimulationXML* of the Simulation
   *  \return a pointer on the SimulationXML of the Simulation
   */
  inline SimulationXML* getSimulationXMLPtr() const
  {
    return simulationxml;
  }

  /** \fn inline setSimulationXMLPtr(SimulationXML* strxml)
   *  \brief set the SimulationXML of the Simulation
   *  \param SimulationXML* : the pointer to set the SimulationXML
   */
  inline void setSimulationXMLPtr(SimulationXML* strxml)
  {
    simulationxml = strxml;
  }

  /** \fn Model* getModelPtr()
   *  \brief get the Model which contains the Simulation
   *  \return Model* : the Model which the Simulation
   */
  inline Model* getModelPtr() const
  {
    return model;
  }

  /** \fn void setModelPtr(Model* m)
   *  \brief set the Model which contains the Simulation
   *  \param Model* : the Model to set
   */
  inline void setModelPtr(Model* m)
  {
    model = m;
  }

  /** \fn virtual void updateIndexSets() = 0;
   *   \brief update all index sets of the topology, using current y and lambda values of Interactions.
   */
  void updateIndexSets();

  /** \fn virtual void updateIndexSet(const unsigned int i) = 0;
   *  \brief update indexSets[i] of the topology, using current y and lambda values of Interactions.
   *  \param unsigned int: the number of the set to be updated
   */
  virtual void updateIndexSet(const unsigned int) = 0;

  /** \fn void initialize()
   *  \brief executes the complete initialisation of Simulation (OneStepIntegrators, OneStepNSProblem, TImediscretisation) with the XML Object
   */
  virtual void initialize() = 0;

  /** \fn void computeFreeState()
   *  \brief integrates all the DynamicalSystem taking not into account nslaw, reactions ...
   */
  virtual void computeFreeState();

  /** \fn void nextStep()
   *  \brief increments all the Integrators to next step of the simulation
   */
  virtual void nextStep() = 0;

  /** \fn void computeOneStepNSProblem(const string = DEFAULT_OSNS_NAME)
   *  \brief computes the one step NS problem
   *  \param a string, the id of the OneStepNSProblem to be computed
   */
  virtual void computeOneStepNSProblem(const std::string = DEFAULT_OSNS_NAME);

  /** \fn void update(const unsigned int)
   *  \brief update input, state of each dynamical system and output
   *  \param lambda order used to compute input
   */
  virtual void update(const unsigned int) = 0;

  /** \fn void run()
   *  \brief run the simulation, from t0 to T
   */
  virtual void run() = 0;

  /** \fn void computeOneStep()
   *  \brief run one step of the simulation
   */
  virtual void computeOneStep() = 0;

  /** \fn void newtonSolve(const double criterion , const unsigned int maxStep)
   *  \brief newton algorithm
   *  \param double criterion: convergence criterion, int maxStep: maximum number of Newton steps
   */
  void newtonSolve(const double, const unsigned int);

  /** \fn newtonCheckConvergence(const double criterion);
   *  \brief check the convergence of Newton algorithm
   */
  bool Simulation::newtonCheckConvergence(const double criterion);

  // --- XML RELATED FUNCTIONS ---

  /** \fn void saveSimulationToXML()
   *  \brief copys the data of the Simulation to the XML tree
   *  \exception RuntimeException
   */
  virtual void saveSimulationToXML();

  /** \fn void updateInput(int level = levelMin)
   *  \brief compute r thanks to lambda[level]
   *  \param unsigned int: lambda order, default = levelMin
   */
  void updateInput(int = -1);

  /** \fn void updateOutput(const  int = 0, int = levelMax)
   *  \brief compute output for all the interactions
   *  \param  int: y min order to be computed, default = 0
   *  \param  int: y max order to be computed, default = levelMax
   */
  void updateOutput(const  int = 0, int = -1);

};

#endif // SIMULATION_H
