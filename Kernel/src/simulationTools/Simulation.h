/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
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
/*! \file Simulation.h
Global interface for simulation process description (Base for TimeStepping or EventDriven).
 */

#ifndef SIMULATION_H
#define SIMULATION_H

#include "SiconosConst.h"
#include "Tools.h"
#include "RuntimeException.h"
#include "UnitaryRelationsSet.h"
#include "EventsManager.h"
#include <vector>
#include <fstream>

class Model;
class DynamicalSystem;
class EventsManager;
class OneStepIntegrator;
class OneStepNSProblem;
class TimeDiscretisation;
class SimulationXML;

/** vector of OneStepIntegrator */
typedef std::set<OneStepIntegrator*> OSISet;

/** iterator through vector of OSI*/
typedef OSISet::iterator OSIIterator;

/** const iterator through vector of OSI*/
typedef OSISet::const_iterator ConstOSIIterator;

/** return type value for insert function - bool = false if insertion failed. */
typedef std::pair<OSISet::iterator, bool> CheckInsertOSI;

/** vector that contains a sequel of sets of UnitaryRelations*/
typedef std::vector< UnitaryRelationsSet* > VectorOfSetOfUnitaryRelations;

/** map of OSNS */
typedef std::map<std::string, OneStepNSProblem * > OneStepNSProblems;

/** iterator through OneStepNSProblems */
typedef OneStepNSProblems::iterator OSNSIterator;

/** const iterator through OneStepNSProblems */
typedef OneStepNSProblems::const_iterator ConstOSNSIterator;

/** A map that links DynamicalSystems and their OneStepIntegrator. */
typedef std::map<DynamicalSystem*, OneStepIntegrator*> DSOSIMap;

/** Iterator through a DSOSIMap. */
typedef DSOSIMap::iterator DSOSIIterator;

/** Const Iterator through a DSOSIMap. */
typedef DSOSIMap::const_iterator DSOSIConstIterator;

/** default tolerance value, used to update index sets */
const double DEFAULT_TOLERANCE = 10 * MACHINE_PREC;

/** Description of the simulation process (integrators, time discretisation,...) - Base class for TimeStepping or EventDriven.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.1.
 *  \date (Creation) Apr 26, 2004
 *
 * !!! This is an abstract class !!!
 *
 * The available simulations are TimeStepping and EventDriven. See derived classes for more details.
 *
 * "Unitary block matrices" are saved in the simulation (using stl maps). They will be used to save block matrices corresponding
 * to diagonal or extra diagonal terms in OneStepNSProblem "complete" matrices, such as M in LCP formulation.
 * Each block corresponds to one UnitaryRelation (diagonal terms) or two UnitaryRelations (extra diagonal terms);
 * Their dimensions are NSLawSize X NSLawSize, NSLawSize being the dimension of the NonSmoothLaw corresponding to the concerned UnitaryRelation.
 *
 * Rules:
 *     - A Model* must be given to the constructor, else => exception.
 */
class Simulation
{
protected:

  /** name or id of the Simulation */
  std::string name;

  /** the type of the Simulation, ie the derived class type (TimeStepping, EventDriven ...) */
  std::string simulationType;

  /** the time discretisation scheme */
  TimeDiscretisation *timeDiscretisation;

  /** tool to manage all events */
  EventsManager* eventsManager;

  /** Inside-class allocation flags */
  BoolMap isAllocatedIn;

  /** the dynamical systems integrators */
  OSISet allOSI;

  /** Map to link all DynamicalSystems and their OneStepIntegrator*/
  DSOSIMap osiMap;

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

  /** tolerance value used to compute the index sets - Default: equal to machine double precision (from dlamch lapack routine).*/
  double tolerance;

  /** Flag for optional output. True if output display for solver stat required, else false.*/
  bool printStat;

  /**Output file for stats*/
  std::ofstream statOut;

  /** default constructor
   */
  Simulation(const std::string& = "undefined");

  /** copy constructor. Private => no copy nor pass-by value.
   */
  Simulation(const Simulation&);

  /** initialisation for OneStepNSProblem.
   */
  virtual void initOSNS() = 0;

  /** compute LevelMax */
  virtual void initLevelMax() = 0;

public:

  /** default constructor
   *  \param a pointer to a timeDiscretisation (linked to the model that owns this simulation)
   *  \param string: simulation type, default = undefined
   */
  Simulation(TimeDiscretisation*, const std::string& = "undefined");

  /** constructor with XML object of the Simulation
   *  \param SimulationXML* : the XML object corresponding
   *  \param the model which owns this simulation  (optional parameter)
   *  \param string: simulation type, default = undefined
   */
  Simulation(SimulationXML*, Model*, const std::string& = "undefined");

  /** destructor
   */
  virtual ~Simulation();

  // GETTERS/SETTERS

  /** get the name of the Simulation
   *  \return string : the name of the Simulation
   */
  inline const std::string getName() const
  {
    return name;
  }

  /** set the name of the Simulation
   */
  inline void setName(const std::string& newName)
  {
    name = newName;
  }

  /** get the type of the Simulation
   *  \return string : the type of the Simulation
   */
  inline const std::string getType() const
  {
    return simulationType;
  }

  /** set the type of the Simulation
   */
  inline void setType(const std::string& newType)
  {
    simulationType = newType;
  }

  /** get the TimeDiscretisation of the Simulation
   *  \return the TimeDiscretisation
   */
  inline TimeDiscretisation* getTimeDiscretisationPtr() const
  {
    return timeDiscretisation;
  };

  /** set the TimeDiscretisation of the Simulation
   *  \param the TimeDiscretisation
   */
  void setTimeDiscretisationPtr(TimeDiscretisation*);

  /** get the EventsManager
   *  \return a pointer to EventsManager
   */
  inline EventsManager* getEventsManagerPtr() const
  {
    return eventsManager;
  };

  /** get "current time" (ie starting point for current integration, time of currentEvent of eventsManager.)
   *  \return a double.
   */
  inline const double getStartingTime() const
  {
    return eventsManager->getStartingTime();
  };

  /** get "next time" (ie ending point for current integration, time of nextEvent of eventsManager.)
   *  \return a double.
   */
  inline const double getNextTime() const
  {
    return eventsManager->getNextTime();
  };

  /** get the current time step size ("next time"-"current time")
   *  \return a double.
   */
  inline const double getTimeStep() const
  {
    return (getNextTime() - getStartingTime());
  };

  /** check if a future event is to be treated or not (ie if some events remain in the eventsManager).
   *  \return a bool.
   */
  inline const bool hasNextEvent() const
  {
    return eventsManager->hasNextEvent();
  };

  /** get all the Integrators of the Simulation
   *  \return an OSISset
   */
  inline const OSISet getOneStepIntegrators() const
  {
    return allOSI;
  };

  /** set the Integrators of the Simulation
   *  \param an OSISset
   */
  void setOneStepIntegrators(const OSISet&);

  /** searchs the integrator of the DS number "numberDS"
   * \param an int ("numberDS")
   */
  OneStepIntegrator* getIntegratorOfDSPtr(int) const ;

  /** get the integrator of "ds"
   * \param a pointer to DynamicalSystem ("ds")
   */
  OneStepIntegrator* getIntegratorOfDSPtr(DynamicalSystem *) const ;

  /** get the number of OSIs in the Simulation (ie the size of allOSI)
   *  \return an unsigned int
   */
  inline const unsigned int getNumberOfOSI() const
  {
    return allOSI.size();
  }

  /** add an Integrator into allOSI (pointer link, no copy!)
   *  \param a pointer to a OneStepIntegrator
   */
  void addOneStepIntegratorPtr(OneStepIntegrator *);

  /** register a DS and its OSI into the osiMap.
      \param a pointer to a DynamicalSystem.
   *  \param a pointer to a OneStepIntegrator.
   */
  void addInOSIMap(DynamicalSystem*, OneStepIntegrator *);

  /** get a pointer to indexSets[i]
   *  \return a UnitaryRelationsSet
   */
  UnitaryRelationsSet * getIndexSetPtr(unsigned int) ;

  /** get allNSProblems
   *  \return a OneStepNSProblems object (container of OneStepNSProblem*)
   */
  inline const OneStepNSProblems getOneStepNSProblems() const
  {
    return allNSProblems;
  };

  /** get allNSProblems[name], a specific OneStepNSProblem
   *  \param a string, the name of the osns
   *  \return a pointer to OneStepNSProblem
   */
  OneStepNSProblem* getOneStepNSProblemPtr(const std::string&);

  /** set allNSProblems map - Warning: no copy between OneStepNSProblem of each map, pointers links!
   *  \param a OneStepNSProblems object (map of OneStepNSProblem*)
   */
  void setOneStepNSProblems(const OneStepNSProblems&);

  /** remove all OneStepNSProblem of the Simulation
   */
  void clearOneStepNSProblems();

  /** check if a OneStepNSProblem osns is already in the map
   *  \param a pointer to OneStepNSProblem
   *  \return a bool
   */
  const bool hasOneStepNSProblem(OneStepNSProblem*) const ;

  /** check if a OneStepNSProblem named id is already in the map
   *  \param a string ("id")
   *  \return a bool
   */
  const bool hasOneStepNSProblem(const std::string&) const ;

  /** add a OneStepNSProblem in the Simulation (if its not the first, it needs to have an id clearly defined)
   *  \param a pointer to OneStepNSProblem
   */
  virtual void addOneStepNSProblemPtr(OneStepNSProblem*);

  /** get the SimulationXML* of the Simulation
   *  \return a pointer on the SimulationXML of the Simulation
   */
  inline SimulationXML* getSimulationXMLPtr() const
  {
    return simulationxml;
  }

  /** set the SimulationXML of the Simulation
   *  \param SimulationXML* : the pointer to set the SimulationXML
   */
  inline void setSimulationXMLPtr(SimulationXML* strxml)
  {
    simulationxml = strxml;
  }

  /** get the Model which contains the Simulation
   *  \return Model* : the Model which the Simulation
   */
  inline Model* getModelPtr() const
  {
    return model;
  }

  /** set the Model which contains the Simulation
   *  \param Model* : the Model to set
   */
  inline void setModelPtr(Model* m)
  {
    model = m;
  }

  /** get tolerance
   *  \return a double
   */
  const double getTolerance() const
  {
    return tolerance;
  };

  /** set the value of offset for q dof vector in dynamical systems (to avoid events accumulation)
   *  \param a double;
   */
  void setTolerance(double inputVal)
  {
    tolerance = inputVal;
  };

  /** set printStat value: if true, print solver stats.
   * \param a bool
   */
  inline void setPrintStat(const bool& newVal)
  {
    printStat = newVal;
  };

  /** get printStat value
   */
  inline const bool getPrintStat() const
  {
    return printStat;
  };

  /** update all index sets of the topology, using current y and lambda values of Interactions.
   */
  void updateIndexSets();

  /** update indexSets[i] of the topology, using current y and lambda values of Interactions.
   *  \param unsigned int: the number of the set to be updated
   */
  virtual void updateIndexSet(unsigned int) = 0;

  /** Complete initialisation of the Simulation (OneStepIntegrators, OneStepNSProblem, TImediscretisation).
   */
  void initialize();

  /** Set OSI (DS) non-smooth part to zero.
   */
  void reset();

  /** save DynamicalSystems and Interactions states in Memories (through OSI and OSNS).
   */
  void saveInMemory();

  /** computes a one step NS problem
   *  \param a string, the id of the OneStepNSProblem to be computed
   */
  virtual void computeOneStepNSProblem(const std::string&);

  /** update input, state of each dynamical system and output
   *  \param lambda order used to compute input
   */
  virtual void update(unsigned int) = 0;

  /** run the simulation, from t0 to T
   * \param: simulation option. Used only for TimeStepping.
   * \param: Used only for TimeStepping (Newton convergence criterion).
   * \param: Used only for TimeStepping (Newton maximum number of iterations).
   */
  virtual void run(const std::string& = "linear", double = 0.005, unsigned int = 500);

  /** step from current event to next event of EventsManager
   */
  virtual void advanceToEvent() = 0;

  // --- XML RELATED FUNCTIONS ---

  /** copys the data of the Simulation to the XML tree
   *  \exception RuntimeException
   */
  virtual void saveSimulationToXML();

  /** compute r thanks to lambda[level] for all Interactions
   *  \param unsigned int: lambda order, default = levelMin
   */
  void updateInput(int = -1);

  /** compute output for all the interactions
   *  \param  int: y min order to be computed, default = 0
   *  \param  int: y max order to be computed, default = levelMax
   */
  void updateOutput(int = 0, int = -1);

  /** call eventsManager processEvents.
   */
  void processEvents();
};

#endif // SIMULATION_H
