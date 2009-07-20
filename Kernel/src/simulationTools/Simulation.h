/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
  \brief Global interface for simulation process description.
*/

#ifndef SIMULATION_H
#define SIMULATION_H

#include "SiconosConst.h"
#include "Tools.hpp"
#include "SimulationTypeDef.hpp"
#include "TimeDiscretisation.h"
#include "UnitaryRelationsSet.hpp"
#include "EventsManager.h"
#include "SiconosPointers.hpp"
#include "InteractionsSet.hpp"
#include "DynamicalSystemsSet.hpp"
#include <fstream>
#include <boost/bind.hpp>
#include "Model.h"
#include "NonSmoothDynamicalSystem.h"
#include "Topology.h"

class DynamicalSystem;
class EventsManager;
class OneStepIntegrator;
class OneStepNSProblem;
class TimeDiscretisation;
class SimulationXML;

/** Description of the simulation process (integrators, time
    discretisation and so on) - Base class for TimeStepping or
    EventDriven.

    \author SICONOS Development Team - copyright INRIA
    \version 3.0.0.
    \date (Creation) Apr 26, 2004

    !!! This is an abstract class !!!

    The available simulations are TimeStepping and EventDriven. See
    derived classes for more details.



  Rules:
  - A Model must be given to the constructor, else => exception.
*/
class Simulation : public boost::enable_shared_from_this<Simulation>
{
protected:

  /** name or id of the Simulation */
  std::string name;

  /** the type of the Simulation, ie the derived class type
      (TimeStepping, EventDriven ...) */
  std::string simulationType;

  /** the default time discretisation scheme */
  SP::TimeDiscretisation timeDiscretisation;

  /** tool to manage all events */
  SP::EventsManager eventsManager;

  /** current starting time for integration */
  double tinit;

  /** current ending time for integration */
  double tend;

  /** real ending time for integration (different from tend in case of
      stop during integrate, for example when a root is found in
      Lsodar procedur)
  */
  double tout;

  /** the dynamical systems integrators */
  SP::OSISet allOSI;

  /** Map to link all DynamicalSystems and their OneStepIntegrator*/
  DSOSIMap osiMap;

  /** index sets vector (indexSets[0] is the set where y[0]=0,
      indexSets[1] where y[0] = 0 and y[1]=0 and so on */
  //std::vector<UnitaryRelationsGraph> indexSets;

  /** the non smooth problems (each problem is identified thanks to
      its id) */
  SP::OneStepNSProblems allNSProblems;

  /** the XML object linked to the Simulation to read XML data */
  SP::SimulationXML simulationxml;

  /** A link to the Model which contains the Simulation */
  SP::Model model;

  /** int used to set the minimal derivative order used in the OSNS
      variables */
  unsigned int levelMin;

  /** int used to set the maximal derivative order used in the OSNS
      variables */
  unsigned int levelMax;

  /** tolerance value used to compute the index sets - Default: equal
      to machine double precision (from dlamch lapack routine).*/
  double tolerance;

  /** Flag for optional output. True if output display for solver stat
      required, else false.*/
  bool printStat;

  /**Output file for stats*/
  std::ofstream statOut;

  /**
   * bool, option specifying if a critere of relative convergence is used. Default value is false.
   *
   *
   */
  bool mUseRelativeConvergenceCriterion;
  /**
   * bool used to remind if the relative convergence held(useful for the newton-check-convergence). Modified only if mUseRelativeConvergenceCriterion is true.
   *
   *
   */
  bool mRelativeConvergenceCriterionHeld;
  /**
   *double, relative tolerence. Used only if mUseRelativeConvergenceCriterion is true.
   *
   *
   */
  double mRelativeConvergenceTol;


  /** initialisation for OneStepNSProblem.
   */
  virtual void initOSNS() = 0;

  /** compute LevelMax */
  virtual void initLevelMax() = 0;

  /** default constructor.
   */
  Simulation() {};

private:

  /** copy constructor. Private => no copy nor pass-by value.
   */
  Simulation(const Simulation&);

  /** assignment of simulation. Private => no copy nor pass-by value.
   */
  Simulation& operator=(const Simulation&);

public:

  /** default constructor
   *  \param a pointer to a timeDiscretisation (linked to the model
   *  that owns this simulation)
   *  \param string: simulation type, default = undefined
   */
  Simulation(SP::TimeDiscretisation, const std::string& = "undefined");

  /** constructor with XML object of the Simulation
      \param SimulationXML* : the XML object corresponding
      \param initial time
      \param final time
      \param the set of all DS in the NSDS
      \param the set of all interactions in the NSDS
      \param string: simulation type, default = undefined
  */
  Simulation(SP::SimulationXML, double, double, SP::DynamicalSystemsSet,
             SP::InteractionsSet, const std::string& = "undefined");

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
  inline SP::TimeDiscretisation getTimeDiscretisationPtr()
  {
    return timeDiscretisation;
  };

  /** set the TimeDiscretisation of the Simulation
   *  \param the TimeDiscretisation
   */
  inline void setTimeDiscretisationPtr(SP::TimeDiscretisation td)
  {
    timeDiscretisation = td;
  }

  /** get time instant k of the time discretisation
   *  \return a double.
   */
  inline const double getTk() const
  {
    return timeDiscretisation->getCurrentTime();
  };

  /** get time instant k+1 of the time discretisation - Warning: this
      instant may be different from getNextTime(), if for example some
      non-smooth events or some sensor events are present
      \return a double.
  */
  inline const double getTkp1() const
  {
    return timeDiscretisation->getNextTime();
  };

  /** get the EventsManager
   *  \return a pointer to EventsManager
   */
  inline SP::EventsManager getEventsManagerPtr() const
  {
    return eventsManager;
  };

  /** get "current time" (ie starting point for current integration,
      time of currentEvent of eventsManager.)
   *  \return a double.
   */
  inline const double getStartingTime() const
  {
    return eventsManager->getStartingTime();
  };

  /** get "next time" (ie ending point for current integration, time
      of nextEvent of eventsManager.)
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

  /** check if a future event is to be treated or not (ie if some
      events remain in the eventsManager).
   *  \return a bool.
   */
  inline const bool hasNextEvent() const
  {
    return eventsManager->hasNextEvent();
  };

  /** get all the Integrators of the Simulation
   *  \return an OSISset
   */
  inline const SP::OSISet getOneStepIntegrators() const
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
  SP::OneStepIntegrator getIntegratorOfDSPtr(int) const ;

  /** get the integrator of "ds"
   * \param a pointer to DynamicalSystem ("ds")
   */
  SP::OneStepIntegrator getIntegratorOfDSPtr(SP::DynamicalSystem) const ;

  /** get the number of OSIs in the Simulation (ie the size of allOSI)
   *  \return an unsigned int
   */
  inline const unsigned int getNumberOfOSI() const
  {
    return allOSI->size();
  }

  /** add an Integrator into the simulation list of integrators
      (pointer link, no copy!)
   *  \param a smart pointer to a OneStepIntegrator
   */
  void recordIntegrator(SP::OneStepIntegrator);

  /** register a DS and its OSI into the osiMap.
      \param a pointer to a DynamicalSystem.
   *  \param a pointer to a OneStepIntegrator.
   */
  void addInOSIMap(SP::DynamicalSystem, SP::OneStepIntegrator);


  /** get a pointer to indexSets[i]
   *  \return a UnitaryRelationsSet
   */
  SP::UnitaryRelationsGraph getIndexSetPtr(unsigned int i)
  {
    return (model->getNonSmoothDynamicalSystemPtr()->getTopologyPtr()->getIndexSetPtr(i)) ;
  };


  /** get allNSProblems
   *  \return a pointer to OneStepNSProblems object (container of
   *  SP::OneStepNSProblem)
   */
  inline const SP::OneStepNSProblems getOneStepNSProblems() const
  {
    return allNSProblems;
  };

  /** get levelMin
   *  \return the value of LevelMin
   */
  inline const int getLevelMin() const
  {
    return levelMin;
  };

  /** get levelMax
   *  \return the value of LevelMax
   */
  inline const int getLevelMax() const
  {
    return levelMax;
  };

  /** get allNSProblems[name], a specific OneStepNSProblem
   *  \param a string, the name of the osns
   *  \return a pointer to OneStepNSProblem
   */
  SP::OneStepNSProblem getOneStepNSProblemPtr(const std::string&);

  /** set allNSProblems map - Warning: no copy between
      OneStepNSProblem of each map, pointers links!
   *  \param a OneStepNSProblems object (map of SP::OneStepNSProblem)
   */
  void setOneStepNSProblems(const OneStepNSProblems&);

  /** remove all OneStepNSProblem of the Simulation
   */
  void clearOneStepNSProblems();

  /** check if a OneStepNSProblem osns is already in the map
   *  \param a pointer to OneStepNSProblem
   *  \return a bool
   */
  const bool hasOneStepNSProblem(SP::OneStepNSProblem) const ;

  /** check if a OneStepNSProblem named id is already in the map
   *  \param a string ("id")
   *  \return a bool
   */
  const bool hasOneStepNSProblem(const std::string&) const ;

  /** add a OneStepNSProblem in the Simulation (if its not the first,
      it needs to have an id clearly defined)
   *  \param a smart pointer to OneStepNSProblem
   */
  virtual void recordNonSmoothProblem(SP::OneStepNSProblem);

  /** get the SimulationXML* of the Simulation
   *  \return a pointer on the SimulationXML of the Simulation
   */
  inline SP::SimulationXML getSimulationXMLPtr() const
  {
    return simulationxml;
  }

  /** set the SimulationXML of the Simulation
   *  \param SimulationXML* : the pointer to set the SimulationXML
   */
  inline void setSimulationXMLPtr(SP::SimulationXML strxml)
  {
    simulationxml = strxml;
  }

  /** get the Model which contains the Simulation
   *  \return SP::Model : the Model which the Simulation
   */
  inline SP::Model getModelPtr() const
  {
    return model;
  }

  /** set the Model which contains the Simulation
   *  \param SP::Model : the Model to set
   */
  inline void setModelPtr(SP::Model m)
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

  /** set the value of offset for q dof vector in dynamical systems
      (to avoid events accumulation)
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

  /** update all index sets of the topology, using current y and
      lambda values of Interactions.
   */
  void updateIndexSets();

  /** update indexSets[i] of the topology, using current y and lambda
      values of Interactions.
   *  \param unsigned int: the number of the set to be updated
   */
  virtual void updateIndexSet(unsigned int) = 0;

  /** Complete initialisation of the Simulation (OneStepIntegrators,
      OneStepNSProblem, TImediscretisation).
      \param the model, which will own the Simulation
      \param optional flag for partial initialisation
  */
  void initialize(SP::Model, bool = true);

  /** Update simulation if some interactions are modified */
  void updateInteractions();

  /** Set OSI (DS) non-smooth part to zero.
   */
  void reset();

  /** save DynamicalSystems and Interactions states in Memories
      (through OSI and OSNS).
   */
  void saveInMemory();

  /** computes a one step NS problem
   *  \param a string, the id of the OneStepNSProblem to be computed
   *  \return an int, information about the solver convergence.
   */
  int computeOneStepNSProblem(const std::string&);

  /** update input, state of each dynamical system and output
   *  \param lambda order used to compute input
   */
  virtual void update(unsigned int) = 0;

  /** update input, state of each dynamical system and output for all
      levels between levelMin and levelMax
   */
  virtual void update();

  /** run the simulation, from t0 to T
   * \param: simulation option. Used only for TimeStepping.
   * \param: Used only for TimeStepping (Newton convergence criterion).
   * \param: Used only for TimeStepping (Newton maximum number of iterations).
   */
  virtual void run(const std::string& = "linear", double = 0.005, unsigned int = 500);

  /** step from current event to next event of EventsManager
   */
  virtual void advanceToEvent() = 0;


  /**
   * set the option to specify if a relative convergence citeron must be used to stop the Newton iterations.
   *
   */
  inline void setUseRelativeConvergenceCriteron(bool use)
  {
    mUseRelativeConvergenceCriterion = use;
  };
  /**
   * return true iff the relative convergence criteron is activated.
   *
   */
  inline bool getUseRelativeConvergenceCriteron()
  {
    return mUseRelativeConvergenceCriterion;
  };

  /**
   * set the relative convergence tolerence
   *
   */
  inline void setRelativeConvergenceTol(double v)
  {
    mRelativeConvergenceTol = v;
  };

  /**
   * get the relative convergence tolerence.
   *
   */
  inline const double & getRelativeConvergenceTol()
  {
    return mRelativeConvergenceTol;
  };

  /**
   * To specify if the relative convergence criteron held.
   *
   */
  inline void setRelativeConvergenceCriterionHeld(bool newVal)
  {
    mRelativeConvergenceCriterionHeld = newVal;
  };
  /**
   * return true iff the relative convergence criteron held.
   *
   */
  inline bool getRelativeConvergenceCriterionHeld()
  {
    return mRelativeConvergenceCriterionHeld;
  };


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

  /** visitors hook
   */
  VIRTUAL_ACCEPT_VISITORS(Simulation);

};

#endif // SIMULATION_H
