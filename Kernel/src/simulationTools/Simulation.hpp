/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
/*! \file Simulation.h
  \brief Global interface for simulation process description.
*/

#ifndef SIMULATION_H
#define SIMULATION_H

#include "SiconosConst.hpp"
#include "Tools.hpp"
#include "SimulationTypeDef.hpp"
#include "TimeDiscretisation.hpp"
#include "UnitaryRelationsSet.hpp"
#include "EventsManager.hpp"
#include "SiconosPointers.hpp"
#include "InteractionsSet.hpp"
#include "DynamicalSystemsSet.hpp"
#include <fstream>
#include <boost/bind.hpp>
#include "Model.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "Topology.hpp"
#include <boost/weak_ptr.hpp>

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
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(Simulation);


  /** name or id of the Simulation */
  std::string _name;

  /** the default time discretisation scheme */
  SP::TimeDiscretisation _timeDiscretisation;

  /** tool to manage all events */
  SP::EventsManager _eventsManager;

  /** current starting time for integration */
  double _tinit;

  /** current ending time for integration */
  double _tend;

  /** real ending time for integration (different from tend in case of
      stop during integrate, for example when a root is found in
      Lsodar procedur)
  */
  double _tout;

  /** the dynamical systems integrators */
  SP::OSISet _allOSI;

  /** Map to link all DynamicalSystems and their OneStepIntegrator*/
  DSOSIMap _osiMap;

  /** index sets vector (indexSets[0] is the set where y[0]=0,
      indexSets[1] where y[0] = 0 and y[1]=0 and so on */
  //std::vector<UnitaryRelationsGraph> indexSets;

  /** the non smooth problems (each problem is identified thanks to
      its id) */
  SP::OneStepNSProblems _allNSProblems;

  /** the XML object linked to the Simulation to read XML data */
  SP::SimulationXML _simulationxml;

  /** A link to the Model which contains the Simulation */
  boost::weak_ptr<Model> _model;

  /** _levelMinForOutput is the minimum level for the output
   * (Interaction::_lowerlevelForOutput) for all the interactions
   */
  unsigned int _levelMinForOutput;

  /** _levelMaxForOutput is the maximunm level for the output
   * (Interaction::_upperlevelForOutput) for all the interactions
   */
  unsigned int _levelMaxForOutput;

  /** _levelMinForInput is the minimum level for the input
   * (Interaction::_lowerlevelForInput) for all the interactions
   */
  unsigned int _levelMinForInput;

  /** _levelMaxForInput is the maximum level for the input
   * (Interaction::_upperlevelForInput) for all the interactions
   */
  unsigned int _levelMaxForInput;

  /** tolerance value used to compute the index sets - Default: equal
      to machine double precision (from dlamch lapack routine).*/
  double _tolerance;

  /** Flag for optional output. True if output display for solver stat
      required, else false.*/
  bool _printStat;

  /**Output file for stats*/
  std::ofstream statOut;

  /**
   * bool, option specifying if a critere of relative convergence is
   * used. Default value is false.
   */
  bool _useRelativeConvergenceCriterion;
  /**
   * bool used to remind if the relative convergence held(useful for
   * the newton-check-convergence). Modified only if
   * _useRelativeConvergenceCriterion is true.
   */
  bool _relativeConvergenceCriterionHeld;
  /**
   *double, relative tolerance. Used only if
   *_useRelativeConvergenceCriterion is true.
   */
  double _relativeConvergenceTol;

  /** initializations of levels
   *
   */
  struct SetupLevels;
  friend class Simulation::SetupLevels;

  /** initialisation for OneStepNSProblem.
   */
  virtual void initOSNS() = 0;

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
   */
  Simulation(SP::TimeDiscretisation);

  /** constructor with XML object of the Simulation
      \param SimulationXML* : the XML object corresponding
      \param initial time
      \param final time
      \param the set of all DS in the NSDS
      \param the set of all interactions in the NSDS
  */
  Simulation(SP::SimulationXML, double, double, SP::DynamicalSystemsSet,
             SP::InteractionsSet);

  /** destructor
   */
  virtual ~Simulation();

  // GETTERS/SETTERS

  /* type name of the instance */
  virtual std::string typeName() PURE_DEF;

  /** get the name of the Simulation
   *  \return string : the name of the Simulation
   */
  inline const std::string name() const
  {
    return _name;
  }

  /** set the name of the Simulation
   */
  inline void setName(const std::string& newName)
  {
    _name = newName;
  }

  /** get the TimeDiscretisation of the Simulation
   *  \return the TimeDiscretisation
   */
  inline SP::TimeDiscretisation timeDiscretisation()
  {
    return _timeDiscretisation;
  };

  /** set the TimeDiscretisation of the Simulation
   *  \param the TimeDiscretisation
   */
  inline void setTimeDiscretisationPtr(SP::TimeDiscretisation td)
  {
    _timeDiscretisation = td;
  }

  /** get time instant k of the time discretisation
   *  \return a double.
   */
  inline double getTk() const
  {
    return _timeDiscretisation->currentTime();
  };

  /** get time instant k+1 of the time discretisation - Warning: this
      instant may be different from nextTime(), if for example some
      non-smooth events or some sensor events are present
      \return a double.
  */
  inline double getTkp1() const
  {
    return _timeDiscretisation->nextTime();
  };

  /** get the EventsManager
   *  \return a pointer to EventsManager
   */
  inline SP::EventsManager eventsManager() const
  {
    return _eventsManager;
  };

  /** get "current time" (ie starting point for current integration,
      time of currentEvent of eventsManager.)
   *  \return a double.
   */
  inline double startingTime() const
  {
    return _eventsManager->startingTime();
  };

  /** get "next time" (ie ending point for current integration, time
      of nextEvent of eventsManager.)
   *  \return a double.
   */
  inline double nextTime() const
  {
    return _eventsManager->nextTime();
  };

  /** get the current time step size ("next time"-"current time")
   *  \return a double.
   */
  inline double timeStep() const
  {
    return (nextTime() - startingTime());
  };

  /** check if a future event is to be treated or not (ie if some
      events remain in the eventsManager).
   *  \return a bool.
   */
  inline bool hasNextEvent() const
  {
    return _eventsManager->hasNextEvent();
  };

  /** get all the Integrators of the Simulation
   *  \return an OSISset
   */
  inline const SP::OSISet oneStepIntegrators() const
  {
    return _allOSI;
  };

  /** set the Integrators of the Simulation
   *  \param an OSISset
   */
  void setOneStepIntegrators(const OSISet&);

  /** searchs the integrator of the DS number "numberDS"
   * \param an int ("numberDS")
   */
  SP::OneStepIntegrator integratorOfDS(int) const ;

  /** get the integrator of "ds"
   * \param a pointer to DynamicalSystem ("ds")
   */
  SP::OneStepIntegrator integratorOfDS(SP::DynamicalSystem) const ;

  /** get the number of OSIs in the Simulation (ie the size of allOSI)
   *  \return an unsigned int
   */
  inline unsigned int numberOfOSI() const
  {
    return _allOSI->size();
  }

  /** insert an Integrator into the simulation list of integrators
   *  \param a smart pointer to a OneStepIntegrator
   */
  void insertIntegrator(SP::OneStepIntegrator);

  /** register a DS and its OSI into the osiMap.
      \param a pointer to a DynamicalSystem.
   *  \param a pointer to a OneStepIntegrator.
   */
  void addInOSIMap(SP::DynamicalSystem, SP::OneStepIntegrator);


  /** get a pointer to indexSets[i]
   *  \return a UnitaryRelationsSet
   */
  SP::UnitaryRelationsGraph indexSet(unsigned int i)
  {
    return (_model.lock()->nonSmoothDynamicalSystem()->topology()->indexSet(i)) ;
  };


  /** get allNSProblems
   *  \return a pointer to OneStepNSProblems object (container of
   *  SP::OneStepNSProblem)
   */
  inline const SP::OneStepNSProblems oneStepNSProblems() const
  {
    return _allNSProblems;
  };

  /** get allNSProblems[name], a specific OneStepNSProblem
   *  \param a string, the name of the osns
   *  \return a pointer to OneStepNSProblem
   */
  SP::OneStepNSProblem oneStepNSProblem(int);

  //   /** set allNSProblems map - Warning: no copy between
  //       OneStepNSProblem of each map, pointers links!
  //    *  \param a OneStepNSProblems object (map of SP::OneStepNSProblem)
  //    */
  //   void setOneStepNSProblems(const OneStepNSProblems&);

  //   /** remove all OneStepNSProblem of the Simulation
  //    */
  //   void clearOneStepNSProblems();

  /** check if a OneStepNSProblem osns is already in the map
   *  \param a pointer to OneStepNSProblem
   *  \return a bool
   */
  //bool hasOneStepNSProblem(SP::OneStepNSProblem) const ;

  /** check if a OneStepNSProblem named id is already in the map
   *  \param a string ("id")
   *  \return a bool
   */
  //bool hasOneStepNSProblem(const std::string&) const ;

  /** add a OneStepNSProblem in the Simulation (if its not the first,
      it needs to have an id clearly defined)
   *  \param a smart pointer to OneStepNSProblem
   */
  virtual void insertNonSmoothProblem(SP::OneStepNSProblem, int Id = SICONOS_OSNSP_DEFAULT);

  /** get the SimulationXML* of the Simulation
   *  \return a pointer on the SimulationXML of the Simulation
   */
  inline SP::SimulationXML simulationXML() const
  {
    return _simulationxml;
  }

  /** set the SimulationXML of the Simulation
   *  \param SimulationXML* : the pointer to set the SimulationXML
   */
  inline void setSimulationXMLPtr(SP::SimulationXML strxml)
  {
    _simulationxml = strxml;
  }

  /** get the Model which contains the Simulation
   *  \return SP::Model : the Model which the Simulation
   */
  inline SP::Model model() const
  {
    return _model.lock();
  }

  /** set the Model which contains the Simulation
   *  \param SP::Model : the Model to set
   */
  inline void setModelPtr(SP::Model m)
  {
    _model = boost::weak_ptr<Model>(m);
  }

  /** get tolerance
   *  \return a double
   */
  double tolerance() const
  {
    return _tolerance;
  };

  /** set the value of offset for q dof vector in dynamical systems
      (to avoid events accumulation)
   *  \param a double;
   */
  void setTolerance(double inputVal)
  {
    _tolerance = inputVal;
  };

  /** set printStat value: if true, print solver stats.
   * \param a bool
   */
  inline void setPrintStat(const bool& newVal)
  {
    _printStat = newVal;
  };

  /** get printStat value
   */
  inline bool getPrintStat() const
  {
    return _printStat;
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
  int computeOneStepNSProblem(int);

  /** update input, state of each dynamical system and output
   *  \param lambda order used to compute input
   */
  virtual void update(unsigned int) = 0;

  /** update input, state of each dynamical system and output for all
      levels between levelMin and levelMax
   */
  virtual void update();

  /** run the simulation, from t0 to T
   * with default parameters if any particular settings has been done
   */
  virtual void run();

  /** step from current event to next event of EventsManager
   */
  virtual void advanceToEvent() = 0;


  /**
   * set the option to specify if a relative convergence citeron must
   * be used to stop the Newton iterations.
   *
   */
  inline void setUseRelativeConvergenceCriteron(bool use)
  {
    _useRelativeConvergenceCriterion = use;
  };
  /**
   * return true iff the relative convergence criteron is activated.
   *
   */
  inline bool useRelativeConvergenceCriteron()
  {
    return _useRelativeConvergenceCriterion;
  };

  /**
   * set the relative convergence tolerence
   *
   */
  inline void setRelativeConvergenceTol(double v)
  {
    _relativeConvergenceTol = v;
  };

  /**
   * get the relative convergence tolerence.
   *
   */
  inline const double & relativeConvergenceTol()
  {
    return _relativeConvergenceTol;
  };

  /**
   * To specify if the relative convergence criteron held.
   *
   */
  inline void setRelativeConvergenceCriterionHeld(bool newVal)
  {
    _relativeConvergenceCriterionHeld = newVal;
  };
  /**
   * return true iff the relative convergence criteron held.
   *
   */
  inline bool relativeConvergenceCriterionHeld()
  {
    return _relativeConvergenceCriterionHeld;
  };


  // --- XML RELATED FUNCTIONS ---

  /** copys the data of the Simulation to the XML tree
   *  \exception RuntimeException
   */
  virtual void saveSimulationToXML();

  /** compute r thanks to lambda[level] for all Interactions
   *  \param unsigned int: lambda level
   */
  void updateInput(unsigned int);

  /** compute output for all the interactions
   *  \param  unsigned int: y min order to be computed default 0
   */
  void updateOutput(unsigned int = 0);

  /** call eventsManager processEvents.
   */
  void processEvents();

  /**
   */
  void ComputeLevelsForInputAndOutput();

  /** visitors hook
   */
  VIRTUAL_ACCEPT_VISITORS(Simulation);

};

#endif // SIMULATION_H
