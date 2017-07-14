/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
/*! \file Simulation.hpp
  \brief Global interface for simulation process description.
*/

#ifndef SIMULATION_H
#define SIMULATION_H

#include "SiconosConst.hpp"
#include "SimulationTypeDef.hpp"
#include "SiconosFwd.hpp"
// #include "EventsManager.hpp"
// #include "SiconosPointers.hpp"
#include <fstream>
// #include "Model.hpp"
// #include "NonSmoothDynamicalSystem.hpp"
// #include "Topology.hpp"
#include "InteractionManager.hpp"


/** Description of the simulation process (integrators, time
    discretisation and so on) - Base class for TimeStepping or
    EventDriven.

    \author SICONOS Development Team - copyright INRIA
    \version 3.8.0.
    \date (Creation) Apr 26, 2004

    !!! This is an abstract class !!!

    The available simulations are TimeStepping and EventDriven. See
    derived classes for more details.

  Rules:
  - A Model must be given to the constructor, otherwise an exception is thrown.
*/
class Simulation : public std11::enable_shared_from_this<Simulation>
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(Simulation);


  /** name or id of the Simulation */
  std::string _name;

  /** tool to manage all events */
  SP::EventsManager _eventsManager;

  /** current starting time for integration */
  double _tinit;

  /** current ending time for integration */
  double _tend;

  /** real ending time for integration (different from tend in case of
      stop during integrate, for example when a root is found in
      LsodarOSI procedure)
  */
  double _tout;

  double _T;

  /** the dynamical systems integrators */
  SP::OSISet _allOSI;

  /** the non smooth problems (each problem is identified thanks to
      its id) */
  SP::OneStepNSProblems _allNSProblems;

  /** A pointer to the simulated nonsmooth dynamical system
   */
  SP::NonSmoothDynamicalSystem _nsds;

  /** An interaction manager
   */
  SP::InteractionManager _interman;

  /** _numberOfIndexSets is the number of index sets that we need for
   * simulation. It corresponds for most of the simulation to
   * _numberOfIndexSets = _levelMaxForOutput + 1
   * Nevetheless, some simulation needs more sets of indices that the number
   * of output that we considered.
   */
  unsigned int _numberOfIndexSets;

  /** tolerance value used to compute the index sets - Default: equal
      to machine double precision (from dlamch lapack routine).*/
  double _tolerance;

  /** Flag for optional output. True if output display for solver stat
      required, else false.*/
  bool _printStat;

  /** _staticLevels : do not recompute levels once they have been
   * initialized */
  bool _staticLevels;

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

  /**
   * track whether link() or unlink() has been called
   */
  bool _linkOrUnlink;

  /** default constructor.
   */
  Simulation() {};

  /** Call the interaction manager one if is registered, otherwise do nothing. */
  void updateInteractions();

  /** Call the interaction manager one if is registered, otherwise do nothing. */
  void updateInteractionsNewtonIteration();

  /*TS set the ds->q memory, the world (CAD model for example) must be updated.
    Overload this method to update user model.*/
  virtual void updateWorldFromDS()
  {
    ;
  };

private:

  /** copy constructor. Private => no copy nor pass-by value.
   */
  Simulation(const Simulation&);

  /* assignment of simulation. Private => no copy nor pass-by value.
   */
  Simulation& operator=(const Simulation&);


public:


  /** default constructor
   *  \param td the timeDiscretisation for this Simulation
   */
  Simulation(SP::TimeDiscretisation td);

  /** destructor
   */
  virtual ~Simulation();

  /** clear all maps. This function should not exist, but there is a cycle
   * with the shared_ptr: the OneStepIntegrator and OneStepNSProblem have
   * both a link to the Simulation, and here we have all the OneStepIntegrator
   * and OneStepNSProblem in maps. Then the memory is never freed. The clumsy
   * way to deal with it is to call this function from the Model destructor
   * to free the maps and then the cycle is broken
   * \warning do not call this yourself, it is meant to be called from
   * the desctructor of the Model
   */
  void clear();

  // GETTERS/SETTERS

  /** get the name of the Simulation
   *  \return std::string : the name of the Simulation
   */
  inline const std::string name() const
  {
    return _name;
  }

  /** set the name of the Simulation
      \param newName the new name
   */
  inline void setName(const std::string& newName)
  {
    _name = newName;
  }

  /** get time instant k of the time discretisation
   *  \return the time instant t_k
   */
  double getTk() const;

  /** get time instant k+1 of the time discretisation
   * \warning: this instant may be different from nextTime(), if for example some
      non-smooth events or some sensor events are present
     \return a double. If the simulation is near the end (t_{k+1} > T), it returns NaN.
   */
  double getTkp1() const;

  /** get time instant k+2 of the time discretisation
   * \warning: this instant may be different from nextTime(), if for example some
      non-smooth events or some sensor events are present
   * \return a double. If the simulation is near the end (t_{k+2} > T), it returns NaN.
   */
  double getTkp2() const;

  /** Get current timestep
   * \return the current timestep
   */
  double currentTimeStep() const;

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
  double startingTime() const;

  /** get "next time" (ie ending point for current integration, time
      of nextEvent of eventsManager.)
   *  \return a double.
   */
  double nextTime() const;

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
  bool hasNextEvent() const;

  /** get all the Integrators of the Simulation
   *  \return an OSISset
   */
  inline const SP::OSISet oneStepIntegrators() const
  {
    return _allOSI;
  };

  /** get the number of OSIs in the Simulation (ie the size of allOSI)
   *  \return an unsigned int
   */
  inline unsigned int numberOfOSI() const
  {
    return _allOSI->size();
  }

  /** insert an Integrator into the simulation list of integrators
   *  \param osi the OneStepIntegrator to add
   */
  virtual void insertIntegrator(SP::OneStepIntegrator osi);

  /** get a pointer to indexSets[i]
      \param i number of the required index set
      \return a graph of interactions
   */
  SP::InteractionsGraph indexSet(unsigned int i);

  /** get allNSProblems
   *  \return a pointer to OneStepNSProblems object (container of
   *  SP::OneStepNSProblem)
   */
  inline const SP::OneStepNSProblems oneStepNSProblems() const
  {
    return _allNSProblems;
  };

  /** get the number of OSNSP in the Simulation (ie the size of allNSProblems)
   *  \return an unsigned int
   */
  inline unsigned int numberOfOSNSProblems() const
  {
    return _allNSProblems->size();
  }

  /** get a OneStep nonsmooth problem of the simulation, identify with its number.
      \param id number of the required osnspb
      \return a pointer to OneStepNSProblem
   */
  SP::OneStepNSProblem oneStepNSProblem(int id);

  /** add a OneStepNSProblem in the Simulation
      \param osns the OneStepNSProblem to insert
      \param Id its id: default is SICONOS_OSNSP_DEFAULT,
      at impact level SICONOS_OSNSP_ED_IMPACT, at acceleration level
      SICONOS_OSNSP_ED_ACCELERATION
   */
  virtual void insertNonSmoothProblem(SP::OneStepNSProblem osns, int Id = SICONOS_OSNSP_DEFAULT);


  /** get the NonSmoothDynamicalSystem
   *  \return NonSmoothDynamicalSystem
   */
  inline SP::NonSmoothDynamicalSystem nonSmoothDynamicalSystem() const
  {
    return _nsds;
  }

  /** set the NonSmoothDynamicalSystem of the Simulation
   *  \param newPtr a pointer on NonSmoothDynamicalSystem
   */
  void setNonSmoothDynamicalSystemPtr(SP::NonSmoothDynamicalSystem newPtr)
  {
    _nsds = newPtr;
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
      \param inputVal new tolerance
   */
  void setTolerance(double inputVal)
  {
    _tolerance = inputVal;
  };

  /** set printStat value: if true, print solver stats.
      \param newVal true to activate stats
   */
  inline void setPrintStat(const bool& newVal)
  {
    _printStat = newVal;
  };

  /** get printStat value
      \return true if stats are activated
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
      \param level the number of the set to be updated
   */
  virtual void updateIndexSet(unsigned int level) = 0;

  /** Complete initialisation of the Simulation (OneStepIntegrators,
      OneStepNSProblem, TImediscretisation).
      \param m the model to be linked to this Simulation
      \param init optional flag for partial initialization
  */
  virtual void initialize(SP::Model m, bool init = true);

  /** Initialize a single Interaction for this Simulation, used for dynamic
   *  topology updates. */
  virtual void initializeInteraction(double time, SP::Interaction inter);

  /** Associate an OSI with a DynamicalSystem in the graph and
   *  initialize any necessary graph properties (e.g. work vectors).
   *  Inserts the integrator into the set if not already present.
   *
   *  \param osi The OneStepIntegrator to associate with the DynamicalSystem.
   *  \param ds The DynamicalSystem, which must be already inserted
   *            into the NonSmoothDynamicalSystem.
   *  \param m The Model for initializing the OSI.
   *  \param time The current time for initializing the OSI. */
  void prepareIntegratorForDS(SP::OneStepIntegrator osi, SP::DynamicalSystem ds,
                              SP::Model m, double time);

  /** Set an object to automatically manage interactions during the simulation */
  void insertInteractionManager(SP::InteractionManager manager)
    { _interman = manager; }

  /** computes a one step NS problem
   *  \param nb the id of the OneStepNSProblem to be computed
   *  \return information about the solver convergence.
   */
  int computeOneStepNSProblem(int nb);

  /** update input
   *  \param level lambda order used to compute input
   * level is set to 0 by default since in all time-stepping schemes we update all the state
   */
  virtual void updateInput(unsigned int level=0);

  /** update state of each dynamical system
   */
  virtual void updateState(unsigned int level=0);

  /** update output
   *  \param level lambda order used to compute output
   * level is set to 0 by default since in all time-stepping schemes we update all the state
   */
  virtual void updateOutput(unsigned int level=0);

  /** update output, state, and input
   *  \param level lambda order used to compute input
   * level is set to 0 by default since in all time-stepping schemes we update all the state
   */
  void update(unsigned int level=0)
    { updateInput(level); updateState(level); updateOutput(level); }

  /** run the simulation, from t0 to T
   * with default parameters if any particular settings has been done
   */
  virtual void run();

  /** initialisation for OneStepNSProblem.
   */
  virtual void initOSNS() = 0;


  /** step from current event to next event of EventsManager
   */
  virtual void advanceToEvent() = 0;


  /* Set the option to specify if a relative convergence criterion must
   be used to stop the Newton iterations.
   \param use true if relative critarion activated
  */
  inline void setUseRelativeConvergenceCriteron(bool use)
  {
    _useRelativeConvergenceCriterion = use;
  };
  /**
     \return true if the relative convergence criterion is activated.
   */
  inline bool useRelativeConvergenceCriteron()
  {
    return _useRelativeConvergenceCriterion;
  };

  /** Set the relative convergence tolerance
      \param v tolerance value
   */
  inline void setRelativeConvergenceTol(double v)
  {
    _relativeConvergenceTol = v;
  };

  /**
   \return the relative convergence tolerence.
   *
   */
  inline double relativeConvergenceTol()
  {
    return _relativeConvergenceTol;
  };

  /**
   \param newVal a new relative convergence criterion
   *
   */
  inline void setRelativeConvergenceCriterionHeld(bool newVal)
  {
    _relativeConvergenceCriterionHeld = newVal;
  };
  /**
     \return true if the relative convergence criterion held.
   *
   */
  inline bool relativeConvergenceCriterionHeld()
  {
    return _relativeConvergenceCriterionHeld;
  };

  /** return input lambda[level](coor) for all the interactions
      \param level lambda min order to be computed
      \param coor the coordinate of interest
      \return a SP::SiconosVector that contains the concatenated value
   */
  SP::SiconosVector lambda(unsigned int level = 0, unsigned int coor=0 );

  /** return output y[level](coor) for all the interactions
      \param level y min order to be computed
      \param coor the coordinate of interest
      \return a SP::SiconosVector that contains the concatenated value
   */
  SP::SiconosVector y(unsigned int level = 0, unsigned int coor=0 );

  /** call eventsManager processEvents.
   */
  void processEvents();


  /** set staticLevels
   * \param b decides whether levels should be computed at each iteration
   */
  void setStaticLevels(bool b)
  {
    _staticLevels = b;
  }

  /** This updates the end of the Simulation.
   * \warning this should be called only from the Model, to synchronise the 2 values
   * \param T the new final time
   */
  void updateT(double T);

  virtual bool computeResiduY() { return false; };

  virtual bool computeResiduR() { return false; };

  /** Add a new Interaction between one or a pair of DSs.
   * \param inter the SP::Interaction to add
   * \param ds1 the first SP::DynamicalSystem in the Interaction
   * \param ds2 the second SP::DynamicalSystem in the Interaction, if any
   */
  void link(SP::Interaction inter,
            SP::DynamicalSystem ds1,
            SP::DynamicalSystem ds2 = SP::DynamicalSystem());

  /** Remove an Interaction from the simulation.
   * \param inter the SP::Interaction to remove
   */
  void unlink(SP::Interaction inter);

  /** visitors hook
   */
  VIRTUAL_ACCEPT_VISITORS(Simulation);
};

#endif // SIMULATION_H
