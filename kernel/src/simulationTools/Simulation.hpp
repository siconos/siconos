/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
  Global interface for simulation process description.
*/

#ifndef SIMULATION_H
#define SIMULATION_H

#include <fstream>
#include <map>
#include <set>

#include "SiconosVector.hpp"
#include "NonSmoothDynamicalSystem.hpp"  // To get ChangeLogIter
#include "SiconosConst.hpp"              // for DEFAULT_TOLERANCE
#include "SiconosSerialization.hpp"
#include "SimulationTypes.hpp"  // SICONOS_OSNSP_DEFAULT
// #include "InteractionManager.hpp"

namespace siconos::algebra {
// class SiconosVector;
}

namespace siconos::integrators {

class OneStepIntegrator;
}

namespace siconos::nonsmooth_formulations {
class OneStepNSProblem;
}

namespace siconos::simulation {

class EventsManager;
class InteractionManager;
class TimeDiscretisation;

/**
   Description of the simulation process (integrators, time
   discretisation and so on).

   This is an abstract class.

   The available simulations are TimeStepping, EventDriven and TimeSteppingD1Minus.
*/

class Simulation : public std::enable_shared_from_this<Simulation> {
 protected:
  ACCEPT_SERIALIZATION(Simulation);

  /** name or id of the Simulation */
  std::string _name{"unamed"};

  /** tool to manage all events */
  std::shared_ptr<EventsManager> _eventsManager{nullptr};

  /** current starting time for integration */
  double _tinit{0.};

  /** current ending time for integration */
  double _tend{0.};

  /** real ending time for integration (different from tend in case of
      stop during integrate, for example when a root is found in
      an EventDriven strategy)
  */
  double _tout{0.};

  double _T{0.};

  /** the dynamical systems integrators */
  std::shared_ptr<std::set<std::shared_ptr<siconos::integrators::OneStepIntegrator>>> _allOSI{
      nullptr};

  /** the non smooth problems (each problem is identified thanks to
      its id) */
  std::shared_ptr<
      std::vector<std::shared_ptr<siconos::nonsmooth_formulations::OneStepNSProblem>>>
      _allNSProblems{nullptr};

  /** A pointer to the simulated nonsmooth dynamical system
   */
  std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> _nsds{nullptr};

  /** An interaction manager
   */
  std::shared_ptr<siconos::simulation::InteractionManager> _interman{nullptr};

  /** _numberOfIndexSets is the number of index sets that we need for
   * simulation. It corresponds for most of the simulations to levelMaxForOutput + 1.
   * Nevertheless, some simulations need more sets of indices that the number
   * of outputs that we considered.
   */
  unsigned int _numberOfIndexSets{0};

  /** tolerance value used to compute the index sets.
      Default: equal to 10 x machine double precision (std::numeric_limits<double>::epsilon)
  */
  double _tolerance{siconos::internal::DEFAULT_TOLERANCE};

  /** Output setup: if true, display solver stats */
  bool _printStat{false};

  /** _staticLevels : do not recompute levels once they have been
   * initialized */
  bool _staticLevels{false};

  /** File id for stats outputs.*/
  std::ofstream statOut;

  /**
   * bool, option specifying if a critere of relative convergence is
   * used. Default value is false.
   */
  bool _useRelativeConvergenceCriterion{false};
  /**
   * bool used to remind if the relative convergence held(useful for
   * the newton-check-convergence). Modified only if
   * _useRelativeConvergenceCriterion is true.
   */
  bool _relativeConvergenceCriterionHeld{false};
  /**
   *double, relative tolerance. Used only if
   *_useRelativeConvergenceCriterion is true.
   */
  double _relativeConvergenceTol{10e-3};

  bool _isInitialized{false};

  /** current NSDS changelog position */
  std::list<siconos::modeling::NonSmoothDynamicalSystem::Change>::const_iterator
      _nsdsChangeLogPosition;

  /** map of not-yet-initialized DS variables for each OS. This an internal and temporary map
      that must be empty after initialize. */
  std::map<std::shared_ptr<siconos::integrators::OneStepIntegrator>,
           std::list<std::shared_ptr<siconos::modeling::DynamicalSystem>>>
      _OSIDSmap;

 private:
  // Rule of five
  Simulation(const Simulation&) = delete;
  Simulation(Simulation&&) = delete;
  Simulation& operator=(const Simulation&) = delete;
  Simulation& operator=(Simulation&&) = delete;

 public:
  /** default constructor, for serialization
   */
  Simulation() = default;

  /** standard constructor
   *
   *  \param nsds current nonsmooth dynamical system
   *  \param td the timeDiscretisation for this Simulation
   */
  Simulation(std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> nsds,
             std::shared_ptr<TimeDiscretisation> td);

  /** destructor
   */
  virtual ~Simulation() noexcept;

  /** clear all maps. This function should not exist, but there is a cycle
   *  with the shared_ptr: the OneStepIntegrator and OneStepNSProblem have
   *  both a link to the Simulation, and here we have all the OneStepIntegrator
   *  and OneStepNSProblem in maps. Then the memory is never freed. The clumsy
   *  way to deal with it is to call this function from the Model destructor
   *  to free the maps and then the cycle is broken
   *  \warning do not call this yourself, it is meant to be called from
   *  the desctructor of the Model
   */
  void clear();

  // GETTERS/SETTERS

  /** get the name of the Simulation
   *
   *  \return std::string : the name of the Simulation
   */
  inline const std::string name() const { return _name; }

  /** set the name of the Simulation
      \param newName the new name
   */
  inline void setName(const std::string& newName) { _name = newName; }

  /** returns time instant k of the time discretisation  */
  double getTk() const;

  /** get time instant k+1 of the time discretisation
   *  \warning: this instant may be different from nextTime(), if for example some
   *  non-smooth events or some sensor events are present
   *
   *  \return a double. If the simulation is near the end (t_{k+1} > T), it returns NaN.
   */
  double getTkp1() const;

  /** get time instant k+2 of the time discretisation
   *  \warning: this instant may be different from nextTime(), if for example some
   *  non-smooth events or some sensor events are present
   *
   *  \return a double. If the simulation is near the end (t_{k+2} > T), it returns NaN.
   */
  double getTkp2() const;

  /** returns current timestep  */
  double currentTimeStep() const;

  /** returns a pointer to the EventsManager
   */
  inline std::shared_ptr<EventsManager> eventsManager() const { return _eventsManager; };

  /** get "current time" (ie starting point for current integration,
      time of currentEvent of eventsManager.)

      \return a double.
  */
  double startingTime() const;

  /** get "next time" (ie ending point for current integration, time
      of nextEvent of eventsManager.)

      \return a double.
  */
  double nextTime() const;

  /** get the current time step size ("next time"-"current time")
   *
   *  \return a double.
   */
  inline double timeStep() const { return (nextTime() - startingTime()); };

  /** true if a future event is to be treated or not (ie if some
      events remain in the eventsManager).
  */
  bool hasNextEvent() const;

  /** get all the Integrators of the Simulation
   *
   *  \return a container of OneStepIntegrator
   */
  inline const std::shared_ptr<
      std::set<std::shared_ptr<siconos::integrators::OneStepIntegrator>>>
  oneStepIntegrators() const {
    return _allOSI;
  };

  /** get the number of OSIs in the Simulation (ie the size of allOSI)
   *
   *  \return an unsigned int
   */
  inline size_t numberOfOSI() const { return _allOSI->size(); }

  /** insert an Integrator into the simulation list of integrators
   *
   *  \param osi the OneStepIntegrator to add
   */
  virtual void insertIntegrator(std::shared_ptr<siconos::integrators::OneStepIntegrator> osi);

  /** associate an OSI with a DS */
  void associate(std::shared_ptr<siconos::integrators::OneStepIntegrator> osi,
                 std::shared_ptr<siconos::modeling::DynamicalSystem> ds);

  /**
      get a pointer to indexSets[i]

      \param i number of the required index set
      \return a graph of interactions
   */
  std::shared_ptr<siconos::graphs::InteractionsGraph> indexSet(unsigned int i);

  /** get allNSProblems
   *
   *  \return a set of pointers to OneStepNSProblem
   */
  inline auto oneStepNSProblems() const { return _allNSProblems; };

  /** get the number of OSNSP in the Simulation (ie the size of allNSProblems)
   *
   *  \return an unsigned int
   */
  inline size_t numberOfOSNSProblems() const { return _allNSProblems->size(); }

  /**
     get a OneStep nonsmooth problem of the simulation, identify with its number.

     \param id number of the required osnspb
     \return a pointer to OneStepNSProblem
  */
  std::shared_ptr<siconos::nonsmooth_formulations::OneStepNSProblem> oneStepNSProblem(int id);

  /**
      add a OneStepNSProblem in the Simulation

      \param osns the OneStepNSProblem to insert
      \param Id its id: default is SICONOS_OSNSP_DEFAULT,
      at impact level SICONOS_OSNSP_ED_IMPACT, at acceleration level
      SICONOS_OSNSP_ED_ACCELERATION
   */
  virtual void insertNonSmoothProblem(
      std::shared_ptr<siconos::nonsmooth_formulations::OneStepNSProblem> osns,
      int Id = SICONOS_OSNSP_DEFAULT);

  /** get the NonSmoothDynamicalSystem
   *
   *  \return NonSmoothDynamicalSystem
   */
  inline std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem>
  nonSmoothDynamicalSystem() const {
    return _nsds;
  }
  /** set the NonSmoothDynamicalSystem of the Simulation
   *
   *  \param newPtr a pointer on NonSmoothDynamicalSystem
   */
  void setNonSmoothDynamicalSystemPtr(
      std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> newPtr) {
    _nsdsChangeLogPosition = _nsds->changeLog().begin();
    _nsds = newPtr;
  }

  /** get tolerance
   *
   *  \return a double
   */
  double tolerance() const { return _tolerance; };

  /** set the value of offset for q dof vector in dynamical systems
      (to avoid events accumulation)

      \param inputVal new tolerance
   */
  void setTolerance(double inputVal) { _tolerance = inputVal; };

  /**
     set printStat value: if true, print solver stats.

     \param newVal true to activate stats
  */
  inline void setPrintStat(const bool& newVal) { _printStat = newVal; };

  /** \return true if stats are activated
   */
  inline bool getPrintStat() const { return _printStat; };

  virtual void computeInitialStateOfTheStep(){};

  /**
      update all index sets of the topology, using current y and lambda values of Interactions
   */
  void updateIndexSets();

  /** update indexSets[i] of the topology, using current y and lambda
      values of Interactions.

      \param level the number of the set to be updated
   */
  virtual void updateIndexSet(unsigned int level) = 0;

  /** Complete initialisation of the Simulation (OneStepIntegrators,
      siconos::nonsmooth_formulations::OneStepNSProblem, TImediscretisation).
  */
  virtual void initialize() = 0;

  /** Initialize a single Interaction for this Simulation, used for dynamic
   *  topology updates. */
  virtual void initializeInteraction(double time,
                                     std::shared_ptr<siconos::modeling::Interaction> inter);

  /** Set an object to automatically manage interactions during the simulation
   *
   *  \param manager
   */
  void insertInteractionManager(std::shared_ptr<InteractionManager> manager) {
    _interman = manager;
  }

  /** Compute the residu of all OSI
   */
  void computeResidu();

  /** computes a one step NS problem
   *
   *  \param nb the id of the siconos::nonsmooth_formulations::OneStepNSProblem to be computed
   *  \return information about the solver convergence.
   */
  int computeOneStepNSProblem(int nb);

  /** update the plugins of the DS
   *
   *  \param time to be used for plugins
   */
  virtual void updateDSPlugins(double time);

  /** update input
   *  \param level lambda order used to compute input
   */
  virtual void updateInput(unsigned int level);

  /** update all input terms
   */
  virtual void updateAllInput();

  /** update state of each dynamical system
   */
  virtual void updateState(unsigned int level = 0);

  /** update output
   *
   *  \param level lambda order used to compute output
   *  level is set to 0 by default since in all time-stepping schemes we update all the state
   */
  virtual void updateOutput(unsigned int level = 0);

  /** update output, state, and input
   *
   *  \param level lambda order used to compute input
   *  level is set to 0 by default since in all time-stepping schemes we update all the state
   */
  void update(unsigned int level = 0) {
    updateInput(level);
    updateState(level);
    updateOutput(level);
  }

  /** run the simulation, from t0 to T
   *  with default parameters if any particular settings has been done
   */
  virtual void run();

  /** initialisation for siconos::nonsmooth_formulations::OneStepNSProblem.
   */
  virtual void initOSNS() = 0;

  /** step from current event to next event of EventsManager
   */
  virtual void advanceToEvent() = 0;

  /** clear the NSDS changelog up to current position.  If you have a
   *  particularly dynamic simulation (DS and Interactions created and
   *  destroyed frequently), then it is important to call this
   *  periodically.
   */
  void clearNSDSChangeLog();

  /** Set the option to specify if a relative convergence criterion must
      be used to stop the Newton iterations.

      \param use true if relative critarion activated
  */
  inline void setUseRelativeConvergenceCriteron(bool use) {
    _useRelativeConvergenceCriterion = use;
  };
  /**\return true if the relative convergence criterion is activated.
   */
  inline bool useRelativeConvergenceCriteron() { return _useRelativeConvergenceCriterion; };

  /**
     Set the relative convergence tolerance

     \param v tolerance value
  */
  inline void setRelativeConvergenceTol(double v) { _relativeConvergenceTol = v; };

  /** \return the relative convergence tolerence.
   *
   */
  inline double relativeConvergenceTol() { return _relativeConvergenceTol; };

  /**
   \param newVal a new relative convergence criterion
   *
   */
  inline void setRelativeConvergenceCriterionHeld(bool newVal) {
    _relativeConvergenceCriterionHeld = newVal;
  };
  /**
     \return true if the relative convergence criterion held.
   *
   */
  inline bool relativeConvergenceCriterionHeld() { return _relativeConvergenceCriterionHeld; };

  /**
     return input lambda[level](coor) for all the interactions

     \param level lambda min order to be computed
     \param coor the coordinate of interest
     \return a std::shared_ptr<siconos::algebra::SiconosVector> that contains the concatenated
     value
  */
  std::shared_ptr<siconos::algebra::SiconosVector> lambda(unsigned int level = 0,
                                                          unsigned int coor = 0);

  /**
      return output y[level](coor) for all the interactions

      \param level y min order to be computed
      \param coor the coordinate of interest
      \return a std::shared_ptr<siconos::algebra::SiconosVector> that contains the concatenated
     value
  */
  std::shared_ptr<siconos::algebra::SiconosVector> y(unsigned int level = 0,
                                                     unsigned int coor = 0);

  /** call eventsManager processEvents.
   */
  void processEvents();

  /** set staticLevels
   *
   *  \param b decides whether levels should be computed at each iteration
   */
  void setStaticLevels(bool b) { _staticLevels = b; }

  /** This updates the end of the Simulation.
   *  \warning this should be called only from the Model, to synchronise the 2 values
   *
   *  \param T the new final time
   */
  void updateT(double T);

  virtual bool computeResiduY() { return false; };

  virtual bool computeResiduR() { return false; };

  /** Add a new Interaction between one or a pair of DSs.
   *
   *  \param inter the std::shared_ptr<siconos::modeling::Interaction> to add
   *  \param ds1 the first std::shared_ptr<siconos::modeling::DynamicalSystem> in the
   * Interaction \param ds2 the second std::shared_ptr<siconos::modeling::DynamicalSystem> in
   * the Interaction, if any
   */
  void link(std::shared_ptr<siconos::modeling::Interaction> inter,
            std::shared_ptr<siconos::modeling::DynamicalSystem> ds1,
            std::shared_ptr<siconos::modeling::DynamicalSystem> ds2 =
                std::shared_ptr<siconos::modeling::DynamicalSystem>());

  /** Remove an Interaction from the simulation.
   *
   *  \param inter the std::shared_ptr<siconos::modeling::Interaction> to remove
   */
  void unlink(std::shared_ptr<siconos::modeling::Interaction> inter);

  /** Call the interaction manager if one is registered, otherwise do nothing. */
  void updateInteractions();

  /*TS set the ds->q memory, the world (CAD model for example) must be updated.
    Overload this method to update user model.*/
  virtual void updateWorldFromDS() { ; };

  /** initialize OSI-DS links in the NSDS graph. */
  void initializeOSIAssociations();

  /** initialize objects (DSs and Interations) found in the NSDS
   *  Changelog and update the changelog iterator.
   */
  void initializeNSDSChangelog();

  void applyNSDSChangelogForDS();

  /** initialize index sets for OSIs */
  void initializeIndexSets();

  /** Complete initialisation of the Simulation (OneStepIntegrators,
      siconos::nonsmooth_formulations::OneStepNSProblem, TImediscretisation).
  */
  virtual void firstInitialize();
};
}  // namespace siconos::simulation
#endif  // SIMULATION_H
